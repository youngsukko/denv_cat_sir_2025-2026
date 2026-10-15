
data {
  int<lower=1> S;
  int<lower=1> A5;
  int<lower=1> A1;
  int<lower=1> T;

  array[A1] int<lower=1, upper=A5> band_of_age;

  array[A5, T, S] int<lower=0> cases;
  array[A5, T, S] real<lower=0> pop;
}

parameters {
  matrix<lower=0>[S, T] lambda;
  vector<lower=0>[S]    lambda_pre;

  // 1번 그룹은 고정, 2~17 그룹(16개)만 추정
  vector<lower=1>[A5 - 1] kappa_age_rest;
  vector<lower=0, upper=1>[A5] omega_age;

  real<lower=0, upper=1> delta_primary;
  real<lower=0, upper=1> delta_secondary;
  real<lower=0, upper=1> delta_tert_frac;
  real<lower=0, upper=1> delta_quat_frac;

  real          mu_rho_global;
  real<lower=0> sigma_rho_state;
  vector[S]     z_rho_state;
  real<lower=0> sigma_rho_rw;
  matrix[S, T]  z_rho_rw;

  vector<lower=0>[S] phi_cases;
}

transformed parameters {
  vector[A5] kappa_age;
  kappa_age[1] = 1.0;                     // 1번 그룹 1.0 고정
  kappa_age[2:A5] = kappa_age_rest;       // 나머지 채움
  
  vector[A1] kappa_full;
  for (i in 1:A1) kappa_full[i] = kappa_age[band_of_age[i]];

  vector[4] delta;
  delta[1] = delta_primary;
  delta[2] = delta_secondary;
  delta[3] = delta_primary * delta_tert_frac;
  delta[4] = delta[3] * delta_quat_frac;

  matrix[S, T] logit_rho;
  for (s in 1:S) {
    real base = mu_rho_global + sigma_rho_state * z_rho_state[s];
    logit_rho[s, 1] = base + sigma_rho_rw * z_rho_rw[s, 1];
    for (t in 2:T)
      logit_rho[s, t] = logit_rho[s, t - 1] + sigma_rho_rw * z_rho_rw[s, t];
  }
  matrix<lower=0, upper=1>[S, T] rho;
  for (s in 1:S)
    for (t in 1:T)
      rho[s, t] = inv_logit(logit_rho[s, t]);
      
  // *주의: 모델 속도 최적화를 위해 Lambda 계산은 여기서 하지 않고 model 블록 내 지역 변수로 옮겼습니다.
}

model {
  to_vector(lambda) ~ exponential(10);
  lambda_pre        ~ exponential(10);
  
  kappa_age_rest    ~ normal(1, 1);
  omega_age         ~ normal(0.5, 0.5);

  delta_primary    ~ beta(2, 8);
  delta_secondary  ~ beta(5, 5);
  delta_tert_frac  ~ beta(1, 1);
  delta_quat_frac  ~ beta(1, 1);

  mu_rho_global       ~ normal(logit(0.15), 0.7);
  sigma_rho_state     ~ normal(0, 1);
  z_rho_state         ~ std_normal();
  sigma_rho_rw        ~ normal(0, 0.3);
  to_vector(z_rho_rw) ~ std_normal();

  phi_cases ~ normal(0, 10);

  for (s in 1:S) {
    // 1. 메모리 절약을 위한 지역변수 할당 및 누적 감염력 계산
    array[A1, T + 1] real Lam_local;
    Lam_local[1, 1] = 0.0;
    for (i in 2:A1) {
      Lam_local[i, 1] = Lam_local[i - 1, 1] + lambda_pre[s] * kappa_full[i - 1];
    }
    
    for (t in 1:T) {
      Lam_local[1, t + 1] = kappa_full[1] * lambda[s, t];
      for (i in 2:(A1 - 1)) {
        Lam_local[i, t + 1] = Lam_local[i - 1, t] + kappa_full[i] * lambda[s, t];
      }
      Lam_local[A1, t + 1] = Lam_local[A1, t] + kappa_full[A1] * lambda[s, t];
    }
    
    // 2. 우도(Likelihood) 계산 (C++ -expm1 최적화 적용)
    for (t in 1:T) {
      array[A5] real band_sum;
      array[A5] int  band_n;
      for (b in 1:A5) { band_sum[b] = 0.0; band_n[b] = 0; }
      
      vector[A1] delta_haz_vec = kappa_full * lambda[s, t];

      for (i in 1:A1) {
        real Lam_before;
        if (i == 1)       Lam_before = 0.0;
        else if (i == A1) Lam_before = Lam_local[A1, t];
        else              Lam_before = Lam_local[i - 1, t];

        real p_before  = -expm1(-Lam_before); // 1-exp(-x) 대체
        real delta_haz = delta_haz_vec[i];

        real q_before = 1 - p_before;
        real p1 = p_before;
        real p2 = p1 * p_before;
        real p3 = p2 * p_before;
        real q1 = q_before;
        real q2 = q1 * q_before;
        real q3 = q2 * q_before;
        real q4 = q3 * q_before;

        real pk1 = q4;
        real pk2 = 4  * p1 * q3;
        real pk3 = 6  * p2 * q2;
        real pk4 = 4  * p3 * q1;

        real apparent_i =
          delta[1] * pk1 * (-expm1(-4 * delta_haz)) +
          delta[2] * pk2 * (-expm1(-3 * delta_haz)) +
          delta[3] * pk3 * (-expm1(-2 * delta_haz)) +
          delta[4] * pk4 * (-expm1(-1 * delta_haz));
          
        int b = band_of_age[i];
        band_sum[b] += apparent_i;
        band_n[b]   += 1;
      }

      for (b in 1:A5) {
        real mu_cases = fmax(rho[s, t] * omega_age[b] * (band_sum[b] / band_n[b]) * pop[b, t, s], 1e-8);
        cases[b, t, s] ~ neg_binomial_2(mu_cases, phi_cases[s]);
      }
    }
  }
}

generated quantities {
  array[A5, T, S] int cases_pred;
  array[A5, T, S] real<lower=0, upper=1> pi_exposed;

  for (s in 1:S) {
    // 생성된 변수 블록에서도 Lam_local 재계산 (메모리 부담 해소)
    array[A1, T + 1] real Lam_local;
    Lam_local[1, 1] = 0.0;
    for (i in 2:A1) Lam_local[i, 1] = Lam_local[i - 1, 1] + lambda_pre[s] * kappa_full[i - 1];
    
    for (t in 1:T) {
      Lam_local[1, t + 1] = kappa_full[1] * lambda[s, t];
      for (i in 2:(A1 - 1)) Lam_local[i, t + 1] = Lam_local[i - 1, t] + kappa_full[i] * lambda[s, t];
      Lam_local[A1, t + 1] = Lam_local[A1, t] + kappa_full[A1] * lambda[s, t];
    }
    
    for (t in 1:T) {
      array[A5] real band_sum_app;
      array[A5] int  band_n_app;
      array[A5] real band_sum_p0;
      array[A5] int  band_n_p0;
      
      for (b in 1:A5) { 
        band_sum_app[b] = 0.0; band_n_app[b] = 0; 
        band_sum_p0[b]  = 0.0; band_n_p0[b]  = 0; 
      }
      
      vector[A1] delta_haz_vec = kappa_full * lambda[s, t];
      
      for (i in 1:A1) {
        real Lam_before;
        if (i == 1)       Lam_before = 0.0;
        else if (i == A1) Lam_before = Lam_local[A1, t];
        else              Lam_before = Lam_local[i - 1, t];
        
        real p_before  = -expm1(-Lam_before);
        real delta_haz = delta_haz_vec[i];
        
        real q_before = 1 - p_before;
        real p1 = p_before;
        real p2 = p1 * p_before;
        real p3 = p2 * p_before;
        real q1 = q_before;
        real q2 = q1 * q_before;
        real q3 = q2 * q_before;
        real q4 = q3 * q_before;
        
        real pk1 = q4;            
        real pk2 = 4  * p1 * q3;  
        real pk3 = 6  * p2 * q2;  
        real pk4 = 4  * p3 * q1;  
        
        real apparent_i =
          delta[1] * pk1 * (-expm1(-4 * delta_haz)) +
          delta[2] * pk2 * (-expm1(-3 * delta_haz)) +
          delta[3] * pk3 * (-expm1(-2 * delta_haz)) +
          delta[4] * pk4 * (-expm1(-1 * delta_haz));
        
        int b = band_of_age[i];
        band_sum_app[b] += apparent_i;
        band_n_app[b]   += 1;
        
        real p_end = -expm1(-Lam_local[i, t + 1]);
        real q_end = 1 - p_end;
        real p0    = q_end * q_end * q_end * q_end;   
        band_sum_p0[b] += p0;
        band_n_p0[b]   += 1;
      }
      
      for (b in 1:A5) {
        real mu_cases = fmax(rho[s, t] * omega_age[b] * (band_sum_app[b] / band_n_app[b]) * pop[b, t, s], 1e-8);
        cases_pred[b, t, s]  = neg_binomial_2_rng(mu_cases, phi_cases[s]);
        pi_exposed[b, t, s]  = 1 - band_sum_p0[b] / band_n_p0[b];
      }
    }
  }
}

