functions {
  real partial_sum(array[] int slice_s,
                   int start, int end,
                   int T, int A1, int A5,
                   array[] int band_of_age,
                   array[,,] int cases,
                   array[,,] real pop,
                   matrix lambda,
                   vector lambda_pre,
                   matrix rho,
                   vector phi_cases,
                   vector kappa_full,
                   vector omega_age,
                   vector delta) {
    real target_acc = 0.0;
    for (idx in 1:size(slice_s)) {
      int s = slice_s[idx];
 
      array[A1, T + 1] real Lam;
      Lam[1, 1] = 0.0;
      for (i in 2:A1)
        Lam[i, 1] = Lam[i-1, 1] + lambda_pre[s] * kappa_full[i-1];
      for (t in 1:T) {
        Lam[1, t+1] = kappa_full[1] * lambda[s, t];
        for (i in 2:(A1-1))
          Lam[i, t+1] = Lam[i-1, t] + kappa_full[i] * lambda[s, t];
        Lam[A1, t+1] = Lam[A1, t] + kappa_full[A1] * lambda[s, t];
      }
 
      for (t in 1:T) {
        array[A5] real band_sum;
        array[A5] int  band_n;
        for (b in 1:A5) { band_sum[b] = 0.0; band_n[b] = 0; }
        vector[A1] dh = kappa_full * lambda[s, t];
 
        for (i in 1:A1) {
          real Lb;
          if (i == 1)       Lb = 0.0;
          else if (i == A1) Lb = Lam[A1, t];
          else              Lb = Lam[i-1, t];
 
          real p = -expm1(-Lb);
          real q = 1 - p;
          real p2 = p*p; real p3 = p2*p;
          real q2 = q*q; real q3 = q2*q; real q4 = q3*q;
 
          real app =
            delta[1] *      q4       * (-expm1(-4 * dh[i])) +
            delta[2] * (4*p *q3)     * (-expm1(-3 * dh[i])) +
            delta[3] * (6*p2*q2)     * (-expm1(-2 * dh[i])) +
            delta[4] * (4*p3*q)      * (-expm1(-1 * dh[i]));
 
          int b = band_of_age[i];
          band_sum[b] += app;
          band_n[b]   += 1;
        }
 
        for (b in 1:A5) {
          real mu = fmax(rho[s,t] * omega_age[b] * (band_sum[b]/band_n[b]) * pop[b,t,s], 1e-8);
          target_acc += neg_binomial_2_lpmf(cases[b,t,s] | mu, phi_cases[s]);
        }
      }
    }
    return target_acc;
  }
}
 
data {
  int<lower=1> S;
  int<lower=1> A5;
  int<lower=1> A1;
  int<lower=1> T;
  array[A1] int<lower=1, upper=A5> band_of_age;
  array[A5, T, S] int<lower=0>    cases;
  array[A5, T, S] real<lower=0>   pop;
}
 
transformed data {
  array[S] int seq_S;
  for (s in 1:S) seq_S[s] = s;
}
 
parameters {
  matrix<lower=0, upper=1>[S, T]   lambda;
  vector<lower=0, upper=1>[S]      lambda_pre;
 
  vector<lower=0>[A5]              kappa_age;   // 전 연령군 자유 추정, normal(1,2) prior
  vector<lower=0, upper=1>[A5]     omega_age;   // 전 연령군 자유 추정, normal(0.5,1) prior
 
  real<lower=0, upper=1> delta_primary;
  real<lower=0, upper=1> delta_secondary;
  real<lower=0, upper=1> delta_tert_frac;
  real<lower=0, upper=1> delta_quat_frac;
 
  matrix<lower=0, upper=1>[S, T]   rho;
  vector<lower=0>[S]               phi_cases;
}
 
transformed parameters {
  vector[A1] kappa_full;
  for (i in 1:A1) kappa_full[i] = kappa_age[band_of_age[i]];
 
  vector[4] delta;
  delta[1] = delta_primary;
  delta[2] = delta_secondary;
  delta[3] = delta_primary * delta_tert_frac;
  delta[4] = delta[3]      * delta_quat_frac;
}
 
model {
  to_vector(lambda) ~ exponential(2);
  lambda_pre        ~ exponential(40);
 
  kappa_age ~ normal(1, 0.5);     // <lower=0> 제약 안에서 양수 유지
  omega_age ~ normal(0.5, 0.5);   // <lower=0,upper=1> 제약 안에서
 
  delta_primary   ~ beta(2, 8);
  delta_secondary ~ beta(5, 5);
  delta_tert_frac ~ beta(1, 1);
  delta_quat_frac ~ beta(1, 1);
 
  to_vector(rho) ~ beta(1, 2);
 
  phi_cases ~ normal(0, 10);
 
  target += reduce_sum(partial_sum, seq_S, 1,
                       T, A1, A5, band_of_age, cases, pop,
                       lambda, lambda_pre, rho, phi_cases,
                       kappa_full, omega_age, delta);
}
 
generated quantities {
  array[A5, T, S] int   cases_pred;
  array[A5, T, S] real  pi_exposed;
 
  for (s in 1:S) {
    array[A1, T+1] real Lam;
    Lam[1,1] = 0.0;
    for (i in 2:A1) Lam[i,1] = Lam[i-1,1] + lambda_pre[s] * kappa_full[i-1];
    for (t in 1:T) {
      Lam[1,t+1] = kappa_full[1] * lambda[s,t];
      for (i in 2:(A1-1)) Lam[i,t+1] = Lam[i-1,t] + kappa_full[i] * lambda[s,t];
      Lam[A1,t+1] = Lam[A1,t] + kappa_full[A1] * lambda[s,t];
    }
 
    for (t in 1:T) {
      array[A5] real bapp; array[A5] int napp;
      array[A5] real bp0;  array[A5] int np0;
      for (b in 1:A5) { bapp[b]=0.0; napp[b]=0; bp0[b]=0.0; np0[b]=0; }
      vector[A1] dh = kappa_full * lambda[s,t];
 
      for (i in 1:A1) {
        real Lb;
        if (i==1)       Lb = 0.0;
        else if (i==A1) Lb = Lam[A1,t];
        else            Lb = Lam[i-1,t];
 
        real p = -expm1(-Lb); real q = 1-p;
        real p2=p*p; real p3=p2*p;
        real q2=q*q; real q3=q2*q; real q4=q3*q;
 
        real app =
          delta[1]*q4       *(-expm1(-4*dh[i])) +
          delta[2]*(4*p*q3) *(-expm1(-3*dh[i])) +
          delta[3]*(6*p2*q2)*(-expm1(-2*dh[i])) +
          delta[4]*(4*p3*q) *(-expm1(-1*dh[i]));
 
        int b = band_of_age[i];
        bapp[b] += app; napp[b] += 1;
 
        real pe = -expm1(-Lam[i,t+1]); real qe=1-pe;
        bp0[b] += qe*qe*qe*qe; np0[b] += 1;
      }
 
      for (b in 1:A5) {
        real mu = fmax(rho[s,t]*omega_age[b]*(bapp[b]/napp[b])*pop[b,t,s], 1e-8);
        cases_pred[b,t,s] = neg_binomial_2_rng(mu, phi_cases[s]);
        pi_exposed[b,t,s] = 1 - bp0[b]/np0[b];
      }
    }
  }
}

