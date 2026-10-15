
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
  matrix<lower=0>[S, T] lambda;          // per-serotype FOI, state x year
  vector<lower=0>[S]    lambda_pre;      // constant pre-2000 per-serotype FOI per state

  vector<lower=0>[A5] kappa_age;        // age-FOI multiplier, all 17 bands now free (prior anchors the scale)
  vector<lower=0, upper=1>[A5] omega_age;   // national apparent-illness / care-seeking age multiplier

  real<lower=0, upper=1> delta_primary;
  real<lower=0, upper=1> delta_secondary;
  // CHAINED shrinkage: delta_3 = delta_1 * frac_tert,
  //                     delta_4 = delta_3 * frac_quat
  // Each frac in (0,1) -> structurally guarantees
  // delta_primary > delta_tertiary > delta_quaternary
  // secondary, which stays a free, unconstrained parameter since it's
  // driven by a different mechanism -- ADE-style enhancement -- not
  // progressive shrinkage).
  real<lower=0, upper=1> delta_tert_frac;
  real<lower=0, upper=1> delta_quat_frac;

  real          mu_rho_global;
  real<lower=0> sigma_rho_state;
  vector[S]     z_rho_state;
  real<lower=0> sigma_rho_rw;
  matrix[S, T]  z_rho_rw;

  vector<lower=0>[S] phi_cases;          // NegBin overdispersion, PER STATE
}

transformed parameters {
  vector[A1] kappa_full;
  for (i in 1:A1) kappa_full[i] = kappa_age[band_of_age[i]];

  vector[4] delta;
  delta[1] = delta_primary;
  delta[2] = delta_secondary;
  delta[3] = delta_primary * delta_tert_frac;   // 1차 > 3차
  delta[4] = delta[3] * delta_quat_frac;        // 3차 > 4차 (now chained off delta[3], not delta_primary)

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

  // Lambda[i, slot, s]: slot 1 = pre-2000 baseline (cumulative hazard
  // accumulated under lambda_pre[s] alone, for someone about to turn
  // single-year-age i). slot (t+1) = cumulative hazard as of the END
  // of model year t (t = 1..T <-> calendar year years[t]).
  array[A1, T + 1, S] real<lower=0> Lambda;
  print(kappa_age);
  for (s in 1:S) {
    Lambda[1, 1, s] = 0.0;
    for (i in 2:A1)
      Lambda[i, 1, s] = Lambda[i - 1, 1, s] + lambda_pre[s] * kappa_full[i - 1];

    for (t in 1:T) {
      Lambda[1, t + 1, s] = kappa_full[1] * lambda[s, t];
      for (i in 2:(A1 - 1))
        Lambda[i, t + 1, s] = Lambda[i - 1, t, s] + kappa_full[i] * lambda[s, t];
      Lambda[A1, t + 1, s] = Lambda[A1, t, s] + kappa_full[A1] * lambda[s, t];  // 80+ self-loop
    }
  }

  // expected reported cases, band x year x state
  array[A5, T, S] real<lower=0> mu_cases;
  {
    for (s in 1:S) {
      for (t in 1:T) {
        array[A5] real band_sum;
        array[A5] int  band_n;
        for (b in 1:A5) { band_sum[b] = 0.0; band_n[b] = 0; }

        for (i in 1:A1) {
          real Lam_before;
          if (i == 1)       Lam_before = 0.0;
          else if (i == A1) Lam_before = Lambda[A1, t, s];
          else              Lam_before = Lambda[i - 1, t, s];

          real p_before  = 1 - exp(-Lam_before);
          real delta_haz = kappa_full[i] * lambda[s, t];

          // n=4 fixed -> binomial coefficients are constants (1,4,6,4).
          // Closed form avoids binomial_lpmf's internal lgamma calls,
          // which get evaluated S*T*A1*4 times per gradient step.
          real q_before = 1 - p_before;
          real p1 = p_before;
          real p2 = p1 * p_before;
          real p3 = p2 * p_before;
          real q1 = q_before;
          real q2 = q1 * q_before;
          real q3 = q2 * q_before;
          real q4 = q3 * q_before;

          real pk1 = q4;            // P(k=0) = C(4,0) p^0 q^4
          real pk2 = 4  * p1 * q3;  // P(k=1) = C(4,1) p^1 q^3
          real pk3 = 6  * p2 * q2;  // P(k=2) = C(4,2) p^2 q^2
          real pk4 = 4  * p3 * q1;  // P(k=3) = C(4,3) p^3 q^1

          real apparent_i =
            delta[1] * pk1 * (1 - exp(-4 * delta_haz)) +
            delta[2] * pk2 * (1 - exp(-3 * delta_haz)) +
            delta[3] * pk3 * (1 - exp(-2 * delta_haz)) +
            delta[4] * pk4 * (1 - exp(-1 * delta_haz));
          int b = band_of_age[i];
          band_sum[b] += apparent_i;
          band_n[b]   += 1;
        }

        for (b in 1:A5)
          mu_cases[b, t, s] = fmax(
            rho[s, t] * omega_age[b] * (band_sum[b] / band_n[b]) * pop[b, t, s],
            1e-8
          );
      }
    }
  }
}

model {
  to_vector(lambda) ~ exponential(10);   // mean 0.1; 1차감염 hazard 4*lambda ~ exponential(2.5), mean 0.4
  lambda_pre        ~ exponential(10);
  kappa_age         ~ normal(1, 1);      // <lower=0>라 사실상 0에서 잘린 normal(1,1), 17개 밴드 전부
  omega_age         ~ normal(0.5, 0.5);  // <lower=0,upper=1>이라 [0,1]에서 잘린 normal(0.5,0.5)

  delta_primary    ~ beta(2, 8);    // prior mean 0.20 (matches supplement)
  delta_secondary  ~ beta(5, 5);    // prior mean 0.50 (matches supplement)
  delta_tert_frac  ~ beta(1, 1);    // flat on (0,1) -- delta_tertiary   = delta_primary   * this
  delta_quat_frac  ~ beta(1, 1);    // flat on (0,1) -- delta_quaternary = delta_tertiary * this (chained)

  mu_rho_global       ~ normal(logit(0.15), 0.7);  // PLACEHOLDER central reporting-rate prior
  sigma_rho_state     ~ normal(0, 1);
  z_rho_state         ~ std_normal();
  sigma_rho_rw        ~ normal(0, 0.3);
  to_vector(z_rho_rw) ~ std_normal();

  phi_cases ~ normal(0, 10);

  for (s in 1:S)
    for (t in 1:T)
      for (b in 1:A5)
        cases[b, t, s] ~ neg_binomial_2(mu_cases[b, t, s], phi_cases[s]);
}

generated quantities {
  array[A5, T, S] int cases_pred;
  array[A5, T, S] real<lower=0, upper=1> pi_exposed;   // P(>=1 prior infection), band-averaged

  for (s in 1:S) {
    for (t in 1:T) {
      array[A5] real band_sum_p0;
      array[A5] int  band_n2;
      for (b in 1:A5) { band_sum_p0[b] = 0.0; band_n2[b] = 0; }

      for (i in 1:A1) {
        real p_end = 1 - exp(-Lambda[i, t + 1, s]);
        real q_end = 1 - p_end;
        real p0    = q_end * q_end * q_end * q_end;   // P(k=0) = q^4
        int  b     = band_of_age[i];
        band_sum_p0[b] += p0;
        band_n2[b]     += 1;
      }
      for (b in 1:A5) {
        pi_exposed[b, t, s]  = 1 - band_sum_p0[b] / band_n2[b];
        cases_pred[b, t, s]  = neg_binomial_2_rng(mu_cases[b, t, s], phi_cases[s]);
      }
    }
  }
}

