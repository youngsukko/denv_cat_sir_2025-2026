// Bayesian age-structured catalytic model of dengue infection history,
// Brazil, 27 federative units, 2000-2025 (see Methods of the paper).
//
// Paper notation -> code:
//   per-serotype FOI by state and year        lambda~_s(t) -> lambda[s, t]
//   constant per-serotype FOI before 2000     lambda~_s(0) -> lambda_pre[s]
//   age-specific FOI multiplier               kappa_a      -> kappa_age[a]
//   age-specific reporting multiplier         omega_a      -> omega_age[a]
//   state- and year-specific reporting rate   rho~_s(t)    -> rho[s, t]
//   apparent illness ratio, k-th infection    delta_k      -> delta[k]
//   fixed negative-binomial overdispersion    phi_s        -> phi_fixed[s]
//
// Exposure is accumulated over single-year age cohorts (A1 = 81, ages 0..80+)
// that are mapped to 17 five-year age bands (A5) through band_of_age.

functions {
  // Cumulative per-serotype hazard Lam[i, c] of single-year cohort i.
  // Column 1 is the pre-2000 baseline, column t + 1 is the end of model year t.
  array[,] real cumulative_hazard(int T, int A1, vector kappa_full,
                                  row_vector lambda_s, real lambda_pre_s) {
    array[A1, T + 1] real Lam;
    Lam[1, 1] = 0.0;
    for (i in 2:A1)
      Lam[i, 1] = Lam[i - 1, 1] + lambda_pre_s * kappa_full[i - 1];
    for (t in 1:T) {
      Lam[1, t + 1] = kappa_full[1] * lambda_s[t];
      for (i in 2:(A1 - 1))
        Lam[i, t + 1] = Lam[i - 1, t] + kappa_full[i] * lambda_s[t];
      // the 80+ cohort keeps accumulating exposure
      Lam[A1, t + 1] = Lam[A1, t] + kappa_full[A1] * lambda_s[t];
    }
    return Lam;
  }

  // Log-likelihood of a slice of states (used with reduce_sum).
  real partial_sum(array[] int slice_s,
                   int start, int end,
                   int T, int A1, int A5,
                   array[] int band_of_age,
                   array[,,] int cases,
                   array[,,] real pop,
                   matrix lambda,
                   vector lambda_pre,
                   matrix rho,
                   array[] real phi_fixed,
                   vector kappa_full,
                   vector omega_age,
                   vector delta,
                   vector pop_weight) {
    real target_acc = 0.0;
    for (idx in 1:size(slice_s)) {
      int s = slice_s[idx];
      array[A1, T + 1] real Lam = cumulative_hazard(T, A1, kappa_full,
                                                    lambda[s], lambda_pre[s]);

      for (t in 1:T) {
        array[A5] real band_sum;
        array[A5] int  band_n;
        for (b in 1:A5) { band_sum[b] = 0.0; band_n[b] = 0; }
        vector[A1] dh = kappa_full * lambda[s, t];

        for (i in 1:A1) {
          // exposure at the start of year t
          real Lb;
          if (i == 1)       Lb = 0.0;
          else if (i == A1) Lb = Lam[A1, t];
          else              Lb = Lam[i - 1, t];

          // probability of having been infected by a given serotype
          real p = -expm1(-Lb);
          real q = 1 - p;
          real p2 = p * p; real p3 = p2 * p;
          real q2 = q * q; real q3 = q2 * q; real q4 = q3 * q;

          // P(j prior infections) x P(>= 1 new infection this year) x delta_{j+1}
          real app =
            delta[1] *        q4     * (-expm1(-4 * dh[i])) +
            delta[2] * (4 * p  * q3) * (-expm1(-3 * dh[i])) +
            delta[3] * (6 * p2 * q2) * (-expm1(-2 * dh[i])) +
            delta[4] * (4 * p3 * q)  * (-expm1(-1 * dh[i]));

          int b = band_of_age[i];
          band_sum[b] += app;
          band_n[b]   += 1;
        }

        for (b in 1:A5) {
          real mu = fmax(rho[s, t] * omega_age[b] * (band_sum[b] / band_n[b]) * pop[b, t, s], 1e-8);
          target_acc += pop_weight[s] * neg_binomial_2_lpmf(cases[b, t, s] | mu, phi_fixed[s]);
        }
      }
    }
    return target_acc;
  }
}

data {
  int<lower=1> S;                              // number of states
  int<lower=1> A5;                             // number of five-year age bands
  int<lower=1> A1;                             // number of single-year cohorts
  int<lower=1> T;                              // number of years (2000-2025)
  array[A1] int<lower=1, upper=A5> band_of_age;
  array[A5, T, S] int<lower=0>     cases;      // reported cases
  array[A5, T, S] real<lower=0>    pop;        // population
  vector<lower=0>[S]               pop_weight; // log-likelihood weight per state
  array[S] real<lower=0>           phi_fixed;  // overdispersion per state
}

transformed data {
  array[S] int seq_S;
  for (s in 1:S) seq_S[s] = s;
}

parameters {
  matrix<lower=0>[S, T] lambda;
  vector<lower=0>[S]    lambda_pre;

  vector<lower=0>[A5]          kappa_age;
  vector<lower=0, upper=1>[A5] omega_age;

  real<lower=0, upper=1> delta_primary;
  real<lower=0, upper=1> delta_secondary;
  real<lower=0, upper=1> delta_tert_frac;  // c_tert
  real<lower=0, upper=1> delta_quat_frac;  // c_quat

  matrix<lower=0, upper=1>[S, T] rho;
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
  // FOI priors scaled by the mean age multiplier (kappa itself has no prior)
  to_vector(lambda) ~ exponential(4 / mean(kappa_age));
  lambda_pre        ~ exponential(40 / mean(kappa_age));

  to_vector(rho) ~ beta(1, 1.5);

  // nationally shared parameters: prior log-densities scaled by S
  target += S * beta_lpdf(omega_age[1] | 1.5, 1);
  for (i in 2:A5) target += S * normal_lpdf(omega_age[i] | omega_age[i - 1], 1);

  target += S * beta_lpdf(delta_primary   | 4, 16);
  target += S * beta_lpdf(delta_secondary | 5, 5);

  target += reduce_sum(partial_sum, seq_S, 1,
                       T, A1, A5, band_of_age, cases, pop,
                       lambda, lambda_pre, rho, phi_fixed,
                       kappa_full, omega_age, delta,
                       pop_weight);
}

generated quantities {
  array[A5, T, S] int  cases_pred;    // posterior predictive reported cases
  array[A5, T, S] real inf_total;     // expected new infections (reported or not)
  array[A5, T, S] real seropos_0;     // proportion with 0..4 prior infections
  array[A5, T, S] real seropos_1;     // at the end of each year
  array[A5, T, S] real seropos_2;
  array[A5, T, S] real seropos_3;
  array[A5, T, S] real seropos_4;
  array[A5, S]    real seropos_pre_0; // same at the start of 2000
  array[A5, S]    real seropos_pre_1;
  array[A5, S]    real seropos_pre_2;
  array[A5, S]    real seropos_pre_3;
  array[A5, S]    real seropos_pre_4;

  for (s in 1:S) {
    array[A1, T + 1] real Lam = cumulative_hazard(T, A1, kappa_full,
                                                  lambda[s], lambda_pre[s]);

    // pre-2000 baseline (column 1 of Lam)
    {
      array[A5] real bx0; array[A5] real bx1; array[A5] real bx2;
      array[A5] real bx3; array[A5] real bx4; array[A5] int nx;
      for (b in 1:A5) { bx0[b] = 0.0; bx1[b] = 0.0; bx2[b] = 0.0; bx3[b] = 0.0; bx4[b] = 0.0; nx[b] = 0; }

      for (i in 1:A1) {
        real pe = -expm1(-Lam[i, 1]); real qe = 1 - pe;
        real pe2 = pe * pe; real pe3 = pe2 * pe;
        real qe2 = qe * qe; real qe3 = qe2 * qe; real qe4 = qe3 * qe;

        int b = band_of_age[i];
        bx0[b] += qe4;
        bx1[b] += 4 * pe * qe3;
        bx2[b] += 6 * pe2 * qe2;
        bx3[b] += 4 * pe3 * qe;
        bx4[b] += pe2 * pe2;
        nx[b]  += 1;
      }
      for (b in 1:A5) {
        seropos_pre_0[b, s] = bx0[b] / nx[b];
        seropos_pre_1[b, s] = bx1[b] / nx[b];
        seropos_pre_2[b, s] = bx2[b] / nx[b];
        seropos_pre_3[b, s] = bx3[b] / nx[b];
        seropos_pre_4[b, s] = bx4[b] / nx[b];
      }
    }

    // 2000-2025: predicted cases, total infections, immune profile
    for (t in 1:T) {
      array[A5] real bapp; array[A5] int napp;
      array[A5] real binf;
      array[A5] real bx0; array[A5] real bx1; array[A5] real bx2;
      array[A5] real bx3; array[A5] real bx4; array[A5] int nx;
      for (b in 1:A5) {
        bapp[b] = 0.0; napp[b] = 0; binf[b] = 0.0;
        bx0[b] = 0.0; bx1[b] = 0.0; bx2[b] = 0.0; bx3[b] = 0.0; bx4[b] = 0.0; nx[b] = 0;
      }
      vector[A1] dh = kappa_full * lambda[s, t];

      for (i in 1:A1) {
        real Lb;
        if (i == 1)       Lb = 0.0;
        else if (i == A1) Lb = Lam[A1, t];
        else              Lb = Lam[i - 1, t];

        real p = -expm1(-Lb); real q = 1 - p;
        real p2 = p * p; real p3 = p2 * p;
        real q2 = q * q; real q3 = q2 * q; real q4 = q3 * q;

        real app =
          delta[1] *        q4     * (-expm1(-4 * dh[i])) +
          delta[2] * (4 * p  * q3) * (-expm1(-3 * dh[i])) +
          delta[3] * (6 * p2 * q2) * (-expm1(-2 * dh[i])) +
          delta[4] * (4 * p3 * q)  * (-expm1(-1 * dh[i]));

        real inf_any =
                 q4     * (-expm1(-4 * dh[i])) +
          (4 * p  * q3) * (-expm1(-3 * dh[i])) +
          (6 * p2 * q2) * (-expm1(-2 * dh[i])) +
          (4 * p3 * q)  * (-expm1(-1 * dh[i]));

        int b = band_of_age[i];
        bapp[b] += app; napp[b] += 1;
        binf[b] += inf_any;

        real pe = -expm1(-Lam[i, t + 1]); real qe = 1 - pe;
        real pe2 = pe * pe; real pe3 = pe2 * pe;
        real qe2 = qe * qe; real qe3 = qe2 * qe; real qe4 = qe3 * qe;

        bx0[b] += qe4;
        bx1[b] += 4 * pe * qe3;
        bx2[b] += 6 * pe2 * qe2;
        bx3[b] += 4 * pe3 * qe;
        bx4[b] += pe2 * pe2;
        nx[b]  += 1;
      }

      for (b in 1:A5) {
        real mu = fmax(rho[s, t] * omega_age[b] * (bapp[b] / napp[b]) * pop[b, t, s], 1e-8);
        cases_pred[b, t, s] = neg_binomial_2_rng(mu, phi_fixed[s]);
        inf_total[b, t, s]  = (binf[b] / napp[b]) * pop[b, t, s];
        seropos_0[b, t, s]  = bx0[b] / nx[b];
        seropos_1[b, t, s]  = bx1[b] / nx[b];
        seropos_2[b, t, s]  = bx2[b] / nx[b];
        seropos_3[b, t, s]  = bx3[b] / nx[b];
        seropos_4[b, t, s]  = bx4[b] / nx[b];
      }
    }
  }
}
