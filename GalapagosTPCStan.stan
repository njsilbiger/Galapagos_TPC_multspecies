// Saved as GalapagosTPCStan.stan
data {
  int<lower=1> N;  // number of observations
  vector[N] Y;  // response variable
  int<lower=1> K_lnc;  // number of population-level effects
  matrix[N, K_lnc] X_lnc;  // population-level design matrix
  int<lower=1> K_Eh;  // number of population-level effects
  matrix[N, K_Eh] X_Eh;  // population-level design matrix
  int<lower=1> K_Th;  // number of population-level effects
  matrix[N, K_Th] X_Th;  // population-level design matrix
  int<lower=1> K_E;  // number of population-level effects
  matrix[N, K_E] X_E;  // population-level design matrix
  // covariate vectors
  vector[N] C_1;
  // data for group-level effects of ID 1
  int<lower=1> N_1;
  int<lower=1> M_1;
  int<lower=1> J_1[N];
  vector[N] Z_1_lnc_1;
  int prior_only;  // should the likelihood be ignored?
}
transformed data {
}
parameters {
  vector<lower=0>[K_lnc] b_lnc;  // population-level effects
  vector<lower=0>[K_Eh] b_Eh;  // population-level effects
  vector<lower=0>[K_Th] b_Th;  // population-level effects
  vector<lower=0>[K_E] b_E;  // population-level effects
  
real<lower=0> sigma;  // residual SD
  vector<lower=0>[M_1] sd_1;  // group-level standard deviations
  vector[N_1] z_1[M_1];  // unscaled group-level effects

}
transformed parameters {
  // group-level effects
  vector[N_1] r_1_lnc_1 = (sd_1[1] * (z_1[1]));
}
model {
  vector[N] nlp_lnc = X_lnc * b_lnc;
  vector[N] nlp_Eh = X_Eh * b_Eh;
  vector[N] nlp_Th = X_Th * b_Th;
  vector[N] nlp_E = X_E * b_E;
  vector[N] mu;
  for (n in 1:N) {
    nlp_lnc[n] += r_1_lnc_1[J_1[n]] * Z_1_lnc_1[n];
    // compute non-linear predictor
    mu[n] = nlp_lnc[n] + log(exp(nlp_E[n] / 0.0000862 * (1 / 299.15 - 1 / C_1[n]))) + log(1 / (1 + exp(nlp_Eh[n] / 0.0000862 * (1 / nlp_Th[n] - 1 / C_1[n]))));
  }
  // priors including all constants
  target += normal_lpdf(b_lnc | 0, 10);
  target += normal_lpdf(b_Eh | 0, 10);
  target += normal_lpdf(b_Th | 320, 10);
  target += normal_lpdf(b_E | 0, 10);
  target += student_t_lpdf(sigma | 3, 0, 10)
  - 1 * student_t_lccdf(0 | 3, 0, 10);
  target += student_t_lpdf(sd_1 | 3, 0, 10)
  - 1 * student_t_lccdf(0 | 3, 0, 10);
  target += normal_lpdf(z_1[1] | 0, 1);
  // likelihood including all constants
  if (!prior_only) {
    target += normal_lpdf(Y | mu, sigma);
  }
}
generated quantities {
 // real<lower=1> Topt;
  real Topt;
   Topt =  (b_Eh[1] * b_Th[1]) / (b_Eh[1] + (0.0000862 * b_Th[1] * log((b_Eh[1] / b_E[1]) - 1)));
}