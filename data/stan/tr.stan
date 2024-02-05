functions {
  real sms(real x) {
    return sqrt(1 - x*x);
  }

  matrix getChol(real a) {
    real A = sms(a);

    matrix[6,6] chol = [
      [1, 0],
      [a, A]
    ];
    return chol;
  }

}

data {
  real lb;
  real ub;

  int<lower=0> N;
  array[N] vector[2] x;

  int<lower=0> Nu;
  int<lower=0> Nl;
  array[Nu] int iu;
  array[Nl] int il;

  real muPriorMean;
  real muPriorSd;
  real sigmaPriorLogMean;
}

parameters {
  array[Nu] real<lower=ub> xu;
  array[Nl] real<upper=lb> xl;

  real mu;
  real<lower=0> sigma;
  real<lower=-1,upper=1> a;
}

model {
  vector[2] Mu = [mu, mu]';
  vector[2] sigmaV = [sigma, sigma]';

  matrix[2,2] chol = getChol(a);

  target += lkj_corr_cholesky_lpdf(chol | 1);

  real sum_;
  for (i in 2:2) {
    for (j in 1:i-1) {
      sum_ = 1;
      for (k in 1:j-1) {
        sum_ -= square(chol[i,k]);
      }
      target += log(sum_) / 2;
    }
  }

  array[N] vector[2] xAll = x;
  for (i in 1:Nu) {
    xAll[(iu[i] - 1) %/% 2 + 1, (iu[i] - 1) % 2 + 1] = xu[i];
  }
  for (i in 1:Nl) {
    xAll[(il[i] - 1) %/% 2 + 1, (il[i] - 1) % 2 + 1] = xl[i];
  }

  target += multi_normal_cholesky_lpdf(xAll | Mu, diag_matrix(sigmaV) * chol);

  target += normal_lpdf(mu | muPriorMean, muPriorSd);
  target += lognormal_lpdf(sigma | sigmaPriorLogMean, 1);
}

// generated quantities {
//   matrix[6,6] corrMat = multiply_lower_tri_self_transpose(getChol(a));
// }
