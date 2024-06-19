functions {
  real sms(real x) {
    return sqrt(1 - x*x);
  }

  matrix getChol(real a) {
    real A = sms(a);

    matrix[2,2] chol = [
      [1, 0],
      [a, A]
    ];
    return chol;
  }

}

data {
  int<lower=0> N;
  array[N] vector[2] x;

  vector[2] muPriorMean;
  vector[2] muPriorSd;
  vector[2] sigmaPriorLogMean;
}

parameters {
  vector[2] mu;
  vector<lower=0>[2] sigma;
  real<lower=-1,upper=1> a;
}

model {
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

  target += multi_normal_cholesky_lpdf(x | mu, diag_matrix(sigma) * chol);

  target += normal_lpdf(mu | muPriorMean, muPriorSd);
  target += lognormal_lpdf(sigma | sigmaPriorLogMean, 1);
}
