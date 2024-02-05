functions {
  real sms(real x) {
    return sqrt(1 - x*x);
  }

  matrix getChol(real a, real b, real d) {
    real A = sms(a);
    real B = sms(b);
    real c = (1-a)*b/(A*B);
    real C = sms(c);
    real D = sms(d);

    matrix[4,4] chol = [
      [1, 0, 0, 0],
      [a, A, 0, 0],
      [b, B*c, B*C, 0],
      [b, B*c, B*C*d, B*C*D]
    ];
    return chol;
  }

}

data {
  array[2] real lb;
  array[2] real ub;

  int<lower=0> N;
  array[N] vector[4] x;

  array[2] int<lower=0> Nu;
  array[2] int<lower=0> Nl;
  array[Nu[1]] int iu1;
  array[Nl[1]] int il1;
  array[Nu[2]] int iu2;
  array[Nl[2]] int il2;

  vector[2] muPriorMean;
  vector[2] muPriorSd;
  vector[2] sigmaPriorLogMean;
}

parameters {
  array[Nu[1]] real<lower=ub[1]> xu1;
  array[Nl[1]] real<upper=lb[1]> xl1;
  array[Nu[2]] real<lower=ub[2]> xu2;
  array[Nl[2]] real<upper=lb[2]> xl2;

  vector[2] mu;
  vector<lower=0>[2] sigma;
  real<lower=-1,upper=1> a;
  real<lower=-1,upper=1> b;
  real<lower=-1,upper=1> d;
}

model {
  vector[4] Mu = [mu[1], mu[1], mu[2], mu[2]]';
  vector[4] sigmaV = [sigma[1], sigma[1], sigma[2], sigma[2]]';

  matrix[4,4] chol = getChol(a, b, d);

  target += lkj_corr_cholesky_lpdf(chol | 1);

  real sum_;
  for (i in 2:4) {
    for (j in 1:i-1) {
      sum_ = 1;
      for (k in 1:j-1) {
        sum_ -= square(chol[i,k]);
      }
      target += log(sum_) / 2;
    }
  }

  array[N] vector[4] xAll = x;
  for (i in 1:Nu[1]) {
    xAll[(iu1[i] - 1) %/% 2 + 1, (iu1[i] - 1) % 2 + 1] = xu1[i];
  }
  for (i in 1:Nl[1]) {
    xAll[(il1[i] - 1) %/% 2 + 1, (il1[i] - 1) % 2 + 1] = xl1[i];
  }

  for (i in 1:Nu[2]) {
    xAll[(iu2[i] - 1) %/% 2 + 1, (iu2[i] - 1) % 2 + 3] = xu2[i];
  }
  for (i in 1:Nl[2]) {
    xAll[(il2[i] - 1) %/% 2 + 1, (il2[i] - 1) % 2 + 3] = xl2[i];
  }

  target += multi_normal_cholesky_lpdf(xAll | Mu, diag_matrix(sigmaV) * chol);

  target += normal_lpdf(mu | muPriorMean, muPriorSd);
  target += lognormal_lpdf(sigma | sigmaPriorLogMean, 1);
}

generated quantities {
  matrix[4,4] corrMat = multiply_lower_tri_self_transpose(getChol(a, b, d));
}
