functions {
  real sms(real x) {
    return sqrt(1 - x*x);
  }

  matrix getChol(real a, real b, real c, real f, real g, real i_) {
    real A = sms(a);
    real B = sms(b);
    real C = sms(c);
    real d = (1-a)*b/(A*B);
    real D = sms(d);
    real e_ = (1-a)*c/(A*C);
    real E = sms(e_);
    real F = sms(f);
    real G = sms(g);
    real h = (1-f)*g/(F*G);
    real H = sms(h);
    real I = sms(i_);

    matrix[6,6] chol = [
      [1, 0, 0, 0, 0, 0],
      [a, A, 0, 0, 0, 0],
      [b, B*d, B*D, 0, 0, 0],
      [b, B*d, B*D*f, B*D*F, 0, 0],
      [c, C*e_, C*E*g, C*E*G*h, C*E*G*H, 0],
      [c, C*e_, C*E*g, C*E*G*h, C*E*G*H*i_, C*E*G*H*I]
    ];
    return chol;
  }

}

data {
  real muOffset;
  real lb;
  real ub;

  int<lower=0> N;
  array[N] vector[6] x;

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
  real<lower=-1,upper=1> b;
  real<lower=-1,upper=1> c;
  real<lower=-1,upper=1> f;
  real<lower=-1,upper=1> g;
  real<lower=-1,upper=1> i_;
}

model {
  vector[6] Mu = [mu, mu, mu - muOffset, mu - muOffset, mu + muOffset, mu + muOffset]';
  vector[6] sigmaV = rep_vector(sigma, 6);

  matrix[6,6] chol = getChol(a, b, c, f, g, i_);

  target += lkj_corr_cholesky_lpdf(chol | 1);

  real sum_;
  for (i in 2:6) {
    for (j in 1:i-1) {
      sum_ = 1;
      for (k in 1:j-1) {
        sum_ -= square(chol[i,k]);
      }
      target += log(sum_) / 2;
    }
  }

  array[N] vector[6] xAll = x;
  for (i in 1:Nu) {
    xAll[(iu[i] - 1) %/% 6 + 1, (iu[i] - 1) % 6 + 1] = xu[i];
  }
  for (i in 1:Nl) {
    xAll[(il[i] - 1) %/% 6 + 1, (il[i] - 1) % 6 + 1] = xl[i];
  }

  target += multi_normal_cholesky_lpdf(xAll | Mu, diag_matrix(sigmaV) * chol);

  target += normal_lpdf(mu | muPriorMean, muPriorSd);
  target += lognormal_lpdf(sigma | sigmaPriorLogMean, 1);
}

generated quantities {
  matrix[6,6] corrMat = multiply_lower_tri_self_transpose(getChol(a, b, c, f, g, i_));
}
