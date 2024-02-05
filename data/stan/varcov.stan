functions {
  real sms(real x) {
    return sqrt(1 - x*x);
  }

  matrix getChol(int oneRho, real a, real b, real c, vector f, real g, vector i1) {
    real A = sms(a);
    real B = sms(b);
    real C = sms(c);
    real d = (1-a)*b/(A*B);
    real D = sms(d);
    real e_ = (1-a)*c/(A*C);
    real E = sms(e_);

    real f_;
    if (oneRho) {
      f_ = (a-b*b-B*B*d*d)/(B*B*D*D);
    } else {
      f_ = f[1];
    }

    real F = sms(f_);
    real G = sms(g);
    real h = (1-f_)*g/(F*G);
    real H = sms(h);

    real i_;
    if (oneRho) {
      i_ = (a-c*c-C*C*e_*e_-C*C*E*E*g*g-C*C*E*E*G*G*h*h)/(C*C*E*E*G*G*H*H);
    } else {
      i_ = i1[1];
    }

    real I = sms(i_);

    matrix[6,6] chol = [
      [1, 0, 0, 0, 0, 0],
      [a, A, 0, 0, 0, 0],
      [b, B*d, B*D, 0, 0, 0],
      [b, B*d, B*D*f_, B*D*F, 0, 0],
      [c, C*e_, C*E*g, C*E*G*h, C*E*G*H, 0],
      [c, C*e_, C*E*g, C*E*G*h, C*E*G*H*i_, C*E*G*H*I]
    ];
    return chol;
  }

}

data {
  int<lower=0,upper=1> oneRho;
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

  vector[3] mu;
  vector<lower=0>[3] sigma;
  real<lower=-1,upper=1> a;
  real<lower=-1,upper=1> b;
  real<lower=-1,upper=1> c;
  vector<lower=-1,upper=1>[oneRho ? 0 : 1] f;
  real<lower=-1,upper=1> g;
  vector<lower=-1,upper=1>[oneRho ? 0 : 1] i_;
}

model {
  vector[6] Mu = [mu[1], mu[1], mu[2], mu[2], mu[3], mu[3]]';
  vector[6] sigmaV = [sigma[1], sigma[1], sigma[2], sigma[2], sigma[3], sigma[3]]';

  matrix[6,6] chol = getChol(oneRho, a, b, c, f, g, i_);

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
  matrix[6,6] corrMat = multiply_lower_tri_self_transpose(getChol(oneRho, a, b, c, f, g, i_));
}
