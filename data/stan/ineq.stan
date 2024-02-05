functions {
  real area(real mu, real sigma, real x1, real x2) {
    if (mu > (x1 + x2) / 2) {
      return log_diff_exp(normal_lcdf(x2 | mu, sigma), normal_lcdf(x1 | mu, sigma));
    } else {
      return log_diff_exp(normal_lcdf(mu * 2 - x1 | mu, sigma), normal_lcdf(mu * 2 - x2 | mu, sigma));
    }
  }

  real likelihood(real s, real l, real k, real b, real x) {
    real mu1 = (l+k)/(1-k);
    real sigma1 = 1/sqrt(2*b*7*(1-k));
    real mu2 = (l-k)/(1+k);
    real sigma2 = 1/sqrt(2*b*7*(1+k));
    real xc = sqrt(3) - 1;
    real ratio;
    real area1;
    real area2;
    real i;

    if (s == 1) {
      ratio = normal_lpdf(xc | mu1, sigma1) - normal_lpdf(xc | mu2, sigma2);
      area1 = area(mu1, sigma1, -2, xc);
      area2 = area(mu2, sigma2, xc, 2);
      if (x <= xc) {
        return normal_lpdf(x | mu1, sigma1) - log_sum_exp(area1, area2 + ratio);
      } else {
        return normal_lpdf(x | mu2, sigma2) - log_sum_exp(area1 - ratio, area2);
      }
    }
    else if (s == 2) {
      return normal_lpdf(x | mu1, sigma1) - area(mu1, sigma1, -2, 2);
    } else {
      return normal_lpdf(x | mu2, sigma2) - area(mu2, sigma2, -2, 2);
    }
  }
}

data {
  int<lower=0> Np; // # of participants
  int<lower=0> Nt; // # of targets
  int<lower=0> N; // # of data points
  vector[N] x; // slider position
  array[N] int<lower=1,upper=Np> p; // participant ID
  array[N] int<lower=1,upper=Nt> t; // target
  array[N] int<lower=1,upper=3> s; // slider, 1 = balanced, 2 = selfMore, 3 = targetMore
}

parameters {
  vector[Nt] muL;
  array[Np] vector[Nt] lambda;

  vector<lower=0,upper=.95>[Np] kappa;

  real<lower=0> beta_;
}

model {
  muL ~ normal(0, 1);
  for (i in 1:Np) {
    lambda[i] ~ normal(muL, .5);
  }

  beta_ ~ lognormal(-1, .5);

  for (i in 1:N) {
    target += likelihood(s[i], lambda[p[i]][t[i]], kappa[p[i]], beta_, x[i]);
  }
}

generated quantities {
  vector[Np] llhd = rep_vector(0, Np);

  for (i in 1:N) {
    llhd[p[i]] += likelihood(s[i], lambda[p[i]][t[i]], kappa[p[i]], beta_, x[i]);
  }
}
