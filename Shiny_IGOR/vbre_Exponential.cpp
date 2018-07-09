#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() () {
  // data:
  DATA_MATRIX(age);
  DATA_VECTOR(len);
  DATA_SCALAR(CV_e);
  DATA_INTEGER(num_reads);
  
  // parameters:
  PARAMETER(linf); // asymptoptic length
  PARAMETER(kappa); // growth rate
  PARAMETER(t0); // Age at length 0
  
  
  PARAMETER(CV_Lt);
  
  PARAMETER(beta);
  
  PARAMETER_VECTOR(age_re);
  
  // procedures:
  Type n = len.size();
  
  Type nll = 0.0; // Initialize negative log-likelihood
  
  Type eps = 1e-5;
  
  if (linf < max(len)) {
    linf *= 1.1;
  }
  
  CV_e = CV_e < 0.05 ? 0.05 : CV_e;
  
  for (int i = 0; i < n; i++) {
    Type len_pred = linf * (1.0 - exp(-kappa * (age_re(i) - t0)));  
    
    Type sigma_e = CV_e * age_re(i) + eps;
    Type sigma_Lt = CV_Lt * (len_pred + eps);
    
    nll -= dnorm(len(i), len_pred, sigma_Lt, true);
    nll -= dexp(age_re(i), beta, true);
    
    for (int j = 0; j < num_reads; j++) {
      if (age(j, i) >= 0) {
        nll -= dnorm(age(j, i), age_re(i), sigma_e, true); 
      }
    }
  }
  
  return nll;
}
