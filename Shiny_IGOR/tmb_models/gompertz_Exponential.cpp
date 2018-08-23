#include <TMB.hpp>

template<class Type>
bool isNA(Type x){
  return R_IsNA(asDouble(x));
}

template<class Type>
bool isFinite(Type x){
  return R_finite(asDouble(x));
}

template<class Type>
Type objective_function<Type>::operator() () {
  // data:
  DATA_MATRIX(age);
  DATA_VECTOR(len);
  DATA_SCALAR(CV_e);
  DATA_INTEGER(num_reads);
  
  // parameters:
  PARAMETER(a); // upper asymptote
  PARAMETER(b); // growth range
  PARAMETER(k); // growth rate
  PARAMETER(CV_Lt);
  
  PARAMETER(beta);
  
  PARAMETER_VECTOR(age_re);
  
  // procedures:
  Type n = len.size();
  
  Type nll = 0.0; // Initialize negative log-likelihood
  
  Type eps = 1e-5;
  
  
  CV_e = CV_e < 0.05 ? 0.05 : CV_e;
  
  for (int i = 0; i < n; i++) {
    Type x = age_re(i);
    if (!isNA(x) && isFinite(x)) {
      Type len_pred = a * exp(-b * exp(-k * x));
      
      Type sigma_e = CV_e * x + eps;
      Type sigma_Lt = CV_Lt * (len_pred + eps);
      
      nll -= dnorm(len(i), len_pred, sigma_Lt, true);
      nll -= dexp(x, beta, true);
      
      for (int j = 0; j < num_reads; j++) {
        if (!isNA(age(j, i)) && isFinite(age(j, i)) && age(j, i) >= 0) {
          nll -= dnorm(age(j, i), x, sigma_e, true); 
        }
      } 
    }
  }
  
  return nll;
}
