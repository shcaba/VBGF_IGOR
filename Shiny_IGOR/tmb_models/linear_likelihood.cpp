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
  DATA_VECTOR(age);
  DATA_VECTOR(len);
  int n = len.size();
  
  // parameters:
  PARAMETER(intercept);
  PARAMETER(slope); 
  PARAMETER(CV_Lt);
  
  // procedures:
  
  Type f = 0.0;
  Type sigma_Lt;
  vector<Type> len_pred(n);
  
  // fit to likelihood
  for (int i = 0; i < n; i++) {
    Type x = age(i);
    if (!isNA(x) && isFinite(x)) {
      len_pred(i) = intercept + slope * x;
      sigma_Lt = CV_Lt * len_pred(i);
      f -= dnorm(len(i), len_pred(i), sigma_Lt, true);
    }
  }
  return f;
}
