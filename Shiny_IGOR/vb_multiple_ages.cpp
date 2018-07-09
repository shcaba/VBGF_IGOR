#include <TMB.hpp>

template<class Type>

Type objective_function<Type>::operator() () {
  // data:
  DATA_MATRIX(age);
  DATA_VECTOR(len);
  DATA_INTEGER(num_reads);
  int n = len.size();
  
  // parameters:
  PARAMETER(linf); // asymptoptic length
  PARAMETER(kappa); // growth rate
  PARAMETER(t0); // Age at length 0
  
  // procedures:
  
  Type f = 0.0;
  vector<Type> len_pred(n);
  
  for (int i = 0; i < num_reads; i++) {
    Type f_each_read = 0.0;
    for (int j = 0; j < n; j++) {
      len_pred(j) = linf * (1.0 - exp(-kappa * (age(i, j) - t0)));
      f_each_read += (len_pred(j) - len(j)) * (len_pred(j) - len(j));
    }
    f += 0.5 * n * log(f_each_read / n);
  }
  return f;
}
