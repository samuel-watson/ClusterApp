#pragma once

#include "general.h"
// extends some of the SparseChol functions and operators to Eigen classes
// I don't think these are required. This header will be removed in future
// version

namespace glmmr {

inline ArrayXd operator*(const sparse& A, const ArrayXd& B){
  ArrayXd AB = ArrayXd::Zero(A.n);
  for(int i = 0; i < A.n; i++){
    for(int j = A.Ap[i]; j < A.Ap[i+1]; j++){
      AB(i) += A.Ax[j]*B(A.Ai[j]);
    }
  }
  return AB;
}

inline sparse submat_sparse(const sparse& A, intvec rows){
  sparse B;
  B.n = rows.size();
  B.m = A.m;
  for(int i = 0; i < rows.size(); i++){
    B.Ap.push_back(B.Ai.size());
    for(int j = A.Ap[rows[i]]; j < A.Ap[rows[i]+1]; j++){
      B.Ai.push_back(A.Ai[j]);
      B.Ax.push_back(A.Ax[j]);
    }
  }
  B.Ap.push_back(B.Ax.size());
  return B;
}

}

