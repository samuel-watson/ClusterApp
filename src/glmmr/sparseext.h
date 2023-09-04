#pragma once

#include "general.h"
// extends some of the SparseChol functions and operators to Eigen classes


namespace glmmr {

inline MatrixXd sparse_to_dense(const sparse& m,
                                bool symmetric = true){
  MatrixXd D = MatrixXd::Zero(m.n,m.m);
  for(int i = 0; i < m.n; i++){
    for(int j = m.Ap[i]; j < m.Ap[i+1]; j++){
      D(i,m.Ai[j]) = m.Ax[j];
      if(symmetric) D(m.Ai[j],i) = m.Ax[j];
    }
  }
  return D;
}

inline sparse dense_to_sparse(const MatrixXd& A,
                              bool symmetric = true){
  sparse As(A.rows(),A.cols(),A.data(),true); // this doesn't account for symmetric yet
  return As;
}

inline void mat_mat_mult(const double* a, const double* b, double* ab,
                         const int* Ai, const int* Ap, 
                         int n, int m){
  double val;
  for(int i = 0; i < n; i++){
    for(int j = Ap[i]; j < Ap[i+1]; j++){
      val = a[j];
      for(int k = 0; k<m; k++){
        ab[i+k*n] += val*b[Ai[j]+k*m];
      }
    }
  }

}

inline MatrixXd operator*(const sparse& A, const MatrixXd& B){
  int m = B.cols();
  // alternative code
  // dblvec ab(A.n*m,0);
  // mat_mat_mult(&A.Ax[0],B.data(),&ab[0],&A.Ai[0],&A.Ap[0],n,m);
  MatrixXd AB(A.n,m);
  // AB = Map<MatrixXd>(ab.data(),A.n,m);
  AB.setZero();
  double val;
    for(int i = 0; i < A.n; i++){
      for(int j = A.Ap[i]; j < A.Ap[i+1]; j++){
        val = A.Ax[j];
        for(int k = 0; k<m; k++){
          AB(i,k) += val*B(A.Ai[j],k);
        }
      }
    }
  return AB;
    
}

inline VectorXd operator*(const sparse& A, const VectorXd& B){
  VectorXd AB = VectorXd::Zero(A.n);
  for(int i = 0; i < A.n; i++){
    for(int j = A.Ap[i]; j < A.Ap[i+1]; j++){
      AB(i) += A.Ax[j]*B(A.Ai[j]);
    }
  }
  return AB;
}

inline ArrayXd operator*(const sparse& A, const ArrayXd& B){
  ArrayXd AB = ArrayXd::Zero(A.n);
  for(int i = 0; i < A.n; i++){
    for(int j = A.Ap[i]; j < A.Ap[i+1]; j++){
      AB(i) += A.Ax[j]*B(A.Ai[j]);
    }
  }
  return AB;
}

// multiplication of sparse and diagonal of a vector
inline sparse operator%(const sparse& A, const VectorXd& x){
  sparse Ax(A);
  for(unsigned int i = 0; i < A.Ax.size(); i++){
    Ax.Ax[i] *= x(Ax.Ai[i]);
  }
  return Ax;
}

inline sparse submat_sparse(const sparse& A, intvec rows){
  sparse B;
  B.n = rows.size();
  B.m = A.m;
  for(unsigned int i = 0; i < rows.size(); i++){
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
