#pragma once

#include "general.h"

namespace glmmr {
namespace algo {
inline double inner_sum(double* li, double* lj, int n)
{
  double s = 0;
  for (int i = 0; i < n; i++) {
    s += li[i] * lj[i];
  }
  return s;
}

inline Eigen::VectorXd forward_sub(const Eigen::MatrixXd& U,
                                   const Eigen::VectorXd& u,
                                   const int& n)
{
  Eigen::VectorXd y(n);
  for (int i = 0; i < n; i++) {
    double lsum = 0;
    for (int j = 0; j < i; j++) {
      lsum += U(i,j) * y(j);
    }
    y(i) = (u(i) - lsum) / U(i,i);
  }
  return y;
}

template<typename T>
inline void combinations(const std::vector<std::vector<T> >& vecs, 
                         unsigned int n,
                         unsigned int m,
                         std::vector<T>& buffer,
                         std::vector<std::vector<T> >& result){
  buffer[n] = vecs[n][m];
  if(n == vecs.size()-1){
    result.push_back(buffer);
  } else {
    for(unsigned int i = 0; i < vecs[n+1].size(); i++){
      combinations(vecs,n+1,i,buffer,result);
    }
  }
}

template<typename T>
inline T prod_vec(std::vector<T> vec){
  T result = 1;
  for(const auto& val: vec)result *= val;
  return result;
}

}

namespace Eigen_ext {

template<class ArgType, class RowIndexType, class ColIndexType>
class index_functor {
  const ArgType &m_arg;
  const RowIndexType &m_rowIndices;
  const ColIndexType &m_colIndices;
public:
  typedef Eigen::Matrix<typename ArgType::Scalar,
                        RowIndexType::SizeAtCompileTime,
                        ColIndexType::SizeAtCompileTime,
                        ArgType::Flags&Eigen::RowMajorBit?Eigen::RowMajor:Eigen::ColMajor,
                        RowIndexType::MaxSizeAtCompileTime,
                        ColIndexType::MaxSizeAtCompileTime> MatrixType;
  
  index_functor(const ArgType& arg, const RowIndexType& row_indices, const ColIndexType& col_indices)
    : m_arg(arg), m_rowIndices(row_indices), m_colIndices(col_indices)
  {}
  
  const typename ArgType::Scalar& operator() (Eigen::Index row, Eigen::Index col) const {
    return m_arg(m_rowIndices[row], m_colIndices[col]);
  }
};

template <class ArgType, class RowIndexType, class ColIndexType>
Eigen::CwiseNullaryOp<index_functor<ArgType,RowIndexType,ColIndexType>, typename index_functor<ArgType,RowIndexType,ColIndexType>::MatrixType>
submat(const Eigen::MatrixBase<ArgType>& arg, const RowIndexType& row_indices, const ColIndexType& col_indices)
{
  typedef index_functor<ArgType,RowIndexType,ColIndexType> Func;
  typedef typename Func::MatrixType MatrixType;
  return MatrixType::NullaryExpr(row_indices.size(), col_indices.size(), Func(arg.derived(), row_indices, col_indices));
}

inline void removeRow(Eigen::MatrixXd& matrix, 
               unsigned int rowToRemove)
{
  unsigned int numRows = matrix.rows()-1;
  unsigned int numCols = matrix.cols();
  
  if(rowToRemove < numRows)
    matrix.block(rowToRemove,0,numRows-rowToRemove,numCols) = 
      matrix.bottomRows(numRows-rowToRemove);
  
  matrix.conservativeResize(numRows,numCols);
}

inline void removeColumn(Eigen::MatrixXd& matrix, 
                  unsigned int colToRemove)
{
  unsigned int numRows = matrix.rows();
  unsigned int numCols = matrix.cols()-1;
  
  if( colToRemove < numCols )
    matrix.block(0,colToRemove,numRows,numCols-colToRemove) = 
      matrix.rightCols(numCols-colToRemove);
  
  matrix.conservativeResize(numRows,numCols);
}

inline void removeElement(Eigen::VectorXd& matrix, 
                  unsigned int elemToRemove)
{
  unsigned int nSize = matrix.size()-1;
  
  if( elemToRemove < nSize )
    matrix.segment(elemToRemove,nSize-elemToRemove) = 
      matrix.tail(nSize-elemToRemove);
  
  matrix.conservativeResize(nSize);
}

inline bool issympd(Eigen::MatrixXd& mat){
  Eigen::LLT<Eigen::MatrixXd> lltOfA(mat);
  return lltOfA.info() == Eigen::NumericalIssue;
}
}

class SigmaBlock {
  public:
    intvec Dblocks;
    intvec RowIndexes;
    SigmaBlock(){};
    SigmaBlock(const intvec& db) : Dblocks(db) {};
    SigmaBlock(const SigmaBlock& x): Dblocks(x.Dblocks), RowIndexes(x.RowIndexes) {};
    
    bool operator==(const intvec& x){
      bool element_is_in = false;
      for(auto i : x){
        auto it = std::find(Dblocks.begin(),Dblocks.end(),i);
        if(it != Dblocks.end()){
          element_is_in = true;
          break;
        } 
      }
      return element_is_in;
    }
    
    void add(const intvec& x){
      intvec xout;
      bool element_is_in = false;
      for(auto i : x){
        auto it = std::find(Dblocks.begin(),Dblocks.end(),i);
        if(it != Dblocks.end()){
          element_is_in = true;
        } else {
          xout.push_back(i); 
        }
      }
      if(element_is_in){
        Dblocks.insert(Dblocks.end(),xout.begin(),xout.end());
        std::sort(Dblocks.begin(),Dblocks.end());
      }
    }
    
    void merge(const SigmaBlock& x){
      RowIndexes.insert(RowIndexes.end(),x.RowIndexes.begin(),x.RowIndexes.end());
      std::sort(RowIndexes.begin(), RowIndexes.end() );
      RowIndexes.erase( std::unique( RowIndexes.begin(), RowIndexes.end() ), RowIndexes.end() );
      Dblocks.insert(Dblocks.end(),x.Dblocks.begin(),x.Dblocks.end());
      std::sort(Dblocks.begin(), Dblocks.end() );
      Dblocks.erase( std::unique( Dblocks.begin(), Dblocks.end() ), Dblocks.end() );
    }
    
    void add_row(int i){
      RowIndexes.push_back(i);
    }
    
    intvec rows(){
      return RowIndexes;
    }
};

}
