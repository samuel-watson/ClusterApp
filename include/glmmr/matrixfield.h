#pragma once

#include "general.h"

namespace glmmr {

// create a simple class to store different sized matrices of arbitrary number but trying to minimise 
// copying with std vector as the vector grows.

template<typename T>
class MatrixField{
public: 
  std::vector<std::unique_ptr<T> > data;
  
  MatrixField(int n){
    data.resize(n);
  };
  
  MatrixField(){};
  
  MatrixField(const glmmr::MatrixField<T> &field) {
    for(const auto& e : field.data)data.push_back(std::make_unique<T>(*e));
  }
  
  void add(T matrix){
    data.push_back(std::make_unique<T>(matrix));
  }

  template<class... Args>
  void add_new(Args&&... args){
    data.push_back(std::unique_ptr<T>(new T(std::forward<Args>(args)...)));
  }
  
  T operator()(int i) const
  {
    if(i >= data.size()) throw std::runtime_error("Accessing index out of range matrix field");
    return *(data[i]);
  }
  
  void sum(int i, const Eigen::MatrixXd& A)
  {
    *(data[i]) += A;
  }
  
  Eigen::RowVectorXd get_row(int n, int i) const{
    return data[n]->row(i);
  }
  
  std::unique_ptr<T> get_ptr(int n) const{
    return data[n];
  }
  
  void replace(int i, T matrix){
    *(data[i]) = matrix;
  }
  
  int mat_size(int i) const{
    return data[i]->size();
  }
  
  int size() const{
    return data.size();
  }
  
  int rows(int i) const{
    return data[i]->rows();
  }
  
  int cols(int i) const{
    return data[i]->cols();
  }
  
  ~MatrixField(){
    data.clear();
  }
  
};

}

