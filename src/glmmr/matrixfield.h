#ifndef MATRIXFIELD_H
#define MATRIXFIELD_H

#include <cmath> 
#include <memory>

namespace glmmr {

// create a class to store different sized matrices, not perfect but it'll do!
template<typename T>
class MatrixField{
public: 
  std::vector<std::unique_ptr<T> > data;
  
  MatrixField(int n){
    data.resize(n);
  };
  
  MatrixField(){};
  
  MatrixField(const glmmr::MatrixField<T> &field) {
    for(auto& e : field.data)data.push_back(std::make_unique<T>(*e));
  }
  
  void add(T matrix){
    data.push_back(std::make_unique<T>(matrix));
  }
  
  T operator()(int i)
  {
    return *(data[i]);
  }
  
  Eigen::RowVectorXd get_row(int n, int i){
    return data[n]->row(i);
  }
  
  std::unique_ptr<T> get_ptr(int n){
    return data[n];
  }
  
  void replace(int i, T matrix){
    *(data[i]) = matrix;
  }
  
  int mat_size(int i){
    return data[i]->size();
  }
  
  int size(){
    return data.size();
  }
  
  int rows(int i){
    return data[i]->rows();
  }
  
  int cols(int i){
    return data[i]->cols();
  }
  
  ~MatrixField(){
    data.clear();
  }
  
};

}

#endif
