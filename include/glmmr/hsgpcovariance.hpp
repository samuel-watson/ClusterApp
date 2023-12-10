#pragma once

#include "covariance.hpp"

namespace glmmr {

using namespace Eigen;

class hsgpCovariance : public Covariance {
public:
  int dim;
  intvec m;
  ArrayXXd hsgp_data;
  ArrayXd L_boundary;
  hsgpCovariance(const std::string& formula,const ArrayXXd& data,const strvec& colnames);
  hsgpCovariance(const glmmr::Formula& formula,const ArrayXXd& data,const strvec& colnames);
  hsgpCovariance(const std::string& formula,const ArrayXXd& data,const strvec& colnames,const dblvec& parameters);
  hsgpCovariance(const glmmr::Formula& formula,const ArrayXXd& data,const strvec& colnames,const dblvec& parameters);
  hsgpCovariance(const glmmr::hsgpCovariance& cov);
  double spd_nD(int i);
  ArrayXd phi_nD(int i);
  MatrixXd ZL() override;
  MatrixXd D(bool chol = true, bool upper = false) override;
  MatrixXd LZWZL(const VectorXd& w) override;
  MatrixXd ZLu(const MatrixXd& u) override;
  MatrixXd Lu(const MatrixXd& u) override;
  VectorXd sim_re() override;
  sparse ZL_sparse() override;
  int Q() const override;
  double log_likelihood(const VectorXd &u) override;
  double log_determinant() override;
  void update_parameters(const dblvec& parameters) override;
  void update_parameters(const ArrayXd& parameters) override;
  void update_parameters_extern(const dblvec& parameters) override;
  void set_function(bool squared_exp);
  MatrixXd PhiSPD(bool lambda = true, bool inverse = false);
  ArrayXd LambdaSPD();
  void update_approx_parameters(intvec m_, ArrayXd L_);
  void update_approx_parameters();
protected:
  int total_m;
  MatrixXd L; // Half-eigen decomposition of Lambda + PhiTPhi m^2 * m^2
  ArrayXd Lambda;
  ArrayXXi indices;
  MatrixXd Phi;
  MatrixXd PhiT;
  bool sq_exp = false;
  void parse_hsgp_data();
  void gen_indices();
  void gen_phi_prod();
  void update_lambda();
};

}


inline glmmr::hsgpCovariance::hsgpCovariance(const std::string& formula,
               const ArrayXXd& data,
               const strvec& colnames) : Covariance(formula, data, colnames),
               dim(this->re_cols_data_[0][0].size()),
               m(dim),
               hsgp_data(data.rows(),dim),
               L_boundary(dim),
               L(data.rows(),1), 
               Lambda(1), 
               indices(1,dim), 
               Phi(data.rows(),1), 
               PhiT(2,2) {
  isSparse = false;
  for(int i = 0; i < dim; i++)L_boundary(i) = 1.5;
  std::fill(m.begin(),m.end(),10);
  parse_hsgp_data();
  update_approx_parameters();
};

inline glmmr::hsgpCovariance::hsgpCovariance(const glmmr::Formula& formula,
               const ArrayXXd& data,
               const strvec& colnames) : Covariance(formula, data, colnames),
               dim(this->re_cols_data_[0][0].size()),
               m(dim),
               hsgp_data(data.rows(),dim),
               L_boundary(dim),
               L(data.rows(),1), 
               Lambda(1), 
               indices(1,dim), 
               Phi(data.rows(),1), 
               PhiT(2,2) {
  isSparse = false;
  for(int i = 0; i < dim; i++)L_boundary(i) = 1.5;
  std::fill(m.begin(),m.end(),10);
  parse_hsgp_data();
  update_approx_parameters();
};

inline glmmr::hsgpCovariance::hsgpCovariance(const std::string& formula,
               const ArrayXXd& data,
               const strvec& colnames,
               const dblvec& parameters) : Covariance(formula, data, colnames, parameters),
               dim(this->re_cols_data_[0][0].size()),
               m(dim),
               hsgp_data(data.rows(),dim),
               L_boundary(dim),
               L(data.rows(),1), 
               Lambda(1),
               indices(1,dim), 
               Phi(data.rows(),1), 
               PhiT(2,2) {
  isSparse = false;
  for(int i = 0; i < dim; i++)L_boundary(i) = 1.5;
  std::fill(m.begin(),m.end(),10);
  parse_hsgp_data();
  update_approx_parameters();
  update_lambda();
};

inline glmmr::hsgpCovariance::hsgpCovariance(const glmmr::Formula& formula,
               const ArrayXXd& data,
               const strvec& colnames,
               const dblvec& parameters) : Covariance(formula, data, colnames, parameters),
               dim(this->re_cols_data_[0][0].size()),
               m(dim),
               hsgp_data(data.rows(),dim),
               L_boundary(dim),
               L(data.rows(),1), 
               Lambda(1),
               indices(1,dim), 
               Phi(data.rows(),1), 
               PhiT(2,2) {
  isSparse = false;
  for(int i = 0; i < dim; i++)L_boundary(i) = 1.5;
  std::fill(m.begin(),m.end(),10);
  parse_hsgp_data();
  update_approx_parameters();
  update_lambda();
};

inline glmmr::hsgpCovariance::hsgpCovariance(const glmmr::hsgpCovariance& cov) : Covariance(cov.form_, cov.data_, cov.colnames_, cov.parameters_), 
dim(cov.dim),m(cov.m), hsgp_data(cov.hsgp_data),
L_boundary(cov.L_boundary), L(cov.L), Lambda(cov.Lambda), 
indices(cov.indices), Phi(cov.Phi), PhiT(cov.PhiT) {
  isSparse = false;
};

inline void glmmr::hsgpCovariance::parse_hsgp_data(){
  for(int i = 0; i < dim; i++){
    hsgp_data.col(i) = this->data_.col(this->re_cols_data_[0][0][i]);
  }
  auto sqexpidx = std::find(this->fn_[0].begin(),this->fn_[0].end(),CovFunc::sqexp);
  if(!(sqexpidx == this->fn_[0].end())){
    sq_exp = true;
  } else {
    auto expidx = std::find(this->fn_[0].begin(),this->fn_[0].end(),CovFunc::fexp);
    if(!(expidx == this->fn_[0].end())){
      sq_exp = false;
    } else {
      #ifdef R_BUILD
      Rcpp::stop("HSGP only allows exp and sqexp currently.");
      #endif
    }
  }
}

inline void glmmr::hsgpCovariance::update_approx_parameters(intvec m_, ArrayXd L_){
  m = m_;
  L_boundary = L_;
  total_m = glmmr::algo::prod_vec(m);
  this->Q_ = total_m;
  indices.conservativeResize(total_m,NoChange);
  Phi.conservativeResize(NoChange,total_m);
  PhiT.conservativeResize(total_m,total_m);
  Lambda.conservativeResize(total_m);
  L.conservativeResize(NoChange,total_m);
  gen_indices();
  gen_phi_prod();
}

inline void glmmr::hsgpCovariance::update_approx_parameters(){
  total_m = glmmr::algo::prod_vec(m);
  this->Q_ = total_m;
  indices.conservativeResize(total_m,NoChange);
  Phi.conservativeResize(NoChange,total_m);
  PhiT.conservativeResize(total_m,total_m);
  Lambda.conservativeResize(total_m);
  L.conservativeResize(NoChange,total_m);
  gen_indices();
  gen_phi_prod();
}

inline double glmmr::hsgpCovariance::spd_nD(int i){
  double wprod = 0;
  for(int d = 0; d < dim; d++) {
    double w = (indices(i,d)*M_PI)/(2*L_boundary(d));
    wprod += w*w;
  }
  double S;
  double phisq = parameters_[1] * parameters_[1];
  if(sq_exp){
    S = parameters_[0] * pow(2 * M_PI, dim/2.0) * pow(parameters_[1],dim) * exp(-0.5 * phisq * wprod);
  } else {
    double S1 = parameters_[0] * pow(4 * M_PI, dim/2.0) * (tgamma(0.5*(dim+1))/(tgamma(0.5)*parameters_[1]));
    double S2 = 1/phisq + wprod;
    S = S1 * pow(S2,-1*(dim+1)/2.0);
  }
  return S;
}

inline ArrayXd glmmr::hsgpCovariance::phi_nD(int i){
  ArrayXd fi1(hsgp_data.rows());
  ArrayXd fi2(hsgp_data.rows());
  fi1 = (1/sqrt(L_boundary(0))) * sin(indices(i,0)*M_PI*(hsgp_data.col(0)+L_boundary(0))/(2*L_boundary(0)));
  if(dim > 1){
    for(int d = 1; d < dim; d++){
      fi2 = (1/sqrt(L_boundary(d))) * sin(indices(i,d)*M_PI*(hsgp_data.col(d)+L_boundary(d))/(2*L_boundary(d)));
      fi1 *= fi2;
    }
  }
  return fi1;
}

inline MatrixXd glmmr::hsgpCovariance::D(bool chol, bool upper){
  MatrixXd As = glmmr::hsgpCovariance::ZL();
  if(chol){
    if(upper){
      return As.transpose();
    } else {
      return As;
    }
  } else {
    return As * As.transpose();
  }
}

inline VectorXd glmmr::hsgpCovariance::sim_re(){
  #ifdef R_BUILD
  if(parameters_.size()==0)Rcpp::stop("no parameters");
  #endif
  VectorXd samps(this->Q_);
  boost::variate_generator<boost::mt19937, boost::normal_distribution<> >
    generator(boost::mt19937(time(0)),
              boost::normal_distribution<>());
  VectorXd zz(this->Q_);      
  randomGaussian(generator, zz);
  samps = glmmr::hsgpCovariance::ZL()*zz;
  return samps;
}

inline void glmmr::hsgpCovariance::update_parameters_extern(const dblvec& parameters){
  parameters_ = parameters;
  update_lambda();
};

inline void glmmr::hsgpCovariance::update_parameters(const dblvec& parameters){
  parameters_ = parameters;
  update_lambda();
};

inline void glmmr::hsgpCovariance::update_parameters(const ArrayXd& parameters){
  if(parameters_.size()==0){
    for(unsigned int i = 0; i < parameters.size(); i++){
      parameters_.push_back(parameters(i));
    }
  } else {
    for(unsigned int i = 0; i < parameters.size(); i++){
      parameters_[i] = parameters(i);
    }
  }
  update_lambda();
};

inline MatrixXd glmmr::hsgpCovariance::ZL(){
  MatrixXd ZL = PhiSPD();
  return ZL;
}

inline MatrixXd glmmr::hsgpCovariance::LZWZL(const VectorXd& w){
  MatrixXd ZL = glmmr::hsgpCovariance::ZL();
  MatrixXd LZWZL = ZL.transpose() * w.asDiagonal() * ZL;
  LZWZL += MatrixXd::Identity(LZWZL.rows(), LZWZL.cols());
  return LZWZL;
}

inline MatrixXd glmmr::hsgpCovariance::ZLu(const MatrixXd& u){
  MatrixXd ZLu = glmmr::hsgpCovariance::ZL() * u;
  return ZLu;
}

inline MatrixXd glmmr::hsgpCovariance::Lu(const MatrixXd& u){
  MatrixXd ZLu = glmmr::hsgpCovariance::ZL() * u;
  return ZLu;
}

inline sparse glmmr::hsgpCovariance::ZL_sparse(){
  sparse dummy;
  return dummy;
}

inline int glmmr::hsgpCovariance::Q() const{
  return this->Q_;
}

inline double glmmr::hsgpCovariance::log_likelihood(const VectorXd &u){
  double ll = 0;
  double logdet = log_determinant();
  VectorXd uquad = u * L;
  ll += (-0.5*hsgp_data.rows() * log(2*M_PI) - 0.5*uquad.transpose()*uquad);
  ll += 0.5*logdet;
  return -1.0*ll;
}

inline double glmmr::hsgpCovariance::log_determinant(){
  double logdet = 0;
  for(int i = 0; i < indices.rows(); i++){
    logdet += log(Lambda(i));
  }
  return logdet;
}

inline void glmmr::hsgpCovariance::set_function(bool squared_exp){
  sq_exp = squared_exp;
}

inline void glmmr::hsgpCovariance::gen_indices(){
  intvec2d indices_vec;
  intvec ind_buffer(dim);
  intvec2d linspace_vec;
  for(int i = 0; i < dim; i++){
    intvec linspaced(m[i]);
    for(int k = 1; k <= m[i]; k++)linspaced[k-1] = k;
    linspace_vec.push_back(linspaced);
  }
  for(unsigned int i = 0; i < linspace_vec[0].size(); i++){
    glmmr::algo::combinations(linspace_vec,0,i,ind_buffer,indices_vec);
  }
  // copy into array
  for(int i = 0; i < indices_vec.size(); i++){
    for(int j = 0; j < indices_vec[0].size(); j++){
      indices(i,j) = indices_vec[i][j];
    }
  }
}

inline void glmmr::hsgpCovariance::gen_phi_prod(){
  for(int i = 0; i < total_m; i++){
    ArrayXd phi = phi_nD(i);
    Phi.col(i) = phi.matrix();
  }
  PhiT = Phi.transpose() * Phi;
}

inline MatrixXd glmmr::hsgpCovariance::PhiSPD(bool lambda, bool inverse){
  MatrixXd pnew = Phi;
  if(lambda){
    if(!inverse){
      pnew *= Lambda.sqrt().matrix().asDiagonal();
    } else {
      pnew *= Lambda.sqrt().inverse().matrix().asDiagonal();
    }
  }
  return pnew;
}

inline ArrayXd glmmr::hsgpCovariance::LambdaSPD(){
  return Lambda;
}

inline void glmmr::hsgpCovariance::update_lambda(){
  for(int i = 0; i < total_m; i++){
    Lambda(i) = spd_nD(i);
  }
  L = PhiSPD(true,true);
}
