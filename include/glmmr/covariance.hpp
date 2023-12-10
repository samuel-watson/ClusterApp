#pragma once

#define _USE_MATH_DEFINES

#include "openmpheader.h"
#include "general.h"
#include "algo.h"
#include "interpreter.h"
#include "formula.hpp"
#include "sparse.h"
#include "calculator.hpp"

using namespace Eigen;

namespace glmmr {

class Covariance {
public:
  glmmr::Formula form_;
  const ArrayXXd data_;
  const strvec colnames_;
  dblvec parameters_;
  dblvec other_pars_;
  Covariance(const str& formula,const ArrayXXd &data,const strvec& colnames);
  Covariance(const glmmr::Formula& form,const ArrayXXd &data,const strvec& colnames);
  Covariance(const str& formula,const ArrayXXd &data,const strvec& colnames,const dblvec& parameters);
  Covariance(const glmmr::Formula& form,const ArrayXXd &data,const strvec& colnames,const dblvec& parameters);
  Covariance(const str& formula,const ArrayXXd &data,const strvec& colnames,const ArrayXd& parameters);
  Covariance(const glmmr::Formula& form,const ArrayXXd &data,const strvec& colnames,const ArrayXd& parameters);
  Covariance(const glmmr::Covariance& cov);
  virtual void update_parameters(const dblvec& parameters);
  virtual void update_parameters_extern(const dblvec& parameters);
  virtual void update_parameters(const ArrayXd& parameters);
  virtual int parse();
  virtual double get_val(int b, int i, int j) const;
  virtual MatrixXd Z();
  virtual MatrixXd D(bool chol = false, bool upper = false);
  virtual VectorXd sim_re();
  virtual double log_likelihood(const VectorXd &u);
  virtual double log_determinant();
  virtual int npar() const;
  virtual int B() const;
  virtual int Q() const;
  virtual int max_block_dim() const;
  virtual int block_dim(int b) const;
  virtual void make_sparse();
  virtual MatrixXd ZL();
  virtual MatrixXd LZWZL(const VectorXd& w);
  virtual MatrixXd ZLu(const MatrixXd& u);
  virtual MatrixXd Lu(const MatrixXd& u);
  virtual void set_sparse(bool sparse);
  bool any_group_re() const;
  intvec parameter_fn_index() const;
  virtual intvec re_count() const;
  virtual sparse ZL_sparse();
  virtual sparse Z_sparse() const;
  strvec parameter_names();
  virtual void derivatives(std::vector<MatrixXd>& derivs,int order = 1);
protected:
  std::vector<glmmr::calculator> calc_;
  intvec3d re_pars_;
  //strvec2d fn_;
  std::vector<std::vector<CovFunc> > fn_;
  dblvec2d par_for_calcs_;
  std::vector<MatrixXd> dists;
  intvec re_fn_par_link_;
  intvec re_count_;
  intvec re_order_;
  intvec block_size;
  intvec block_nvar;
  intvec3d re_cols_data_;
  dblvec3d re_temp_data_;
  intvec z_;
  int Q_;
  sparse matZ;
  int n_;
  int B_;
  int npars_;
  MatrixXd dmat_matrix;
  VectorXd zquad;
  bool isSparse = true;
  sparse mat;
  sparse matL;
  SparseChol spchol;
  void update_parameters_in_calculators();
  MatrixXd get_block(int b);
  MatrixXd get_chol_block(int b,bool upper = false);
  MatrixXd D_builder(int b,bool chol = false,bool upper = false);
  void update_ax();
  void L_constructor();
  void Z_constructor();
  MatrixXd D_sparse_builder(bool chol = false,bool upper = false);
};

}

inline glmmr::Covariance::Covariance(const str& formula,
           const ArrayXXd &data,
           const strvec& colnames) :
  form_(formula), data_(data), colnames_(colnames), Q_(parse()), matZ(),
  dmat_matrix(max_block_dim(),max_block_dim()),
  zquad(max_block_dim()) {
  Z_constructor();
};

inline glmmr::Covariance::Covariance(const glmmr::Formula& form,
           const ArrayXXd &data,
           const strvec& colnames) :
  form_(form), data_(data), colnames_(colnames),  
  Q_(parse()),matZ(), dmat_matrix(max_block_dim(),max_block_dim()),
  zquad(max_block_dim()) {
  Z_constructor();
};

inline glmmr::Covariance::Covariance(const str& formula,
           const ArrayXXd &data,
           const strvec& colnames,
           const dblvec& parameters) :
  form_(formula), data_(data), colnames_(colnames), parameters_(parameters), 
  Q_(parse()), matZ(), dmat_matrix(max_block_dim(),max_block_dim()),
  zquad(max_block_dim()), spchol((make_sparse(),mat)) {
  L_constructor();
  Z_constructor();
};

inline glmmr::Covariance::Covariance(const glmmr::Formula& form,
           const ArrayXXd &data,
           const strvec& colnames,
           const dblvec& parameters) :
  form_(form), data_(data), colnames_(colnames), parameters_(parameters), 
  Q_(parse()), matZ(), dmat_matrix(max_block_dim(),max_block_dim()),
  zquad(max_block_dim()), spchol((make_sparse(),mat)) {
  L_constructor();
  Z_constructor();
};

inline glmmr::Covariance::Covariance(const str& formula,
           const ArrayXXd &data,
           const strvec& colnames,
           const ArrayXd& parameters) :
  form_(formula), data_(data), colnames_(colnames),
  parameters_(parameters.data(),parameters.data()+parameters.size()),Q_(parse()), matZ(),
  dmat_matrix(max_block_dim(),max_block_dim()),
  zquad(max_block_dim()), spchol((make_sparse(),mat)) {
  L_constructor();
  Z_constructor();
};

inline glmmr::Covariance::Covariance(const glmmr::Formula& form,
           const ArrayXXd &data,
           const strvec& colnames,
           const ArrayXd& parameters) :
  form_(form), data_(data), colnames_(colnames),
  parameters_(parameters.data(),parameters.data()+parameters.size()),Q_(parse()), matZ(),
  dmat_matrix(max_block_dim(),max_block_dim()),
  zquad(max_block_dim()), spchol((make_sparse(),mat)) {
  L_constructor();
  Z_constructor();
};

inline glmmr::Covariance::Covariance(const glmmr::Covariance& cov) : form_(cov.form_), data_(cov.data_),
colnames_(cov.colnames_),
parameters_(cov.parameters_), Q_(parse()), matZ(),
dmat_matrix(max_block_dim(),max_block_dim()),
zquad(max_block_dim()), spchol((make_sparse(),mat)) {
  L_constructor();
  Z_constructor();
};

inline int glmmr::Covariance::parse(){
  intvec3d re_cols_;
  strvec2d re_par_names_;
  
  // now process each step of the random effect terms
  #ifdef R_BUILD
  if(colnames_.size()!= (unsigned int)data_.cols())Rcpp::stop("colnames length != data columns");
  #endif 
  
  int nre = form_.re_.size();
  
  for(int i = 0; i < nre; i++){
    strvec fn;
    intvec2d fnvars;
    std::stringstream check1(form_.re_[i]);
    std::string intermediate;
    int iter = 0;
    
    while(getline(check1, intermediate, '*')){
      intvec fnvars1;
      fnvars.push_back(fnvars1);
      std::stringstream check2(intermediate);
      std::string intermediate2;
      getline(check2, intermediate2, '(');
      fn.push_back(intermediate2);
      getline(check2, intermediate2, ')');
      if(intermediate2.find(",") != std::string::npos){
        std::stringstream check3(intermediate2);
        std::string intermediate3;
        while(getline(check3, intermediate3, ',')){
          auto colidx = std::find(colnames_.begin(),colnames_.end(),intermediate3);
          if(colidx == colnames_.end()){
            #ifdef R_BUILD
            Rcpp::stop("variable "+intermediate3+" not in data");
            #endif 
          } else {
            int newidx = colidx - colnames_.begin();
            fnvars[iter].push_back(newidx);
          }
        }
      } else {
        auto colidx = std::find(colnames_.begin(),colnames_.end(),intermediate2);
        if(colidx == colnames_.end()){
          #ifdef R_BUILD
          Rcpp::stop("variable "+intermediate2+" not in data");
          #endif
        } else {
          int newidx = colidx - colnames_.begin();
          fnvars[iter].push_back(newidx);
        }
      }
      iter++;
    }
    
    // if any of the functions are group, then use block functions
    auto idxgr = std::find(fn.begin(),fn.end(),"gr");
    dblvec2d groups;
    dblvec vals;
    bool isgr;
    unsigned int j,k,idx,zcol;
    if(form_.z_[i].compare("1")==0){
      zcol = -1;
    } else {
      auto idxz = std::find(colnames_.begin(),colnames_.end(),form_.z_[i]);
      if(idxz == colnames_.end()){
        #ifdef R_BUILD
        Rcpp::stop("z variable "+form_.z_[i]+" not in column names");
        #endif
      } else {
        zcol = idxz - colnames_.begin();
      }
    }
    
    if(idxgr!=fn.end()){
      idx = idxgr - fn.begin();
      isgr = true;
      vals.resize(fnvars[idx].size());
      
      for(j = 0; j < data_.rows(); j++){
        for(k = 0; k < fnvars[idx].size(); k++){
          vals[k] = data_(j,fnvars[idx][k]);
        }
        if(std::find(groups.begin(),groups.end(),vals) == groups.end()){
          groups.push_back(vals);
        }
      }
    } else {
      isgr = false;
      vals.push_back(0.0);
      groups.push_back(vals);
    }
    
    intvec allcols;
    for(j = 0; j< fnvars.size();j++){
      for(k = 0; k < fnvars[j].size();k++){
        allcols.push_back(fnvars[j][k]);
      }
    }
    
    int total_vars = allcols.size();
    int gridx;
    dblvec allvals;
    intvec2d newrecols;
    allvals.resize(total_vars);
    newrecols.resize(fnvars.size());
    int currresize = calc_.size();
    
    calc_.resize(currresize+groups.size());
    re_temp_data_.resize(currresize+groups.size());
    re_cols_.resize(currresize+groups.size());
    re_cols_data_.resize(currresize+groups.size());
    // for each group, create a new z data including the group
    
    int fn_var_counter = 0;
    for(unsigned int m = 0; m < fnvars.size(); m++){
      intvec iter_fn_var_index;
      for(unsigned int p = 0; p < fnvars[m].size(); p++){
        iter_fn_var_index.push_back(p + fn_var_counter);
      }
      fn_var_counter += fnvars[m].size();
      newrecols[m] = iter_fn_var_index;
    }
    
    for(j = 0; j < groups.size(); j++){
      fn_.emplace_back();
      for(const auto& fnvalue: fn){
#ifdef R_BUILD
        if(glmmr::validate_fn(fnvalue))Rcpp::stop("Function " + fnvalue + " not valid");
#endif
        fn_.back().push_back(str_to_covfunc.at(fnvalue));
      }
      //fn_.push_back(fn);
      z_.push_back(zcol);
      re_order_.push_back(form_.re_order_[i]);
      re_cols_[currresize + j] = newrecols;
    }
    
    for(k = 0; k < data_.rows(); k++){
      if(isgr){
        for(int m = 0; m < vals.size(); m++){
          vals[m] = data_(k,fnvars[idx][m]);
        }
        auto gridx2 = std::find(groups.begin(),groups.end(),vals);
        gridx = gridx2 - groups.begin();
      } else {
        gridx = 0;
      }
      for(int m = 0; m<total_vars; m++){
        allvals[m] = data_(k,allcols[m]);
      }
      // this can probably be sped up using an unordered map
      if(std::find(re_temp_data_[gridx + currresize].begin(),
                   re_temp_data_[gridx + currresize].end(),allvals) == re_temp_data_[gridx + currresize].end()){
        re_temp_data_[gridx + currresize].push_back(allvals);
        re_cols_data_[gridx + currresize].push_back(allcols);
      }
    }
  }
  
  // get parameter indexes
  re_pars_.resize(fn_.size());
  re_par_names_.resize(fn_.size());
  bool firsti;
  npars_ = 0;
  int remidx;
  for(int i = 0; i < nre; i++){
    firsti = true;
    for(int j = 0; j < fn_.size(); j++){
      if(re_order_[j]==i){
        if(firsti){
          intvec2d parcount1;
          strvec parnames2;
          str fn_name = "";
          for(int k = 0; k < fn_[j].size(); k++) fn_name += covfunc_to_str.at(fn_[j][k]);
          for(int k = 0; k < fn_[j].size(); k++){
            auto parget = covfunc_to_nvar.find(fn_[j][k]);
            int npars = parget->second;
            intvec parcount2;
            for(int l = 0; l < npars; l++){
              parcount2.push_back(l+npars_);
              parnames2.push_back(fn_name+"."+std::to_string(i)+".("+covfunc_to_str.at(fn_[j][k])+")."+std::to_string(l));
              re_fn_par_link_.push_back(i);
            }
            parcount1.push_back(parcount2);
            npars_ += npars;
          }
          re_pars_[j] = parcount1;
          re_par_names_[j] = parnames2;
          firsti = false;
          remidx = j;
        } else {
          re_pars_[j] = re_pars_[remidx];
          re_par_names_[j] = re_par_names_[remidx];
        }
      }
    }
  }
  
  //now build the reverse polish notation and add distances
  int nvarfn;
  for(int i =0; i<fn_.size();i++){
    std::vector<Do> fn_instruct;
    intvec fn_par;
    int minvalue = 100;
    int ndata = re_temp_data_[i].size();
    MatrixXd distmat = MatrixXd::Zero(ndata*(ndata-1)/2,fn_[i].size());
    dists.push_back(distmat);
    for(unsigned int j = 0; j<fn_[i].size();j++){
      auto min_value_iterator = std::min_element(re_pars_[i][j].begin(),re_pars_[i][j].end());
      if(*min_value_iterator < minvalue) minvalue = *min_value_iterator;
    }
    for(int j = 0; j<fn_[i].size();j++){
      if(fn_[i][j]!=CovFunc::gr){
        nvarfn = re_cols_[i][j].size();
        double dist_val;
        double dist_ab;
        for(int k = 0; k < (re_temp_data_[i].size()-1); k++){
          for(int l = k+1; l < re_temp_data_[i].size(); l++){
            dist_val = 0;
            for(int m = 0; m < nvarfn; m++){
              dist_ab = re_temp_data_[i][k][re_cols_[i][j][m]]-re_temp_data_[i][l][re_cols_[i][j][m]];
              dist_val += dist_ab * dist_ab;
            }
            int idxval = (ndata-1)*k - ((k-1)*k/2) + (l-k-1);
            dists[i](idxval,j) = sqrt(dist_val);
          }
        }
      }
      std::vector<Do> B = glmmr::interpret_re(fn_[i][j]);
      intvec re_par_less_min_ = re_pars_[i][j];
      for(unsigned int k = 0; k < re_pars_[i][j].size(); k++)re_par_less_min_[k] -= minvalue;
      intvec Bpar = glmmr::interpret_re_par(fn_[i][j],j,re_par_less_min_);
      fn_instruct.insert(fn_instruct.end(),B.begin(),B.end());
      fn_par.insert(fn_par.end(),Bpar.begin(),Bpar.end());
    }
    if(fn_[i].size() > 1){
      for(unsigned int j = 0; j < (fn_[i].size()-1); j++){
        fn_instruct.push_back(Do::Multiply);
      }
    }
    calc_[i].instructions = fn_instruct;
    calc_[i].indexes = fn_par;
    calc_[i].parameter_names = re_par_names_[i];
    calc_[i].parameter_count = re_par_names_[i].size();
  }
  //get the number of random effects
  int Qn = 0;
  block_size.resize(calc_.size());
  block_nvar.resize(calc_.size());
  for(unsigned int i = 0; i < calc_.size(); i++){
    Qn += re_temp_data_[i].size();
    block_size[i] = re_temp_data_[i].size();
    block_nvar[i] = re_temp_data_[i][0].size();
  }
  re_count_.resize(form_.re_terms().size());
  std::fill(re_count_.begin(), re_count_.end(), 0);
  for(unsigned int i = 0; i < calc_.size(); i++){
    re_count_[re_order_[i]] += re_temp_data_[i].size();
  }
  B_ = calc_.size();
  n_ = data_.rows();
  return Qn;
}

inline void glmmr::Covariance::Z_constructor(){
  matZ.n = data_.rows();
  matZ.m = Q_;
  matZ.Ap = intvec(data_.rows()+1,0);
  int zcount = 0;
  double insertval;
  for(int i = 0; i < B_; i++){
    dblvec vals(block_nvar[i]);
    dblvec valscomp(block_nvar[i]);
    for(int j = 0; j < block_size[i]; j++){
      for(int m = 0; m < block_nvar[i]; m++){
        valscomp[m] = re_temp_data_[i][j][m];
      }
      for(int k = 0; k < data_.rows(); k++){
        for(int m = 0; m < block_nvar[i]; m++){
          vals[m] = data_(k,re_cols_data_[i][j][m]);
        }
        if(valscomp==vals){
          insertval = z_[i]==-1 ? 1.0 : data_(k,z_[i]);
          matZ.insert(k,zcount,insertval);
        }
      }
      zcount++;
    }
  }
  re_temp_data_.clear();
}


inline void glmmr::Covariance::update_parameters_in_calculators(){
  if(par_for_calcs_.size()==0)par_for_calcs_.resize(B_);
  for(int i = 0; i < B_; i++){
    dblvec par_for_calc;
    for(unsigned int j = 0; j < re_pars_[i].size(); j++){
      for(unsigned int k = 0; k < re_pars_[i][j].size(); k++){
        par_for_calc.push_back(parameters_[re_pars_[i][j][k]]);
      }
    }
    par_for_calcs_[i] = par_for_calc;
  }
}

inline MatrixXd glmmr::Covariance::D(bool chol, bool upper){
  MatrixXd D(Q_,Q_);
  if(isSparse){
    D = D_sparse_builder(chol,upper);
  } else {
    D = D_builder(0,chol,upper);
  }
  return D;
};

inline int glmmr::Covariance::npar() const{
  return npars_;
};


inline int glmmr::Covariance::B() const{
  return B_;
}

inline int glmmr::Covariance::Q() const{
#ifdef R_BUILD
  if(Q_==0)Rcpp::stop("Random effects not initialised");
#endif
  return Q_;
}

inline int glmmr::Covariance::max_block_dim() const{
  int max = 0;
  for(const auto& b: block_size)if(b > max)max = b;
  return max;
}

inline int glmmr::Covariance::block_dim(int b) const{
  return block_size[b];//re_data_[b].rows();
};

inline intvec glmmr::Covariance::parameter_fn_index() const{
  return re_fn_par_link_;
}

inline intvec glmmr::Covariance::re_count() const{
  return re_count_;
}

inline void glmmr::Covariance::update_parameters(const dblvec& parameters){
  if(parameters_.size()==0){
    parameters_ = parameters;
    update_parameters_in_calculators();
    make_sparse();
    spchol.update(mat);
    L_constructor();
  } else {
    parameters_ = parameters;
    update_parameters_in_calculators();
    update_ax();
  }
};

inline void glmmr::Covariance::update_parameters_extern(const dblvec& parameters){
  #ifdef R_BUILD
  if(parameters.size()!=(unsigned)npar())Rcpp::stop(std::to_string(parameters.size())+" covariance parameters provided, "+std::to_string(npar())+" required");
  #endif
  if(parameters_.size()==0){
    parameters_ = parameters;
    update_parameters_in_calculators();
    make_sparse();
    spchol.update(mat);
    L_constructor();
  } else {
    parameters_ = parameters;
    update_parameters_in_calculators();
    update_ax();
  }
};

inline void glmmr::Covariance::update_parameters(const ArrayXd& parameters){
  if(parameters_.size()==0){
    for(unsigned int i = 0; i < parameters.size(); i++){
      parameters_.push_back(parameters(i));
    }
    update_parameters_in_calculators();
  } else if(parameters_.size() == parameters.size()){
    for(unsigned int i = 0; i < parameters.size(); i++){
      parameters_[i] = parameters(i);
    }
    update_parameters_in_calculators();
    update_ax();
  } else {
#ifdef R_BUILD
    Rcpp::stop(std::to_string(parameters.size())+" covariance parameters provided, "+std::to_string(parameters_.size())+" required");
#endif
  }
};

inline double glmmr::Covariance::get_val(int b, int i, int j) const{
  return calc_[b].calculate<CalcDyDx::None>(i,par_for_calcs_[b],dists[b],j,0,0,block_size[b])[0];
}

inline MatrixXd glmmr::Covariance::get_block(int b){
  
#if defined(R_BUILD) && defined(ENABLE_DEBUG)
  if(b > calc_.size()-1)Rcpp::stop("b larger than number of blocks");
  if(parameters_.size()==0)Rcpp::stop("no parameters");
  if(b > B_-1)Rcpp::stop("b is too large");
#endif
  
  int dim = block_dim(b);
  MatrixXd D(dim,dim);
  D.setZero();
  //diagonal
  for(int k = 0; k < dim; k++){
    D(k,k) = get_val(b,k,k);
  }
  if(dim>1){
    for(int i = 0; i < (dim-1); i++){
      for(int j = (i+1); j < dim; j++){
        D(j,i) = get_val(b,j,i);
        D(i,j) = D(j,i);
      }
    }
  }
  return D;
}

inline MatrixXd glmmr::Covariance::Z(){
  return sparse_to_dense(matZ,false,true);
}

inline MatrixXd glmmr::Covariance::get_chol_block(int b,bool upper){
  int n = block_dim(b);//re_data_[b].rows();
  std::vector<double> L(n * n, 0.0);
  
  for (int j = 0; j < n; j++) {
    double s = glmmr::algo::inner_sum(&L[j * n], &L[j * n], j);
    L[j * n + j] = sqrt(get_val(b, j, j) - s);
    for (int i = j + 1; i < n; i++) {
      double s = glmmr::algo::inner_sum(&L[j * n], &L[i * n], j);
      L[i * n + j] = (1.0 / L[j * n + j] * (get_val(b, j, i) - s));
    }
  }
  MatrixXd M = Map<MatrixXd>(L.data(), n, n);
  if (upper) {
    return M;
  } else {
    return M.transpose();
  }
}

inline VectorXd glmmr::Covariance::sim_re(){
  #ifdef R_BUILD
  if(parameters_.size()==0)Rcpp::stop("no parameters, cannot simulate random effects");
  #endif
  VectorXd samps(Q_);
  int idx = 0;
  int ndim;
  boost::variate_generator<boost::mt19937, boost::normal_distribution<> >
    generator(boost::mt19937(time(0)),
              boost::normal_distribution<>());
  if(!isSparse){
    for(int i=0; i< B(); i++){
      MatrixXd L = get_chol_block(i);
      ndim = block_dim(i);//re_data_[i].rows();
      VectorXd zz(ndim);      
      randomGaussian(generator, zz);
      samps.segment(idx,ndim) = L*zz;
      idx += ndim;
    }
  } else {
    VectorXd zz(Q_);
    randomGaussian(generator, zz);
    samps = matL * zz;
  }
  
  return samps;
}

inline MatrixXd glmmr::Covariance::D_builder(int b,
                                             bool chol,
                                             bool upper){
  if (b == B_ - 1) {
    return chol ? get_chol_block(b,upper) : get_block(b);
  }
  else {
    MatrixXd mat1 = chol ? get_chol_block(b,upper) : get_block(b);
    MatrixXd mat2;
    if (b == B_ - 2) {
      mat2 = chol ? get_chol_block(b+1,upper) : get_block(b+1);
    }
    else {
      mat2 = D_builder(b + 1, chol, upper);
    }
    int n1 = mat1.rows();
    int n2 = mat2.rows();
    MatrixXd dmat = MatrixXd::Zero(n1+n2, n1+n2);
    dmat.block(0,0,n1,n1) = mat1;
    dmat.block(n1, n1, n2, n2) = mat2;
    return dmat;
  }
}

inline sparse glmmr::Covariance::ZL_sparse() {
  return matZ * matL;
}

inline sparse glmmr::Covariance::Z_sparse() const{
  return matZ;
}

inline MatrixXd glmmr::Covariance::ZL() {
  sparse ZD = ZL_sparse();
  MatrixXd ZL = sparse_to_dense(ZD,false);
  return ZL;
}

inline MatrixXd glmmr::Covariance::LZWZL(const VectorXd& w){
  sparse ZL = ZL_sparse();
  sparse ZLt = ZL;
  ZLt.transpose();
  ZLt = ZLt % w;
  ZLt *= ZL;
  // add 1 to diagonal
  for(int i = 0; i < ZLt.n; i++){
    for(int j = ZLt.Ap[i]; j<ZLt.Ap[i+1]; j++){
      if(i == ZLt.Ai[j])ZLt.Ax[j]++;
    }
  }
  return sparse_to_dense(ZLt);
}

inline MatrixXd glmmr::Covariance::ZLu(const MatrixXd& u){
  sparse ZL = ZL_sparse();
#if defined(ENABLE_DEBUG) && defined(R_BUILD)
  if(ZL.m != u.rows())Rcpp::stop("ZL*u bad dimension: "+std::to_string(ZL.m)+" cols "+std::to_string(u.rows())+" rows");
#endif
  return ZL * u;
}

inline MatrixXd glmmr::Covariance::Lu(const MatrixXd& u){
  return matL * u;
}

inline double glmmr::Covariance::log_likelihood(const VectorXd &u){
#ifdef R_BUILD
  if(parameters_.size()==0)Rcpp::stop("no covariance parameters, cannot calculate log likelihood");
#endif
  double logdet_val=0.0;
  double loglik_val=0.0;
  int obs_counter=0;
  ArrayXd size_B_array(B_);
  if(!isSparse){
    int blocksize;
    size_B_array.setZero();
    for(int b=0;b<B_;b++){
      blocksize = block_dim(b);
      if(blocksize==1){
        double var = get_val(b,0,0);
        size_B_array[b] = -0.5*log(var) -0.5*log(2*M_PI) -
          0.5*u(obs_counter)*u(obs_counter)/(var);
      } else {
        zquad.setZero();
        dmat_matrix.block(0,0,blocksize,blocksize) = get_chol_block(b);
        logdet_val = 0.0;
        for(int i = 0; i < blocksize; i++){
          logdet_val += 2*log(dmat_matrix(i,i));
        }
        zquad.segment(0,blocksize) = glmmr::algo::forward_sub(dmat_matrix,u.segment(obs_counter,blocksize),blocksize);
        size_B_array[b] = (-0.5*blocksize * log(2*M_PI) - 0.5*logdet_val - 0.5*zquad.transpose()*zquad);
      }
      obs_counter += blocksize;
    }
    loglik_val = size_B_array.sum();
  } else {
    for (auto k : spchol.D) logdet_val += log(k);
    dblvec v(u.data(), u.data()+u.size());
    spchol.ldl_lsolve(&v[0]);
    spchol.ldl_d2solve(&v[0]);
    double quad = glmmr::algo::inner_sum(&v[0],&v[0],Q_);
    loglik_val = (-0.5*Q_ * log(2*M_PI) - 0.5*logdet_val - 0.5*quad);
  }
  return loglik_val;
}

inline double glmmr::Covariance::log_determinant(){
#ifdef R_BUILD
  if(parameters_.size()==0)Rcpp::stop("no covariance parameters, cannot calculate log determinant");
#endif
  int blocksize;
  double logdet_val = 0.0;
  if(!isSparse){
    for(int b=0;b<B_;b++){
      blocksize = block_dim(b);
      dmat_matrix.block(0,0,blocksize,blocksize) = get_chol_block(b);
      for(int i = 0; i < blocksize; i++){
        logdet_val += 2*log(dmat_matrix(i,i));
      }
    }
  } else {
    for (auto k : spchol.D) logdet_val += log(k);
  }
  
  return logdet_val;
}


inline void glmmr::Covariance::make_sparse(){
#ifdef R_BUILD
  if(parameters_.size()==0)Rcpp::stop("no covariance parameters, cannot make sparse");
#endif
  //isSparse = true;
  int dim;
  double val;
  int col_counter=0;
  mat.Ap.clear();
  mat.Ai.clear();
  mat.Ax.clear();
  // algorithm to generate the sparse matrix representation
  for(int b = 0; b < B(); b++){
    dim = block_dim(b);
    for(int i = 0; i < dim; i++){
      mat.Ap.push_back(mat.Ai.size());
      for(int j = 0; j < (i+1); j++){
        val = get_val(b,i,j);
        if(val!=0){
          mat.Ax.push_back(val);
          mat.Ai.push_back((col_counter+j));
        }
      }
    }
    col_counter += dim;
  }
  mat.n = mat.Ap.size();
  mat.m = mat.Ap.size();
  mat.Ap.push_back(mat.Ax.size());
};

inline void glmmr::Covariance::L_constructor(){
  int d = spchol.ldl_numeric();
  (void)d;
  matL = spchol.LD();
}

inline void glmmr::Covariance::set_sparse(bool sparse){
  isSparse = sparse;
}

inline void glmmr::Covariance::update_ax(){
  int llim = 0;
  int nj = 0;
  int ulim = mat.Ap[nj+block_dim(0)];
  int j = 0;
  
  for(int b=0; b < B(); b++){
    for(int i = llim; i<ulim; i++){
      if(i == mat.Ap[j+1])j++;
      mat.Ax[i] = get_val(b,mat.Ai[i]-nj,j-nj);
    }
    llim = ulim;
    if(b<(B()-1)){
      nj += block_dim(b);
      ulim = mat.Ap[nj+block_dim(b+1)];
    }
    if(b == (B()-1)){
      ulim = mat.Ai.size();
    }
  }
  spchol.A_ = mat;
  int d = spchol.ldl_numeric(); // assumes structure of D doesn't change
  (void)d;
  matL = spchol.LD();
};

inline MatrixXd glmmr::Covariance::D_sparse_builder(bool chol,
                                                    bool upper){
  MatrixXd D = MatrixXd::Zero(Q_,Q_);
  if(!chol){
    D = sparse_to_dense(mat,true);
  } else {
    D = sparse_to_dense(matL,false);
  }
  return D;
}

inline bool glmmr::Covariance::any_group_re() const{
  bool gr = false;
  for(unsigned int i = 0; i < fn_.size(); i++){
    for(unsigned int j = 0; j < fn_[i].size(); j++){
      if(fn_[i][j]==CovFunc::gr){
        gr = true;
        break;
      }
    }
    if(gr)break;
  }
  return gr;
}

inline strvec glmmr::Covariance::parameter_names(){
  strvec parnames;
  for(unsigned int i = 0; i < form_.re_.size(); i++){
    for(int j = 0; j < B_; j++){
      if(re_order_[j]==i){
        parnames.insert(parnames.end(),calc_[j].parameter_names.begin(),calc_[j].parameter_names.end());
        break;
      }
    }
  }
  return parnames;
};

inline void glmmr::Covariance::derivatives(std::vector<MatrixXd>& derivs,
                                           int order){
  // get unique parameters
  strvec pars = parameter_names();
  int R = pars.size();
  int matrix_n = order==2 ? R + R*(R+1)/2 + 1 : R+1;
  // initialise all the matrices to zero
  for(int i = 0; i < matrix_n; i++)derivs.push_back(MatrixXd::Zero(Q_,Q_));
  int block_count = 0;
  // iterate over the blocks and insert if the parameter is in the list.
  for(int b = 0; b < B_; b++){
    int block_dimension = block_dim(b);
    int R_block = calc_[b].parameter_names.size();
    intvec par_index;
    for(int k = 0; k < R_block; k++){
      auto par_pos = std::find(pars.begin(),pars.end(),calc_[b].parameter_names[k]);
      int par_pos_int = par_pos - pars.begin();
      par_index.push_back(par_pos_int);
    }
#if defined(R_BUILD) && defined(ENABLE_DEBUG)
    Rcpp::Rcout << "\nSIGMA DERIVATIVES\nBlock " << b << " with " << R_block << " parameters, " << block_dimension << " size";
    Rcpp::Rcout << "\nPar indexes: ";
    for(const auto& i: par_index)Rcpp::Rcout << i << " ";
#endif
    //added conditional parallelisation for large blocks
    dblvec out(matrix_n);
#pragma omp parallel for if(block_dimension > 50) private(out)
    for(int i = 0; i < block_dimension; i++){
      for(int j = i; j < block_dimension; j++){
        if(order == 1){
          out = calc_[b].calculate<CalcDyDx::BetaFirst>(i,par_for_calcs_[b],dists[b],j,0,0,block_dimension);
        } else {
          out = calc_[b].calculate<CalcDyDx::BetaSecond>(i,par_for_calcs_[b],dists[b],j,0,0,block_dimension);
        }
        derivs[0](block_count+i,block_count+j) = out[0];
        if(i!=j)derivs[0](block_count+j,block_count+i) = out[0];
        int index_count = R_block + 1;
        for(int k = 0; k < R_block; k++){
          derivs[par_index[k]+1](block_count+i,block_count+j) = out[k+1];
          if(i!=j)derivs[par_index[k]+1](block_count+j,block_count+i) = out[k+1];
          //second order derivatives
          if(order >= 2){
            for(int l=k; l < R_block; l++){
              int second_pos = par_index[l]*(R-1) - par_index[l]*(par_index[l]-1)/2 + par_index[k];
              derivs[R+1+second_pos](block_count+i,block_count+j) = out[index_count];
              if(i!=j)derivs[R+1+second_pos](block_count+j,block_count+i) = out[index_count];
              index_count++;
            }
          }
        }
      }
    }
    block_count += block_dimension;
  }
}

