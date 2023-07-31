#ifndef RBOBYQA_H
#define RBOBYQA_H

#include <algorithm>
#include <vector>
#include "bobyqa.h"
#include "functor.h"

namespace rminqa {

template<typename Derived, typename Partype>
class Rbobyqa {
public:
  Partype lower() const { return lower_; }
  Partype upper() const { return upper_; }
  Partype par() const { return par_; }
  double fval() const { return fval_; }
  int feval() const { return feval_; }
  std::string msg() const { return msg_; }
private:
  Partype lower_, upper_;
  Partype par_;

  double fval_;
  int feval_;

  std::string msg_;
  void Update_msg(int res) {
    switch(res) {
    case 0:
      msg_ = "Normal exit from bobyqa";
      break;
    case -1:
      msg_ = "bobyqa -- NPT is not in the required interval";
      break;
    case -2:
      msg_ = "bobyqa -- one of the box constraint ranges is too small (< 2*RHOBEG)";
      break;
    case -3:
      msg_ = "bobyqa detected too much cancellation in denominator";
      break;
    case -4:
      msg_ = "bobyqa -- maximum number of function evaluations exceeded";
      break;
    case -5:
      msg_ = "bobyqa -- a trust region step failed to reduce q";
      break;
    default: ;
    }
  }

public:
  struct RbobyqaControl {
    int npt = 0;
    double rhobeg = 0.0;
    double rhoend = 0.0;
    int iprint = 0;
    int maxfun = 0;
  } control;

  Rbobyqa() {}

  void set_lower(const Partype &lower) { lower_ = lower; }
  void set_upper(const Partype &upper) { upper_ = upper; }
  
  void minimize(Derived &func, Partype &par) {
    func.feval = 0;
    
    std::size_t npar = par.size();
    if (!control.npt) control.npt = npar + 2 <= (npar + 2) * (npar + 1) / 2 ? npar + 2 : (npar + 2) * (npar + 1) / 2; //this caused an error for 1 or two parameters, changed
    
    if(lower_.empty()){
      lower_.resize(npar);
      for(int i = 0; i< npar; i++)lower_[i] = -10000000000;
    }
    
    if(upper_.empty()){
      upper_.resize(npar);
      for(int i = 0; i< npar; i++)upper_[i] = 10000000000;
    }
    double max_par = *max_element(par.begin(),par.end());
    if (!control.rhobeg) control.rhobeg = 0.95 <= 0.2 * max_par ? 0.95 : 0.2 * max_par;
    if (!control.rhoend) control.rhoend = 1.0e-6 * control.rhobeg;
    
    if (!control.maxfun) control.maxfun = 10000;
    std::vector<double> w;
    w.resize((control.npt + 5) * (control.npt + npar) + (3 * npar * (npar + 5))/2);
    
    int res = bobyqa(npar, control.npt, minqa_objfun, &func, par.data(), lower_.data(), upper_.data(),
                     control.rhobeg, control.rhoend, control.iprint, control.maxfun, w.data());
    Update_msg(res);
    
    par_ = par;
    fval_ = func(par_);
    feval_ = func.feval;
  }
};

}

#endif
