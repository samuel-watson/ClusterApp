#pragma once

// #define R_BUILD

#include <vector>
#include <functional>
#include <exception>
#include <map>
#include <cmath>
#include <random>
#include <queue>
#include <type_traits>
#include <memory>
#include <tuple>
#include "bobyqa_algo.h"
#include "newuoa.h"
#include "lbfgs.h"
#include "lbfgsb.h"
#ifdef R_BUILD
#include <Rcpp.h> // for printing to R
#endif

// this will be used as the basis of adding more funcitons at a later date

class optim_algo {};
class DIRECT : public optim_algo {};
class BOBYQA : public optim_algo {};
class NEWUOA : public optim_algo {};
class LBFGS : public optim_algo {};



template<typename Signature, class algo>
class optim;

inline size_t random_index(size_t max_index) {
  std::random_device                      rand_dev;
  std::mt19937                            generator(rand_dev());
  std::uniform_int_distribution<size_t>   distr(0, max_index);
  return distr(generator);
}

// class for bobyqa - I am working on a new version of BOBYQA 
template <typename T>
class optim<T(const std::vector<T>&),BOBYQA> {
  using func = T(*)(long, const T*, void*); 
public:
  struct optimControl {
    int npt = 0;
    double rhobeg = 0.0;
    double rhoend = 0.0;
    int trace = 0;
    int maxfun = 0;
  } control;
  
  optim(){};
  optim(const std::vector<T>& start);
  optim(const optim& x) = default;
  auto operator=(const optim& x) -> optim& = default;
  
  // functions
  // the two functions fn() enable capturing of functions to use in the algorithm
  template<auto Function, typename = std::enable_if_t<std::is_invocable_r_v<T, decltype(Function), const std::vector<T>& > > >
  void    fn();
  template<auto Function, typename Class, typename = std::enable_if_t<std::is_invocable_r_v<T, decltype(Function), Class*, const std::vector<T>& > > >
  void    fn(Class* cls);
  
  void    set_bounds(const std::vector<T>& lower, const std::vector<T>& upper);
  void    set_bounds(const std::vector<T>& bound, bool lower = true);
  auto    operator()(const std::vector<T>& vec) const -> T;
  void    minimise(); // the optim algorithm
  std::vector<T> values() const;
  
private:
  [[noreturn]]
  static auto null_fn(long n, const T* x, void* p) -> T {throw std::exception{};}
  
  void*                   optim_instance = nullptr;  // pointer to the class if a member function
  func                    optim_fn = &null_fn;        // pointer to the function
  size_t                  dim;                       // number of dimensions
  std::vector<T>          lower_bound;               // bounds
  std::vector<T>          upper_bound;   
  T                       min_f;                     // current best value
  int                     fn_counter = 0;
  int                     iter = 0;
  std::vector<T>          current_values;
  std::string             msg_;
  
  //functions
  auto            eval(const std::vector<T>& vec) -> T;
  void            update_msg(int res);
};


template <typename T>
inline optim<T(const std::vector<T>&),BOBYQA>::optim(const std::vector<T>& start) {
  dim = start.size();
  current_values.resize(dim);
  current_values = start;
}

template <typename T>
inline std::vector<T> optim<T(const std::vector<T>&),BOBYQA>::values() const {
  return current_values;
}

template <typename T>
inline void optim<T(const std::vector<T>&),BOBYQA>::set_bounds(const std::vector<T>& lower, const std::vector<T>& upper) {
  lower_bound.resize(dim);
  upper_bound.resize(dim);
  lower_bound = lower; 
  upper_bound = upper;
};

template <typename T>
inline void optim<T(const std::vector<T>&),BOBYQA>::set_bounds(const std::vector<T>& bound, bool lower) {
  if(lower){
    lower_bound.resize(dim);
    lower_bound = bound;
  } else {
    upper_bound.resize(dim);
    upper_bound = bound;
  }
};

template <typename T>
template <auto Function, typename>
inline void optim<T(const std::vector<T>&),BOBYQA>::fn() 
{
  optim_instance = nullptr;
  optim_fn = static_cast<func>([](long n, const T* x, void*) -> T {
    return std::invoke(Function, std::vector<T>(x,x+n));
  });
};

template <typename T>
template<auto Function, typename Class, typename>
inline void optim<T(const std::vector<T>&),BOBYQA>::fn(Class* cls)
{
  optim_instance = cls;
  optim_fn = static_cast<func>([](long n, const T* x, void* p) -> T {
    auto* c = static_cast<Class*>(p);
    return std::invoke(Function,c,std::vector<T>(x,x+n));
  });
}

template <typename T>
inline auto optim<T(const std::vector<T>&),BOBYQA>::operator()(const std::vector<T>& vec) const -> T
{
  return std::invoke(optim_fn,vec.size(),vec.data(),optim_instance);
}    

template <typename T>
inline void optim<T(const std::vector<T>&),BOBYQA>::minimise(){
    fn_counter = 0;
#ifndef R_BUILD
    double R_NegInf = -1.0 * std::numeric_limits<double>::infinity();
    double R_PosInf = std::numeric_limits<double>::infinity();
#endif
    if (!control.npt) control.npt = std::min(dim + 2, (dim+2)*(dim+1)/2);     
    if(lower_bound.empty()){
      lower_bound.resize(dim);
      for(int i = 0; i< dim; i++) lower_bound[i] = R_NegInf;
    }    
    if(upper_bound.empty()){
      upper_bound.resize(dim);
      for(int i = 0; i< dim; i++)upper_bound[i] = R_PosInf;
    }
    double max_par = *std::max_element(current_values.begin(),current_values.end());
    if (!control.rhobeg) control.rhobeg = std::min(0.95, 0.2*max_par);
    if (!control.rhoend) control.rhoend = 1.0e-6 * control.rhobeg;    
    if (!control.maxfun) control.maxfun = 10000;
    std::vector<double> w;
    w.resize((control.npt + 5) * (control.npt + dim) + (3 * dim * (dim + 5))/2);    
    int res = bobyqa(dim, control.npt, optim_fn, optim_instance, current_values.data(), 
                     lower_bound.data(), upper_bound.data(),
                     control.rhobeg, control.rhoend, control.trace, 
                     control.maxfun, w.data());
    update_msg(res);
    min_f = eval(current_values);
#ifdef R_BUILD
    if(control.trace >= 1)
    {
      Rcpp::Rcout << "\nEND BOBYQA | fn: " << fn_counter << " | " << msg_;  
    }
#endif
}

template <typename T>
inline auto optim<T(const std::vector<T>&),BOBYQA>::eval(const std::vector<T>& vec) -> T
{   
  fn_counter++;
  return std::invoke(optim_fn,vec.size(),vec.data(),optim_instance);
}  

template <typename T>
inline void optim<T(const std::vector<T>&),BOBYQA>::update_msg(int res) {
    switch(res) {
    case 0:
      msg_ = "Normal exit from optim";
      break;
    case -1:
      msg_ = "optim -- NPT is not in the required interval";
      break;
    case -2:
      msg_ = "optim -- one of the box constraint ranges is too small (< 2*RHOBEG)";
      break;
    case -3:
      msg_ = "optim detected too much cancellation in denominator";
      break;
    case -4:
      msg_ = "optim -- maximum number of function evaluations exceeded";
      break;
    case -5:
      msg_ = "optim -- a trust region step failed to reduce q";
      break;
    default: ;
    }
  }


  // class for newuoa
template <typename T>
class optim<T(const std::vector<T>&),NEWUOA> {
  using func = T(*)(void*, long, const T*); 
public:
  struct optimControl {
    int npt = 0;
    double rhobeg = 0.0;
    double rhoend = 0.0;
    int trace = 0;
    int maxfun = 0;
  } control;
  
  optim(){};
  optim(const std::vector<T>& start);
  optim(const optim& x) = default;
  auto operator=(const optim& x) -> optim& = default;
  
  // functions
  // the two functions fn() enable capturing of functions to use in the algorithm
  template<auto Function, typename = std::enable_if_t<std::is_invocable_r_v<T, decltype(Function), const std::vector<T>& > > >
  void    fn();
  template<auto Function, typename Class, typename = std::enable_if_t<std::is_invocable_r_v<T, decltype(Function), Class*, const std::vector<T>& > > >
  void    fn(Class* cls);
  
  void    set_bounds(const std::vector<T>& lower, const std::vector<T>& upper);
  void    set_bounds(const std::vector<T>& bound, bool lower = true);
  auto    operator()(const std::vector<T>& vec) const -> T;
  void    minimise(); // the optim algorithm
  std::vector<T> values() const;
  
private:
  [[noreturn]]
  static auto null_fn(void* p, long n, const T* x) -> T {throw std::exception{};}
  
  void*                   optim_instance = nullptr;  // pointer to the class if a member function
  func                    optim_fn = &null_fn;        // pointer to the function
  size_t                  dim;                       // number of dimensions
  std::vector<T>          lower_bound;               // bounds
  std::vector<T>          upper_bound;   
  T                       min_f;                     // current best value
  int                     fn_counter = 0;
  int                     iter = 0;
  std::vector<T>          current_values;
  
  //functions
  auto            eval(const std::vector<T>& vec) -> T;
};


template <typename T>
inline optim<T(const std::vector<T>&),NEWUOA>::optim(const std::vector<T>& start) {
  dim = start.size();
  current_values.resize(dim);
  current_values = start;
}

template <typename T>
inline std::vector<T> optim<T(const std::vector<T>&),NEWUOA>::values() const {
  return current_values;
}

template <typename T>
inline void optim<T(const std::vector<T>&),NEWUOA>::set_bounds(const std::vector<T>& lower, const std::vector<T>& upper) {
  lower_bound.resize(dim);
  upper_bound.resize(dim);
  lower_bound = lower; 
  upper_bound = upper;
};

template <typename T>
inline void optim<T(const std::vector<T>&),NEWUOA>::set_bounds(const std::vector<T>& bound, bool lower) {
  if(lower){
    lower_bound.resize(dim);
    lower_bound = bound;
  } else {
    upper_bound.resize(dim);
    upper_bound = bound;
  }
};

template <typename T>
template <auto Function, typename>
inline void optim<T(const std::vector<T>&),NEWUOA>::fn() 
{
  optim_instance = nullptr;
  optim_fn = static_cast<func>([](void*, long n, const T* x) -> T {
    return std::invoke(Function, std::vector<T>(x,x+n));
  });
};

template <typename T>
template<auto Function, typename Class, typename>
inline void optim<T(const std::vector<T>&),NEWUOA>::fn(Class* cls)
{
  optim_instance = cls;
  optim_fn = static_cast<func>([](void* p, long n, const T* x) -> T {
    auto* c = static_cast<Class*>(p);
    return std::invoke(Function,c,std::vector<T>(x,x+n));
  });
}

template <typename T>
inline auto optim<T(const std::vector<T>&),NEWUOA>::operator()(const std::vector<T>& vec) const -> T
{
  return std::invoke(optim_fn,optim_instance,vec.size(),vec.data());
}    

template <typename T>
inline void optim<T(const std::vector<T>&),NEWUOA>::minimise(){
    fn_counter = 0;
#ifndef R_BUILD
    double R_NegInf = -1.0 * std::numeric_limits<double>::infinity();
    double R_PosInf = std::numeric_limits<double>::infinity();
#endif
    if (!control.npt) control.npt = std::min(dim + 2, (dim+2)*(dim+1)/2);     
    if(lower_bound.empty()){
      lower_bound.resize(dim);
      for(int i = 0; i< dim; i++) lower_bound[i] = R_NegInf;
    }    
    if(upper_bound.empty()){
      upper_bound.resize(dim);
      for(int i = 0; i< dim; i++)upper_bound[i] = R_PosInf;
    }
    double max_par = *std::max_element(current_values.begin(),current_values.end());
    if (!control.rhobeg) control.rhobeg = std::min(0.95, 0.2*max_par);
    if (!control.rhoend) control.rhoend = 1.0e-6 * control.rhobeg;    
    if (!control.maxfun) control.maxfun = 10000;
    std::vector<double> w;
    w.resize((control.npt + 5) * (control.npt + dim) + (3 * dim * (dim + 5))/2);   
    auto closure = NewuoaClosure{optim_instance, optim_fn};
    fn_counter = 0;
    double result = newuoa_closure( &closure, dim, control.npt, current_values.data(), control.rhobeg, control.rhoend, control.maxfun, w.data(), &fn_counter);
    min_f = eval(current_values);
#ifdef R_BUILD
    if(control.trace >= 1)
    {
      Rcpp::Rcout << "\nEND NEWUOA | fn: " << fn_counter ;  
    }
#endif
}

template <typename T>
inline auto optim<T(const std::vector<T>&),NEWUOA>::eval(const std::vector<T>& vec) -> T
{   
  fn_counter++;
  return std::invoke(optim_fn, optim_instance, vec.size(), vec.data());
}

// class for L-BFGS
template <>
class optim<double(const VectorXd&, VectorXd&),LBFGS> {
  using func = double(*)(void*, const VectorXd&, VectorXd&); 
public:
  struct optimControl {
    double g_epsilon = 1.0e-8;
    double past = 3;
    double delta = 1.0e-8;
    int max_linesearch = 64;
    int trace = 0;
  } control;
  
  optim(const VectorXd& start);
  optim(const optim& x) = default;
  auto operator=(const optim& x) -> optim& = default;
  
  // functions
  // the two functions fn() enable capturing of functions to use in the algorithm
  template<auto Function, typename = std::enable_if_t<std::is_invocable_r_v<double, decltype(Function), const VectorXd&, VectorXd& > > >
  void                fn();
  template<auto Function, typename Class, typename = std::enable_if_t<std::is_invocable_r_v<double, decltype(Function), Class*, const VectorXd&, VectorXd& > > >
  void                fn(Class* cls);

  auto                operator()(const VectorXd& vec, VectorXd& g) -> double;
  void                minimise(); // the optim algorithm
  VectorXd            values() const;
  void                set_bounds(const VectorXd& lower, const VectorXd& upper);
  void                set_bounds(const std::vector<double>& lower, const std::vector<double>& upper);
  
private:
  [[noreturn]]
  static auto null_fn(void* p, const VectorXd& x, VectorXd& g) -> double {throw std::exception{};}
  
  void*            optim_instance = nullptr;  // pointer to the class if a member function
  func             optim_fn = &null_fn;        // pointer to the function
  size_t           dim;                       // number of dimensions
  double           min_f = 0;                     // current best value
  VectorXd         current_values;
  VectorXd         lower_bound;               // bounds
  VectorXd         upper_bound;   
  int              fn_counter = 0;
  int              iter = 0;
  bool             bounded = false;

  //functions
  double           eval(const VectorXd& vec, VectorXd& g);
};


inline optim<double(const VectorXd&, VectorXd&),LBFGS>::optim(const VectorXd& start) : dim(start.size()), current_values(start), 
  lower_bound(dim), upper_bound(dim) {};

inline VectorXd optim<double(const VectorXd&, VectorXd&),LBFGS>::values() const 
{
  return current_values;
}

inline void optim<double(const VectorXd&, VectorXd&),LBFGS>::set_bounds(const VectorXd& lower, const VectorXd& upper)
{
  for(int i = 0; i < dim; i++)
  {
    lower_bound(i) = lower(i);
    upper_bound(i) = upper(i);
  }
  bounded = true;
}

inline void optim<double(const VectorXd&, VectorXd&),LBFGS>::set_bounds(const std::vector<double>& lower, const std::vector<double>& upper)
{
  
  for(int i = 0; i < dim; i++)
  {
    lower_bound(i) = lower[i];
    upper_bound(i) = upper[i];
  }
  bounded = true;
}

template <auto Function, typename>
inline void optim<double(const VectorXd&, VectorXd&),LBFGS>::fn() 
{
  optim_instance = nullptr;
  optim_fn = static_cast<func>([](void*, const VectorXd& x, VectorXd& g) -> double {
    return std::invoke(Function, x, g);
  });
};

template<auto Function, typename Class, typename>
inline void optim<double(const VectorXd&, VectorXd&),LBFGS>::fn(Class* cls)
{
  optim_instance = cls;
  optim_fn = static_cast<func>([](void* p, const VectorXd& x, VectorXd& g) -> double {
    auto* c = static_cast<Class*>(p);
    return std::invoke(Function,c,x,g);
  });
}

inline auto optim<double(const VectorXd&, VectorXd&),LBFGS>::operator()(const VectorXd& vec, VectorXd& g) -> double
{
  fn_counter++;
  return std::invoke(optim_fn,optim_instance,vec,g);
}   

inline double optim<double(const VectorXd&, VectorXd&),LBFGS>::eval(const VectorXd& vec, VectorXd& g) 
{
  return std::invoke(optim_fn,optim_instance,vec,g);
}   

inline void optim<double(const VectorXd&, VectorXd&),LBFGS>::minimise()
{
        int niter;
        fn_counter = 0;
        if(!bounded){
          LBFGSpp::LBFGSParam<double> param;
          param.epsilon = control.g_epsilon;
          param.max_linesearch = control.max_linesearch;
          param.delta = control.delta;
          param.past = control.past;
          LBFGSpp::LBFGSSolver<double> solver(param);
          niter = solver.minimize(*this, current_values, min_f, control.trace);
        } else {
          LBFGSpp::LBFGSBParam<double> param;
          param.epsilon = control.g_epsilon;
          param.max_linesearch = control.max_linesearch;
          param.delta = control.delta;
          param.past = control.past;
          LBFGSpp::LBFGSBSolver<double> solver(param);
          niter = solver.minimize(*this, current_values, min_f, lower_bound, upper_bound, control.trace);
        }

        VectorXd g(dim);
        double a = eval(current_values, g);

#ifdef R_BUILD
    if(control.trace >= 1)
    {
      Rcpp::Rcout << "\nL-BFGS END: " << niter << " iterations with " << fn_counter-1 << " function evaluations";
      Rcpp::Rcout << "\nx = " << current_values.transpose();
      Rcpp::Rcout << "\nf(x) = " << min_f;
    }
#endif
}

typedef optim<double(const std::vector<double>&),BOBYQA> bobyqad;
typedef optim<float(const std::vector<float>&),BOBYQA> bobyqaf;
typedef optim<double(const std::vector<double>&),NEWUOA> newuoad;
typedef optim<float(const std::vector<float>&),NEWUOA> newuoaf;
typedef optim<double(const VectorXd&, VectorXd&),LBFGS> lbfgsd;


// below is an implementation of the DIRECT algorithm. It works but is much inferior for the maximum likelihood problems!
// I've left it in here as it might be useful for comparison purposes, or if refinements can be identified.

enum class Position {
  Lower,
  Middle,
  Upper
};

// hyperrectangle class - coordinates should be in [0,1]^D
template <typename T>
class Rectangle {
public:
  int                 dim;
  std::vector<T>      min_x;
  std::vector<T>      max_x;
  T fn_value;
  T max_dim_size;
  bool potentially_optimal = false;

  Rectangle(){};
  Rectangle(const int dim_) : dim(dim_), min_x(dim), max_x(dim) {};
  Rectangle(const std::vector<T>& min_x_, const std::vector<T>& max_x_) : dim(min_x_.size()), min_x(min_x_), max_x(max_x_) {};
  Rectangle(const Rectangle<T>& x) : dim(x.dim), min_x(x.min_x), max_x(x.max_x) {};
  auto operator=(const Rectangle<T>& x) -> Rectangle& = default;

  // functions
  std::vector<T>        centroid(); // returns the centroid
  std::vector<T>        centroid(const size_t& dim, const T& delta); // returns the centroid offset by delta in dimension dim
  void                  unit_hyperrectangle(); //sets the rectangle to be the unit hyper-rectangle
  std::pair<T,size_t>   longest_side(); // returns 0.5 times the longest side
  T                     dim_size(const size_t& dim) const; // returns the size of the dimension
  void                  trim_dimension(const size_t& dim_t, const Position pos); // reduces the dimension to a third of the size - either lower, middle or upper
  void                  set_bounds(const std::vector<T>& min, const std::vector<T>& max); //set min x, max x
};

// class for handling function binding
template <typename T>
class optim<T(const std::vector<T>&),DIRECT> {
  using func = T(*)(const void*, const std::vector<T>&);
public:
  struct optimControl {
    T epsilon = 1e-4;
    int max_iter = 1;
    T tol = 1e-4;
    bool select_one = true; //select only one potentially optimal rectangle on each iteration
    bool trisect_once = false; // trisect only one side per division
    int trace = 0;
    int max_eval = 0;
    bool mrdirect = false; // use a multilevel refinement process
    T l2_tol = 1e-2;
    T l1_tol = 1e-4;
    T l2_epsilon = 1e-5;
    T l1_epsilon = 1e-7;
    T l0_epsilon = 0;
  } control;

  optim(){};
  // if starting vals is true, it assumes x is a vector of starting values, and y sets the bounds by adding +/- either side of x
  // otherwise it assumes x is a lower bound, and y is an upper bound
  optim(const std::vector<T>& x, const std::vector<T>& y, bool starting_vals = true);
  optim(const std::vector<T>& x);
  optim(const optim& x) = default;
  auto operator=(const optim& x) -> optim& = default;

  // functions
  // the two functions fn() enable capturing of functions to use in the algorithm
  template<auto Function, typename = std::enable_if_t<std::is_invocable_r_v<T, decltype(Function), const std::vector<T>& > > >
  void    fn();
  template<auto Function, typename Class, typename = std::enable_if_t<std::is_invocable_r_v<T, decltype(Function), Class*, const std::vector<T>& > > >
  void    fn(Class* cls);

  void            set_bounds(const std::vector<T>& lower, const std::vector<T>& upper, bool starting_vals = true);
  auto            operator()(const std::vector<T>& vec) const -> T;
  void            minimise(); // the optim algorithm
  std::vector<T>  values() const;
  T               rect_size() const;

private:
  [[noreturn]]
  static auto null_fn(const void* p, const std::vector<T>& vec) -> T {throw std::exception{};}

  const void*               optim_instance = nullptr;  // pointer to the class if a member function
  func                      optim_fn = &null_fn;        // pointer to the function
  size_t                    dim;                       // number of dimensions
  std::vector<T>            lower_bound;               // bounds
  std::vector<T>            upper_bound;
  std::vector<T>            dim_size;                  // size of each dimension to transform to unit rectangle
  std::vector<std::unique_ptr<Rectangle<T>>> rects;   // the rectangles
  T                         min_f = 0;               // current best value
  int                       fn_counter = 0;
  int                       iter = 0;
  std::vector<T>            current_values;
  std::pair<T,size_t>       current_largest_dim;
  size_t                    mrdirect_level = 2;
  T                         max_diff = control.tol * 1.1;

  //functions
  auto            eval(const std::vector<T>& vec) -> T;
  std::vector<T>  transform(const std::vector<T>& vec);
  size_t          update_map();
  void            filter_rectangles(size_t n_optimal);
  void            divide_rectangles();
};

template <typename T>
inline std::vector<T> Rectangle<T>::centroid()
{
  std::vector<T> centre(dim);
  for(size_t i = 0; i < dim; i++){
    centre[i] = 0.5*(max_x[i] + min_x[i]);
  }
  return(centre);
};

template <typename T>
inline void Rectangle<T>::unit_hyperrectangle()
{
  std::fill(max_x.begin(),max_x.end(),1.0);
  std::fill(min_x.begin(),min_x.end(),0.0);
};

template <typename T>
inline std::vector<T> Rectangle<T>::centroid(const size_t& dim_ex, const T& delta)
{
  std::vector<T> centre(dim);
  for(size_t i = 0; i < dim; i++){
    centre[i] = 0.5*(max_x[i] + min_x[i]);
    if(i == dim_ex) centre[i] += delta;
  }
  return(centre);
};

template <typename T>
inline void Rectangle<T>::set_bounds(const std::vector<T>& min, const std::vector<T>& max)
{
  dim = min.size();
  min_x = min;
  max_x = max;
}

template <typename T>
inline std::pair<T,size_t> Rectangle<T>::longest_side()
{
  T long_len = 0;
  size_t which_dim;
  for(size_t i = 0; i < dim; i++){
    T diff = max_x[i] - min_x[i];
    if(diff > long_len) {
      long_len = diff;
      which_dim = i;
    }
  }
  return std::pair<T,size_t>{0.5*long_len,which_dim};
};

template <typename T>
inline T Rectangle<T>::dim_size(const size_t& dim_) const
{
  return max_x[dim_] - min_x[dim_];
};

template <typename T>
inline void Rectangle<T>::trim_dimension(const size_t& dim_t, const Position pos)
{
  T dsize = dim_size(dim_t);
  switch(pos)
  {
  case Position::Lower:
    max_x[dim_t] -= 2.0*dsize/3.0;
    break;
  case Position::Middle:
    min_x[dim_t] += dsize/3.0;
    max_x[dim_t] -= dsize/3.0;
    break;
  case Position::Upper:
    min_x[dim_t] += 2.0*dsize/3.0;
    break;
  }
};

template <typename T>
inline optim<T(const std::vector<T>&),DIRECT>::optim(const std::vector<T>& x, const std::vector<T>& y, bool starting_vals)
{
  set_bounds(x,y,starting_vals);
};

template <typename T>
inline optim<T(const std::vector<T>&),DIRECT>::optim(const std::vector<T>& x) : dim(x.size()), current_values(x) {};

template <typename T>
inline std::vector<T> optim<T(const std::vector<T>&),DIRECT>::values() const
{
  return current_values;
}

template <typename T>
inline T optim<T(const std::vector<T>&),DIRECT>::rect_size() const
{
  return max_diff;
}

template <typename T>
inline void optim<T(const std::vector<T>&),DIRECT>::set_bounds(const std::vector<T>& x, const std::vector<T>& y, bool starting_vals)
{
  dim = x.size();
  lower_bound.resize(dim);
  upper_bound.resize(dim);
  dim_size.resize(dim);
  if(!starting_vals)
  {
    lower_bound = x;
    upper_bound = y;
    for(size_t i = 0; i < dim; i++) dim_size[i] = y[i] - x[i];
  } else {
    for(size_t i = 0; i < dim; i++){
      lower_bound[i] = x[i] - y[i];
      upper_bound[i] = x[i] + y[i];
      dim_size[i] = 2*y[i];
    }
  }
  current_values.resize(dim);
  std::fill(current_values.begin(),current_values.end(),0.0);
  rects.push_back(std::make_unique<Rectangle<T>>(dim));
  rects.back()->unit_hyperrectangle();
  rects.back()->max_dim_size = 0.5;
  current_largest_dim = rects.back()->longest_side();
};

template <typename T>
template <auto Function, typename>
inline void optim<T(const std::vector<T>&),DIRECT>::fn()
{
  optim_instance = nullptr;
  optim_fn = static_cast<func>([](const void*, const std::vector<T>& vec) -> T {
    return std::invoke(Function, vec);
  });
};

template <typename T>
template<auto Function, typename Class, typename>
inline void optim<T(const std::vector<T>&),DIRECT>::fn(Class* cls)
{
  optim_instance = cls;
  optim_fn = static_cast<func>([](const void* p, const std::vector<T>& vec) -> T {
    auto* c = const_cast<Class*>(static_cast<const Class*>(p));
    return std::invoke(Function,c,vec);
  });
}

template <typename T>
inline auto optim<T(const std::vector<T>&),DIRECT>::operator()(const std::vector<T>& vec) const -> T
{
  return std::invoke(optim_fn,optim_instance,vec);
}

template <typename T>
inline void optim<T(const std::vector<T>&),DIRECT>::minimise()
{
#ifdef R_BUILD
  if(control.trace >= 1){
    Rcpp::Rcout << "\nSTARTING optim-L";
    Rcpp::Rcout << "\nTolerance: " << control.tol << " | Max iter : " << control.max_iter << "\n Starting values :";
    std::vector<T> vals = transform(rects.front()->centroid());
    for(const auto& val: vals) Rcpp::Rcout << val << " ";
  }
#endif
  current_values = transform(rects.back()->centroid());
  rects.back()->fn_value = eval(current_values);
  min_f = rects.back()->fn_value;
  max_diff = control.tol*1.1;
  iter = 0;
  fn_counter = 0;
  bool min_check = true;
  if(control.mrdirect) control.epsilon = control.l2_epsilon;

  while(max_diff > control.tol && iter <= control.max_iter && min_check){

#ifdef R_BUILD
    if(control.trace >= 1)
    {
      Rcpp::Rcout << "\n---------------------------------------------------------------------------------- ";
      Rcpp::Rcout << "\n| Iter: " << iter << " | Evaluations: " << fn_counter << " | Rectangles: " << rects.size() << " | Dimensions: " << dim << " | Start fn: " << min_f << " |";
    }
#endif

    size_t n_optimal = update_map();
    if(control.select_one) filter_rectangles(n_optimal);
    divide_rectangles();
    max_diff = 2 * current_largest_dim.first * dim_size[current_largest_dim.second];

#ifdef R_BUILD
    if(control.trace >= 1)
    {
      Rcpp::Rcout << "\n| New best fn: " << min_f << " | Max difference: " << max_diff << " | New values: ";
      for(const auto& val: current_values) Rcpp::Rcout << val << " ";
      Rcpp::Rcout << " |\n----------------------------------------------------------------------------------";
    }
#endif
    // erase the rectangles from the size map
    iter++;
    if(control.max_eval > 0 && fn_counter > control.max_eval) min_check = false;
    if(control.mrdirect)
    {
      switch(mrdirect_level)
      {
      case 2:
    {
      if(max_diff < control.l2_tol)
    {
#ifdef R_BUILD
      if(control.trace >= 2)
      {
        Rcpp::Rcout << "\nMRDIRECT Shrinking epsilon (level 1)";
      }
#endif
      control.epsilon = control.l1_tol;
      mrdirect_level--;
      std::sort(rects.begin(), rects.end(), [](const std::unique_ptr<Rectangle<T>>& x, const std::unique_ptr<Rectangle<T>>& y){
        if(x->max_dim_size != y->max_dim_size)
          return (x->max_dim_size < y->max_dim_size);
        return x->fn_value < y->fn_value;
      });
      size_t new_size = 0.9*rects.size();
      rects.resize(new_size);
    }
      break;
    }
      case 1:
    {
      if(max_diff < control.l1_tol)
    {
#ifdef R_BUILD
      if(control.trace >= 2)
      {
        Rcpp::Rcout << "\nMRDIRECT Shrinking epsilon (level 0)";
      }
#endif
      control.epsilon = 0;
      mrdirect_level--;
      std::sort(rects.begin(), rects.end(), [](const std::unique_ptr<Rectangle<T>>& x, const std::unique_ptr<Rectangle<T>>& y){
        if(x->max_dim_size != y->max_dim_size)
          return (x->max_dim_size < y->max_dim_size);
        return x->fn_value < y->fn_value;
      });
      size_t new_size = 0.1*rects.size();
      rects.resize(new_size);
    }
      break;
    }
      case 0:
        break;
      };
    }
  }
}

template <typename T>
inline auto optim<T(const std::vector<T>&),DIRECT>::eval(const std::vector<T>& vec) -> T
{
  fn_counter++;
  return std::invoke(optim_fn,optim_instance,vec);
}

template <typename T>
inline std::vector<T> optim<T(const std::vector<T>&),DIRECT>::transform(const std::vector<T>& vec)
{
  // transform from the unit hyperrectangle
  std::vector<T> transformed_vec(dim);
  for(int i = 0; i < dim; i++){
    transformed_vec[i] = vec[i]*dim_size[i] + lower_bound[i];
  }
  return transformed_vec;
};

// after running this function size_fn should contain the potentially optimal rectangles
template <typename T>
inline size_t optim<T(const std::vector<T>&),DIRECT>::update_map()
{

  std::sort(rects.begin(), rects.end(), [](const std::unique_ptr<Rectangle<T>>& x, const std::unique_ptr<Rectangle<T>>& y){
    if(x->max_dim_size != y->max_dim_size)
      return (x->max_dim_size < y->max_dim_size);
    return x->fn_value > y->fn_value;
  });


  // variables
  std::pair<T,T>      coord = {0.0, min_f - control.epsilon*abs(min_f)};
  size_t              index = 0;
  size_t              end = rects.size();
  T                   x, y, angle, new_angle;
  size_t              iter, min_index;
  size_t              n_potentially_optimal = 0;

  // function body
  while(index < end)
  {
    if(index == end - 1)
    {
      rects[index]->potentially_optimal = true;
      n_potentially_optimal++;
      index++;
    } else {
      iter = index;
      min_index = index;
      angle = M_PI*0.5;
      while(iter < end)
      {
        y = abs(rects[iter]->fn_value - coord.second);
        x = abs(rects[iter]->max_dim_size - coord.first);
        new_angle = atan(y/x);
        if(new_angle < angle){
          min_index = iter;
          angle = new_angle;
        }
        iter++;
      }
#ifdef R_BUILD
      if(control.trace >= 2)
      {
        Rcpp::Rcout << "\nNEXT POTENTIALLY OPTIMAL: (" << coord.first << ", " << coord.second << ") => (" << min_index << ": " << rects[min_index]->max_dim_size << ", " << rects[min_index]->fn_value << ")";
      }
#endif
      index = min_index;
      rects[index]->potentially_optimal = true;
      coord.second = rects[index]->fn_value;
      coord.first = rects[index]->max_dim_size;
      index++;
      n_potentially_optimal++;
    }
  }

  return n_potentially_optimal;
};


// selects just one potentially optimal rectangle
template <typename T>
inline void optim<T(const std::vector<T>&),DIRECT>::filter_rectangles(size_t n_optimal)
{
  if(n_optimal > 1){
    size_t keep_index = random_index(n_optimal - 1);
    size_t counter = 0;
    for(auto& r: rects){
      if(r->potentially_optimal){
        if(counter != keep_index)r->potentially_optimal = false;
        counter++;
      }
    }
  }
}

// divides up the potentially optimal rectangles
// it adds the rectangles to rects, and then erases the original
// rectangles, while also clearing size_fn
template <typename T>
inline void optim<T(const std::vector<T>&),DIRECT>::divide_rectangles()
{
  //identify the largest dimensions
  typedef std::pair<std::pair<T,T>,size_t> dimpair;

  struct compare_pair {
    bool operator()(const dimpair& elt1, const dimpair& elt2) const {
      return std::min(elt1.first.first,elt1.first.second) < std::min(elt2.first.first,elt2.first.second);
    };
  };

  std::vector<size_t> largest_dims;
  std::priority_queue< dimpair, std::vector<dimpair>, compare_pair > pq;
  std::vector<std::unique_ptr<Rectangle<T>>> new_rectangles;
  size_t counter = 0;

  for(auto& r: rects){
    if(r->potentially_optimal)
    {
      largest_dims.clear();
      for(size_t i = 0; i < dim; i++)
      {
        if(r->dim_size(i) == 2*r->max_dim_size) largest_dims.push_back(i);
      }

      T delta = 2 * r->max_dim_size / 3.0;

#ifdef R_BUILD
      if(largest_dims.size()==0)Rcpp::stop("No dimension data");
      if(control.trace >= 2)
      {
        Rcpp::Rcout << "\nDIVIDING RECTANGLE " << counter << " | Largest dim size: " << r->max_dim_size << " delta: " << delta << " in dimensions: ";
        for(const auto& val: largest_dims) Rcpp::Rcout << val << " ";
        Rcpp::Rcout << " | fn : " << r->fn_value;
      }
#endif

      T fn1, fn2;
      for(const auto& dd: largest_dims)
      {
        fn1 = eval(transform(r->centroid(dd, delta)));
        fn2 = eval(transform(r->centroid(dd, -delta)));
        pq.push(dimpair(std::pair<T,T>(fn1,fn2), dd));
      }

      bool dim_control = true;

      while(!pq.empty() && dim_control){
        size_t dim_vvv = pq.top().second;
        new_rectangles.push_back(std::make_unique<Rectangle<T>>(dim));
        new_rectangles.back()->set_bounds(r->min_x,r->max_x);
        new_rectangles.back()->trim_dimension(dim_vvv,Position::Lower);
        new_rectangles.back()->fn_value = pq.top().first.second;
        new_rectangles.back()->max_dim_size = new_rectangles.back()->longest_side().first;

        if(new_rectangles.back()->fn_value <= min_f)
        {
          current_values = transform(new_rectangles.back()->centroid());
          current_largest_dim = new_rectangles.back()->longest_side();
          min_f = new_rectangles.back()->fn_value;
        }

        new_rectangles.push_back(std::make_unique<Rectangle<T>>(dim));
        new_rectangles.back()->set_bounds(r->min_x,r->max_x);
        new_rectangles.back()->trim_dimension(dim_vvv,Position::Upper);
        new_rectangles.back()->fn_value = pq.top().first.first;
        new_rectangles.back()->max_dim_size = new_rectangles.back()->longest_side().first;

        if(new_rectangles.back()->fn_value <= min_f)
        {
          current_values = transform(new_rectangles.back()->centroid());
          current_largest_dim = new_rectangles.back()->longest_side();
          min_f = new_rectangles.back()->fn_value;
        }

        r->trim_dimension(dim_vvv, Position::Middle);
        pq.pop();
        if(control.trisect_once)dim_control = false;
      }
      r->potentially_optimal = false;
      r->max_dim_size = r->longest_side().first;

      if(r->fn_value <= min_f){
        current_values = transform(r->centroid());
        current_largest_dim = r->longest_side();
        min_f = r->fn_value;
      }
    }
    counter++;
  }

  // insert new rectangles
  for(int i = 0; i < new_rectangles.size(); i++) rects.push_back(std::move(new_rectangles[i]));
  new_rectangles.clear();
};


typedef optim<double(const std::vector<double>&),DIRECT> directd;
typedef optim<float(const std::vector<float>&),DIRECT> directf;
