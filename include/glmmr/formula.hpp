#pragma once

#include "general.h"
#include "interpreter.h"
#include "calculator.hpp"
#include "formulaparse.h"

namespace glmmr{

class Formula {
public:
  str                 formula_;
  std::vector<char>   linear_predictor_;
  strvec              re_;
  strvec              z_;
  intvec              re_order_;
  bool                RM_INT;
  strvec              fe_parameter_names_;
  Formula(const str& formula) : formula_(formula) {tokenise();};
  Formula(const glmmr::Formula& formula) : formula_(formula.formula_), fe_parameter_names_(formula.fe_parameter_names_) {tokenise();};
  Formula& operator= (const glmmr::Formula& formula){
    formula_ = formula.formula_;
    fe_parameter_names_ = formula.fe_parameter_names_;
    tokenise();
    return *this;
  };
  void    tokenise();
  void    formula_validate();
  void    calculate_linear_predictor(glmmr::calculator& calculator,const ArrayXXd& data,const strvec& colnames, MatrixXd& Xdata);
  strvec  re() const;
  strvec  z() const;
  strvec  re_terms() const;
private:
  strvec  re_terms_;
};
}

inline void glmmr::Formula::tokenise(){
  formula_validate();
  RM_INT = false;
  formula_.erase(std::remove_if(formula_.begin(), formula_.end(), [](unsigned char x) { return std::isspace(x); }), formula_.end());
  
  #if defined(ENABLE_DEBUG) && defined(R_BUILD)
  Rcpp::Rcout << "\nINITIALISE FORMULA: " << formula_;
  #endif
  
  // if there is a -1 then remove it and set no intercept
  std::vector<char> formula_as_chars(formula_.begin(),formula_.end());
  int bracket_count = 0;
  int cursor = 0;
  int nchar = formula_as_chars.size();
  bool has_found_symbol=false;
  while(!has_found_symbol && cursor < (nchar-1)){
    if(formula_as_chars[cursor]=='(')bracket_count++;
    if(formula_as_chars[cursor]==')')bracket_count--;
    if((formula_as_chars[cursor]=='-' && formula_as_chars[cursor+1]=='1') && bracket_count == 0){
      has_found_symbol = true;
      break;
    } 
    cursor++;
  }
  
  if(has_found_symbol){
    RM_INT = true;
    formula_.erase(cursor,2);
    formula_as_chars = std::vector<char>(formula_.begin(),formula_.end());
  } else {
    str add_intercept = "b_intercept+";
    std::vector<char> add_intercept_vec(add_intercept.begin(),add_intercept.end());
    linear_predictor_.insert(linear_predictor_.end(),add_intercept_vec.begin(),add_intercept_vec.end());
  }
  
  #if defined(ENABLE_DEBUG) && defined(R_BUILD)
    if(has_found_symbol)Rcpp::Rcout << "\nfound remove intercept: " << formula_;
  #endif
  
  if(formula_as_chars[0]=='+' && !has_found_symbol)throw std::runtime_error("Cannot start a formula with +");
  if(formula_as_chars[0]=='+' && has_found_symbol)formula_as_chars.erase(formula_as_chars.begin());
  std::vector<char> temp_token;
  int idx = 0;
  // split fixed and random effects in the formula
  // random effects terms are defined as (*|*)
  bracket_count = 0;
  cursor = 0;
  has_found_symbol = false;
  nchar = formula_as_chars.size();
  while(cursor <= nchar){
    idx = cursor == nchar ? nchar - 1 : cursor;
    if(cursor == nchar || (formula_as_chars[idx]=='+' && bracket_count == 0 && cursor < nchar)){
      if(temp_token[0]!='(' || std::find(temp_token.begin(),temp_token.end(),'|') == temp_token.end()){
        linear_predictor_.insert(linear_predictor_.end(),temp_token.begin(),temp_token.end());
        linear_predictor_.push_back('+');
      } else {
        if(temp_token.back()!=')')throw std::runtime_error("Invalid formula, no closing bracket");
        int mm = temp_token.size();
        int cursor_re = 1;
        std::vector<char> temp_token_re;
        while(cursor_re < mm){
          if(temp_token[cursor_re]=='|'){
            str re_new(temp_token_re.begin(),temp_token_re.end());
            z_.push_back(re_new);
            temp_token_re.clear();
          } else if(cursor_re == mm-1){
            str re_new(temp_token_re.begin(),temp_token_re.end());
            re_.push_back(re_new);
          } else {
            temp_token_re.push_back(temp_token[cursor_re]);
          }
          cursor_re++;
        }
      }
      temp_token.clear();
    } else {
      if(formula_as_chars[idx]=='(')bracket_count++;
      if(formula_as_chars[idx]==')')bracket_count--;
      temp_token.push_back(formula_as_chars[idx]);
    }
    cursor++;
  }
  if(linear_predictor_.back()=='+')linear_predictor_.pop_back();
  for(int i =0; i< static_cast<int>(re_.size()); i++){
    re_order_.push_back(i);
  }
  re_terms_ = re_;
  
  #if defined(ENABLE_DEBUG) & defined(R_BUILD)
    Rcpp::Rcout << "\nLinear predictor: ";
    for(const auto& i: linear_predictor_)Rcpp::Rcout << i;
    Rcpp::Rcout << "\nRandom effects: ";
    glmmr::print_vec_1d<strvec>(re_terms_);
  #endif
}

inline void glmmr::Formula::formula_validate(){
  //currently only checks if addition inside re term
  int open = 0;
  bool has_a_plus = false;
  bool has_a_vert = false;
  for(auto ch: formula_){
    if(ch=='(')open++;
    if(ch==')'){
      open--;
      if(open == 0){
        has_a_plus = false;
        has_a_vert = false;
      }
    }
    if(ch=='+' && open > 0)has_a_plus = true;
    if(ch=='|' && open > 0)has_a_vert = true;
    if(has_a_plus && has_a_vert)throw std::runtime_error("Addition inside re term not currently supported");
  }
}

inline void glmmr::Formula::calculate_linear_predictor(glmmr::calculator& calculator,const ArrayXXd& data,const strvec& colnames, MatrixXd& Xdata){
  bool outparse = glmmr::parse_formula(linear_predictor_,
                                       calculator,
                                       data,
                                       colnames,
                                       Xdata);
  (void)outparse;
  std::reverse(calculator.instructions.begin(),calculator.instructions.end());
  std::reverse(calculator.indexes.begin(),calculator.indexes.end());
}

inline strvec glmmr::Formula::re() const{
  return re_;
}

inline strvec glmmr::Formula::z() const{
  return z_;
}

inline strvec glmmr::Formula::re_terms() const{
  return re_terms_;
}
