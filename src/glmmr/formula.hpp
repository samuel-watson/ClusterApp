#pragma once

#include "general.h"
#include "interpreter.h"
#include "calculator.hpp"
#include "formulaparse.h"

namespace glmmr{

class Formula {
public:
  str formula_;
  std::vector<char> linear_predictor_;
  strvec re_;
  strvec z_;
  intvec re_order_;
  bool RM_INT;
  Formula(const str& formula) : formula_(formula) {tokenise();};
  Formula(const glmmr::Formula& formula) : formula_(formula.formula_){tokenise();};
  Formula& operator= (const glmmr::Formula& formula){
    formula_ = formula.formula_;
    tokenise();
    return *this;
  };
  void tokenise();
  void formula_validate();
  void calculate_linear_predictor(glmmr::calculator& calculator,const ArrayXXd& data,const strvec& colnames, MatrixXd& Xdata);
  strvec re();
  strvec z();
  strvec re_terms();
private:
  strvec re_terms_;
};
}

inline void glmmr::Formula::tokenise(){
  formula_validate();
  RM_INT = false;
  formula_.erase(std::remove_if(formula_.begin(), formula_.end(), [](unsigned char x) { return std::isspace(x); }), formula_.end());
  auto minone = formula_.find("-1");
  if(minone != str::npos){
    RM_INT = true;
    formula_.erase(minone,2);
  } else {
    str add_intercept = "b_intercept*1+";
    std::vector<char> add_intercept_vec(add_intercept.begin(),add_intercept.end());
    linear_predictor_.insert(linear_predictor_.end(),add_intercept_vec.begin(),add_intercept_vec.end());
  }
  std::vector<char> formula_as_chars(formula_.begin(),formula_.end());
  int nchar = formula_as_chars.size();
  int cursor = 0;
  int bracket_count = 0; // to deal with opening of brackets
  std::vector<char> temp_token;
  int idx = 0;
  while(cursor <= nchar){
    idx = cursor == nchar ? nchar - 1 : cursor;
    if(cursor == nchar || (formula_as_chars[idx]=='+' && bracket_count == 0 && cursor < nchar)){
      if(temp_token[0]!='('){
        linear_predictor_.insert(linear_predictor_.end(),temp_token.begin(),temp_token.end());
        linear_predictor_.push_back('+');
      } else {
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
  for(unsigned int i =0; i<re_.size(); i++){
    re_order_.push_back(i);
  }
  re_terms_ = re_;
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

inline strvec glmmr::Formula::re(){
  return re_;
}

inline strvec glmmr::Formula::z(){
  return z_;
}

inline strvec glmmr::Formula::re_terms(){
  return re_terms_;
}
