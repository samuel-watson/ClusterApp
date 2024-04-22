#pragma once

#include "general.h"
#include "calculator.hpp"

namespace glmmr{

inline bool check_data(str& formula,
                       glmmr::calculator& calc,
                       const ArrayXXd& data,
                       const strvec& colnames,
                       MatrixXd& Xdata,
                       bool push = true){
  bool variable_in_data = false;
  auto colidx = std::find(colnames.begin(),colnames.end(),formula);
  if(colidx != colnames.end()){
    variable_in_data = true;
    // token is the name of a variable
    if(push) calc.instructions.push_back(Do::PushData);
    auto dataidx = std::find(calc.data_names.begin(),calc.data_names.end(),formula);
    // check if the data has already been added to Xdata
    if(dataidx != calc.data_names.end()){
      int data_index = dataidx - calc.data_names.begin();
      calc.indexes.push_back(data_index);
    } else {
      calc.data_names.push_back(formula);
      int column_index = colidx - colnames.begin();
      calc.indexes.push_back(calc.data_count);
      if(Xdata.cols() <= calc.data_count) Xdata.conservativeResize(NoChange,calc.data_count+1);
      Xdata.col(calc.data_count) = data.col(column_index);
      calc.data_count++;
    }
  }
  return variable_in_data;
}

inline bool check_parameter(str& token_as_str,
                            glmmr::calculator& calc,
                            bool bracket_flag = false){
  
  bool added_a_parameter = true;
#if defined(ENABLE_DEBUG) && defined(R_BUILD)
  Rcpp::Rcout << " parameter name: " << token_as_str;
#endif
  // interpret any other string as the name of a parameter
  // check if the parameter already exists
  calc.instructions.push_back(Do::PushParameter);
  auto find_parameter = std::find(calc.parameter_names.begin(),calc.parameter_names.end(),token_as_str);
  if(find_parameter != calc.parameter_names.end()){
    int param_position = find_parameter - calc.parameter_names.begin();
    
#if defined(ENABLE_DEBUG) && defined(R_BUILD)
    Rcpp::Rcout << "\nRepeated parameter " << token_as_str << " position: " << param_position << "\nCurrent parameters: ";
    for(const auto& p: calc.parameter_names) Rcpp::Rcout << p << " ";
#endif
    
    calc.indexes.push_back(param_position);
  } else {
    calc.parameter_names.push_back(token_as_str);
    calc.indexes.push_back(calc.parameter_count);
    calc.parameter_count++;
  }
  if(bracket_flag)calc.any_nonlinear = true;
  return added_a_parameter;
}

inline bool check_number(str& token_as_str,
                         glmmr::calculator& calc){
  bool added_a_number = glmmr::is_number(token_as_str);
  
  if(added_a_number) {
    double a = std::stod(token_as_str); 
    
    #if defined(ENABLE_DEBUG) && defined(R_BUILD)
    Rcpp::Rcout << " double = " << a << " push index = " << instruction_str.at(static_cast<Do>(calc.user_number_count));
    #endif
    
    #ifdef R_BUILD
    if(calc.user_number_count >=20)Rcpp::stop("Only ten user numbers currently permitted.");
    #endif
    
    calc.instructions.push_back(static_cast<Do>(calc.user_number_count));
    calc.numbers[calc.user_number_count] = a;
    calc.user_number_count++;
  }
  return added_a_number;
}

inline void add_factor(std::vector<char>& s2,
                       glmmr::calculator& calc,
                       const ArrayXXd& data,
                       const strvec& colnames,
                       MatrixXd& Xdata){
  
  str token_as_str2(s2.begin(),s2.end());
  auto colidx = std::find(colnames.begin(),colnames.end(),token_as_str2);
  if(colidx != colnames.end()){
    int column_index = colidx - colnames.begin();
    dblvec unique_values(data.col(column_index).data(),data.col(column_index).data()+data.rows());
    std::sort(unique_values.begin(),unique_values.end());
    auto last = std::unique(unique_values.begin(),unique_values.end());
    unique_values.erase(last,unique_values.end());
    //check for intercept
    auto findintercept = std::find(calc.parameter_names.begin(),calc.parameter_names.end(),"b_intercept");
    int factorrange = findintercept == calc.parameter_names.end() ? unique_values.size() : unique_values.size()-1;
    for(int i = 0; i < factorrange; i++){
      if(i < (factorrange - 1))calc.instructions.push_back(Do::Add);
      calc.instructions.push_back(Do::Multiply);
      if(Xdata.cols()<=calc.data_count)Xdata.conservativeResize(NoChange,calc.data_count+1);
      for(int j = 0; j < data.rows(); j++){
        Xdata(j,calc.data_count) = data(j,column_index)==unique_values[i] ? 1.0 : 0.0;
      }
      calc.indexes.push_back(calc.data_count);
      calc.data_count++;
      calc.instructions.push_back(Do::PushData);
      // parameter
      calc.instructions.push_back(Do::PushParameter);
      str dataname = token_as_str2 + "_" + std::to_string(unique_values[i])[0];
      str parname = "b_" + dataname;
      calc.parameter_names.push_back(parname);
      calc.data_names.push_back(dataname);
      calc.indexes.push_back(calc.parameter_count);
      calc.parameter_count++;
    }
    // added_a_parameter = true;
  } else {
    #ifdef R_BUILD
    Rcpp::stop("Factor variable " + token_as_str2 + " not in data");
    #endif
  }
}

inline void sign_fn(std::vector<char>& formula,
                    glmmr::calculator& calc,
                    const ArrayXXd& data,
                    const strvec& colnames,
                    MatrixXd& Xdata,
                    int type){
  
  str token_as_str(formula.begin(),formula.end());
  if(type == 0){
    calc.instructions.push_back(Do::SignNoZero);
  } else {
    calc.instructions.push_back(Do::Sign);
  }
  bool variable_in_data = check_data(token_as_str,calc,data,colnames,Xdata,false);
  if(!variable_in_data){
#ifdef R_BUILD
    Rcpp::stop("Syntax error in sign: "+token_as_str+" not in data");
#endif
  }
}

inline void two_way_fn(std::vector<char>& formula,
                       glmmr::calculator& calc,
                       const ArrayXXd& data,
                       const strvec& colnames,
                       MatrixXd& Xdata,
                       int type){

    // split at ,
  std::vector<char> f_s2;
  std::vector<char> nu_val;
  std::vector<char> l_val;
  std::vector<char> k_val;
  int comma_count = 0;
  for(int i = 0; i < formula.size(); i++)
  {
    if(formula[i]==','){
      comma_count++;
    } else {
      if(comma_count == 0){
        f_s2.push_back(formula[i]);
      } else if(comma_count == 1){
        nu_val.push_back(formula[i]);
      } else if(comma_count == 2){
        k_val.push_back(formula[i]);
      } else if(comma_count == 3){
        l_val.push_back(formula[i]);
      } else {
        #ifdef R_BUILD
        Rcpp::stop("Syntax error in twoway: too many commas");
        #endif
      }
    }
  }
#ifdef R_BUILD
  if(comma_count != 3)Rcpp::stop("Syntax error in twoway: incorrect number of options specified");
#endif
  double l = 0;
  str token_as_str(f_s2.begin(),f_s2.end());
  str nu_as_str(nu_val.begin(),nu_val.end());
  str l_as_str(l_val.begin(),l_val.end());
  str k_as_str(k_val.begin(),k_val.end());
  
  if(glmmr::is_number(l_as_str)) {
    l = std::stod(l_as_str); 
  } else {
    #ifdef R_BUILD
    Rcpp::stop("Syntax error in twoway: l is not a number");
    #endif
  }
  
  // str par1 = "b_twoway_k";
  str par2 = "b_twoway_del_i";
  str par3 = "b_twoway_del_e";
  str par4 = "b_twoway_eff";
  bool add_check, variable_in_data;
  
  // needs to all be in reverse ;-(
  calc.instructions.push_back(Do::Multiply);
  add_check = check_parameter(par4,calc,true);
  calc.instructions.push_back(Do::Power);
  calc.instructions.push_back(Do::Subtract);
  calc.instructions.push_back(Do::Int1);
  calc.instructions.push_back(Do::Power);
  calc.instructions.push_back(Do::Multiply);
  calc.instructions.push_back(static_cast<Do>(calc.user_number_count));
  calc.numbers[calc.user_number_count] = -1.0 / l;
  calc.user_number_count++;
  if(type > 0){
    calc.instructions.push_back(Do::Multiply);
    sign_fn(f_s2,calc,data,colnames,Xdata,0);
  }
  calc.instructions.push_back(Do::Log);
  calc.instructions.push_back(Do::Add);
  calc.instructions.push_back(Do::Exp);
  calc.instructions.push_back(Do::Multiply);
  calc.instructions.push_back(static_cast<Do>(calc.user_number_count));
  if(type == 0){
    calc.numbers[calc.user_number_count] = -1.0 * l;
    calc.user_number_count++;
    calc.instructions.push_back(Do::Divide);
    variable_in_data = check_data(token_as_str,calc,data,colnames,Xdata);
    add_check = check_parameter(par3,calc,true);
  } else if(type == 1) {
    calc.numbers[calc.user_number_count] = -0.5 * l;
    calc.user_number_count++;
    calc.instructions.push_back(Do::Multiply);
    sign_fn(f_s2,calc,data,colnames,Xdata,0);
    calc.instructions.push_back(Do::Add);
    calc.instructions.push_back(Do::Divide);
    variable_in_data = check_data(token_as_str,calc,data,colnames,Xdata);
    add_check = check_parameter(par3,calc,true);
    calc.instructions.push_back(Do::Int1);
  } else if(type == 2) {
    calc.numbers[calc.user_number_count] = -1.0 * l;
    calc.user_number_count++;
    calc.instructions.push_back(Do::Multiply);
    sign_fn(f_s2,calc,data,colnames,Xdata,0);
    calc.instructions.push_back(Do::Divide);
    calc.instructions.push_back(Do::Add);
    variable_in_data = check_data(token_as_str,calc,data,colnames,Xdata);
    add_check = check_parameter(par2,calc,true);
    calc.instructions.push_back(Do::Add);
    add_check = check_parameter(par3,calc,true);
    add_check = check_parameter(par2,calc,true);
  }
  calc.instructions.push_back(Do::Exp);
  if(type > 0){
    calc.instructions.push_back(Do::Multiply);
    calc.instructions.push_back(static_cast<Do>(calc.user_number_count));
    calc.numbers[calc.user_number_count] = -0.5 * l;
    calc.user_number_count++;
    calc.instructions.push_back(Do::Add);
    sign_fn(f_s2,calc,data,colnames,Xdata,0);
    calc.instructions.push_back(Do::Int1);
  } else {
    calc.instructions.push_back(static_cast<Do>(calc.user_number_count));
    calc.numbers[calc.user_number_count] = -1.0 * l;
    calc.user_number_count++;
  }
  // add_check = check_parameter(par1,calc,true);
  add_check = check_number(k_as_str, calc);
  if(!add_check){
#ifdef R_BUILD
    Rcpp::stop("Syntax error in twoway: k is not a number");
#endif
  }
  add_check = check_number(nu_as_str, calc);
  if(!add_check){
    #ifdef R_BUILD
    Rcpp::stop("Syntax error in twoway: nu is not a number");
    #endif
  }
}




// this is the main function to parse the formulae
// it is recursive
inline bool parse_formula(std::vector<char>& formula,
                          glmmr::calculator& calc,
                          const ArrayXXd& data,
                          const strvec& colnames,
                          MatrixXd& Xdata,
                          bool bracket_flag = false){
  
  #if defined(ENABLE_DEBUG) && defined(R_BUILD)
  Rcpp::Rcout << "\nFORMULA PARSE\nFormula: ";
  for(const auto& i: formula) Rcpp::Rcout << i;
  #endif
  
  bool added_a_parameter = false;
  bool s1_check, s2_check;
  int bracket_count = 0;
  int cursor = 0;
  int nchar = formula.size();
  bool has_found_symbol=false;
  std::vector<char> s1;
  std::vector<char> s2;
  // step 1: split at first +
  while(!has_found_symbol && cursor < nchar){
    #ifdef R_BUILD
    if(cursor==0 && (formula[cursor]=='+'))Rcpp::stop("Error in formula, + symbol in wrong place");
    #endif
    if(formula[cursor]=='(')bracket_count++;
    if(formula[cursor]==')')bracket_count--;
    if((formula[cursor]=='+' || (formula[cursor]=='-' && cursor > 0)) && bracket_count == 0){
      has_found_symbol = true;
      break;
    } else {
      s1.push_back(formula[cursor]);
    }
    cursor++;
  }
  if(has_found_symbol){
    #if defined(ENABLE_DEBUG) && defined(R_BUILD)
    Rcpp::Rcout << " Split at +/-: " << formula[cursor];
    #endif
    // split at +/-
    if(formula[cursor]=='+'){
      calc.instructions.push_back(Do::Add);
    } else if(formula[cursor]=='-'){
      calc.instructions.push_back(Do::Subtract);
    } else {
      #ifdef R_BUILD
      Rcpp::stop("Oops, something has gone wrong (f1)");
      #endif
    }
    cursor++;
    while(cursor < nchar){
      s2.push_back(formula[cursor]);
      cursor++;
    }
    // check first whether s1 or s2 is the name of a data column
    str s1_as_str(s1.begin(),s1.end());
    auto col_idx = std::find(colnames.begin(),colnames.end(),s1_as_str);
    if(col_idx != colnames.end()){
      str s1_parname = "b_" + s1_as_str + "*";
      s1.insert(s1.begin(),s1_parname.begin(),s1_parname.end());
    }
    str s2_as_str(s2.begin(),s2.end());
    col_idx = std::find(colnames.begin(),colnames.end(),s2_as_str);
    if(col_idx != colnames.end()){
      str s2_parname = "b_" + s2_as_str + "*";
      s2.insert(s2.begin(),s2_parname.begin(),s2_parname.end());
    }
    s1_check = parse_formula(s1,calc,data,colnames,Xdata,bracket_flag);
    s2_check = parse_formula(s2,calc,data,colnames,Xdata,bracket_flag);
    added_a_parameter = s1_check || s2_check;
  } else {
    // no +/- to split at, try *//
    s1.clear();
    s2.clear();
    cursor=0;
    bracket_count = 0;
    has_found_symbol = false;
    while(!has_found_symbol && cursor < nchar){
      #ifdef R_BUILD
      if(cursor==0 && (formula[cursor]=='*' || formula[cursor]=='/'))Rcpp::stop("Error in formula, multiply/divide symbol in wrong place");
      #endif
      if(formula[cursor]=='(')bracket_count++;
      if(formula[cursor]==')')bracket_count--;
      if((formula[cursor]=='*' || formula[cursor]=='/') && bracket_count == 0){
        has_found_symbol = true;
        break;
      } else {
        s1.push_back(formula[cursor]);
      }
      cursor++;
    }
    if(has_found_symbol){
      #if defined(ENABLE_DEBUG) && defined(R_BUILD)
      Rcpp::Rcout << " Split at */ : " << formula[cursor];
      #endif
      // split at *//
      if(formula[cursor]=='*'){
        calc.instructions.push_back(Do::Multiply);
      } else if(formula[cursor]=='/'){
        calc.instructions.push_back(Do::Divide);
      } else {
        #ifdef R_BUILD
        Rcpp::stop("Oops, something has gone wrong (f2)");
        #endif
      }
      cursor++;
      while(cursor < nchar){
        s2.push_back(formula[cursor]);
        cursor++;
      }
      s2.insert(s2.begin(),'('); // hacky way of avoiding a crash if s2 is a parameter - no idea why it happens
      s2.insert(s2.end(),')');
      s1_check = parse_formula(s1,calc,data,colnames,Xdata,bracket_flag);
      s2_check = parse_formula(s2,calc,data,colnames,Xdata,bracket_flag);
      if((s1_check && s2_check) || (s2_check && calc.instructions.back() == Do::Divide)) calc.any_nonlinear = true;
      added_a_parameter = s1_check || s2_check;
    } else {
      // no * to split at, try pow
      s1.clear();
      s2.clear();
      cursor=0;
      bracket_count = 0;
      while(!has_found_symbol && cursor < nchar){
        #ifdef R_BUILD
        if(cursor==0 && formula[cursor]=='^')Rcpp::stop("Error in formula, ^ symbol in wrong place");
        #endif
        
        if(formula[cursor]=='(')bracket_count++;
        if(formula[cursor]==')')bracket_count--;
        if(formula[cursor]=='^' && bracket_count == 0){
          has_found_symbol = true;
          break;
        } else {
          s1.push_back(formula[cursor]);
        }
        cursor++;
      }
      if(has_found_symbol){
        // split at ^
        #if defined(ENABLE_DEBUG) && defined(R_BUILD)
        Rcpp::Rcout << " Split at ^";
        #endif
        if(formula[cursor]=='^'){
          calc.instructions.push_back(Do::Power);
        }  else {
          #ifdef R_BUILD
          Rcpp::stop("Oops, something has gone wrong (f3)");
          #endif
        }
        cursor++;
        while(cursor < nchar){
          s2.push_back(formula[cursor]);
          cursor++;
        }
        s1_check = parse_formula(s1,calc,data,colnames,Xdata,bracket_flag);
        s2_check = parse_formula(s2,calc,data,colnames,Xdata,bracket_flag);
        if(s1_check || s2_check)calc.any_nonlinear = true;
        added_a_parameter = s1_check || s2_check;
      } else {
        // no pow, try brackets
        s1.clear();
        s2.clear();
        cursor=0;
        bracket_count = 0;
        while(!has_found_symbol && cursor < nchar){
          if(formula[cursor]=='('){
            has_found_symbol = true;
            break;
          } else {
            s1.push_back(formula[cursor]);
          }
          cursor++;
        }
        str token_as_str(s1.begin(),s1.end());
        if(has_found_symbol){
          #if defined(ENABLE_DEBUG) && defined(R_BUILD)
          Rcpp::Rcout << " Brackets";
          #endif
          cursor++;
          while(!(bracket_count == 0 && formula[cursor]==')') && cursor < nchar){
            s2.push_back(formula[cursor]);
            if(formula[cursor]=='(')bracket_count++;
            if(formula[cursor]==')')bracket_count--;
            cursor++;
          }
          #ifdef R_BUILD
          if(formula[cursor]!=')')Rcpp::stop("Matching bracket missing");
          #endif
          // process s1 as function (if size > 0)
          if(s1.size()>0){
            #if defined(ENABLE_DEBUG) && defined(R_BUILD)
            Rcpp::Rcout << " function";
            #endif
            if(token_as_str == "factor"){
              add_factor(s2,calc,data,colnames,Xdata);
            } else {
              bool check_fn_arg = true;
              if(token_as_str == "exp"){
                calc.instructions.push_back(Do::Exp);
              } else if(token_as_str == "log"){
                calc.instructions.push_back(Do::Log);
              } else if(token_as_str == "sqrt"){
                calc.instructions.push_back(Do::Sqrt);
              } else if(token_as_str == "sin"){
                calc.instructions.push_back(Do::Sin);
              } else if(token_as_str == "cos"){
                calc.instructions.push_back(Do::Cos);
              } else if(token_as_str == "erf"){
                calc.instructions.push_back(Do::ErrorFunc);
              } else if(token_as_str == "sign"){
                sign_fn(s2,calc,data,colnames,Xdata,1);
                check_fn_arg = false;
              } else if(token_as_str == "sign0"){
                sign_fn(s2,calc,data,colnames,Xdata,0);
                check_fn_arg = false;
              } else if(token_as_str == "twoway0"){
                two_way_fn(s2,calc,data,colnames,Xdata,0);
                check_fn_arg = false;
              } else if(token_as_str == "twoway1"){
                two_way_fn(s2,calc,data,colnames,Xdata,1);
                check_fn_arg = false;
              } else if(token_as_str == "twoway2"){
                two_way_fn(s2,calc,data,colnames,Xdata,2);
                check_fn_arg = false;
              } else {
                #ifdef R_BUILD
                Rcpp::stop("String " + token_as_str + " is not a recognised function");
                #endif
              }
              if(check_fn_arg){
                s2_check = parse_formula(s2,calc,data,colnames,Xdata,bracket_flag);
                if(s2_check) calc.any_nonlinear = true;
              }
            }
          } else {
            #if defined(ENABLE_DEBUG) && defined(R_BUILD)
            Rcpp::Rcout << " evaluate interior function";
            #endif
            // else evaluate the inside of the brackets
            s2_check = parse_formula(s2,calc,data,colnames,Xdata,true);
            if(s2_check)calc.any_nonlinear = true;
          }
        } else {
          // no brackets - now check data
          bool variable_in_data = check_data(token_as_str,calc,data,colnames,Xdata);
          if(!variable_in_data){
            // check numbers
            bool added_a_number = check_number(token_as_str,calc);
            if(!added_a_number){
              // anything else is interpreted as a parameter
              added_a_parameter = check_parameter(token_as_str,calc,bracket_flag);
            } 
          }
        }
      }
    }
  }
  
  #if defined(ENABLE_DEBUG) && defined(R_BUILD)
  Rcpp::Rcout << "\nOriginal formula: ";
  for(const auto& j: formula)Rcpp::Rcout << j;
  Rcpp::Rcout << "\ncalc:\n";
  for(const auto& i: calc.instructions)Rcpp::Rcout << instruction_str.at(i) << " ";
  Rcpp::Rcout << "\nParameter count: " << calc.parameter_count << " indexes:\n";
  glmmr::print_vec_1d<intvec>(calc.indexes);
  Rcpp::Rcout<< "\nnumbers: ";
  if(calc.user_number_count>0)for(const auto& i: calc.numbers)Rcpp::Rcout << i << " ";
  Rcpp::Rcout << "\nAny nonlinear: " << calc.any_nonlinear << "\nEND\n";
  #endif
  
  return added_a_parameter;
}


}

