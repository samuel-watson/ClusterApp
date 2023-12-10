#pragma once

#include "general.h"
#include "calculator.hpp"

namespace glmmr{

inline bool parse_formula(std::vector<char>& formula,
                          glmmr::calculator& calc,
                          const ArrayXXd& data,
                          const strvec& colnames,
                          MatrixXd& Xdata,
                          bool bracket_flag = false){
  #ifdef ENABLE_DEBUG
  Rcpp::Rcout << "\nFORMULA PARSING\nFormula: ";
  for(const auto& i: formula) Rcpp::Rcout << i;
  Rcpp::Rcout << "\nAny nonlinear: " << calc.any_nonlinear;
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
    #ifdef ENABLE_DEBUG
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
      str s1_parname = "b_" + s1_as_str;
      s1.push_back('*');
      for(unsigned int j = 0; j < s1_parname.size(); j++){
        s1.push_back(s1_parname[j]);
      }
    }
    
    str s2_as_str(s2.begin(),s2.end());
    col_idx = std::find(colnames.begin(),colnames.end(),s2_as_str);
    if(col_idx != colnames.end()){
      str s2_parname = "b_" + s2_as_str;
      s2.push_back('*');
      for(unsigned int j = 0; j < s2_parname.size(); j++){
        s2.push_back(s2_parname[j]);
      }
    }
    parse_formula(s1,calc,data,colnames,Xdata,bracket_flag);
    parse_formula(s2,calc,data,colnames,Xdata,bracket_flag);
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
      #ifdef ENABLE_DEBUG
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
        #ifdef ENABLE_DEBUG
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
      } else {
        // no pow, try brackets
        s1.clear();
        s2.clear();
        cursor=0;
        bracket_count = 0;
        while(!has_found_symbol && cursor < nchar){
          //if(cursor==0 && formula[cursor]=='(')break;
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
          #ifdef ENABLE_DEBUG
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
            #ifdef ENABLE_DEBUG
            Rcpp::Rcout << " function";
            #endif
            if(token_as_str == "factor"){
              // rewrite s2 to have all the unique values of s2 variable
              // 1. check
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
                added_a_parameter = true;
              } else {
                #ifdef R_BUILD
                Rcpp::stop("Factor variable " + token_as_str + " not in data");
                #endif
              }
            } else {
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
              } else {
                #ifdef R_BUILD
                Rcpp::stop("String " + token_as_str + " is not a recognised function");
                #endif
              }
              s2_check = parse_formula(s2,calc,data,colnames,Xdata,bracket_flag);
              if(s2_check)calc.any_nonlinear = true;
            }
          } else {
            #ifdef ENABLE_DEBUG
            Rcpp::Rcout << " evaluate interior function";
            #endif
            // else evaluate the inside of the brackets
            s2_check = parse_formula(s2,calc,data,colnames,Xdata,true);
            if(s2_check)calc.any_nonlinear = true;
          }
        } else {
          // no brackets - now check data
          auto colidx = std::find(colnames.begin(),colnames.end(),token_as_str);
          if(colidx != colnames.end()){
            
            #ifdef ENABLE_DEBUG
            Rcpp::Rcout << " data";
            #endif
            // token is the name of a variable
            calc.instructions.push_back(Do::PushData);
            calc.data_names.push_back(token_as_str);
            int column_index = colidx - colnames.begin();
            calc.indexes.push_back(calc.data_count);
            if(Xdata.cols()<=calc.data_count)Xdata.conservativeResize(NoChange,calc.data_count+1);
            Xdata.col(calc.data_count) = data.col(column_index);
            calc.data_count++;
          } else if(glmmr::is_number(token_as_str)) {
            double a = std::stod(token_as_str); 
            
            #ifdef ENABLE_DEBUG
            Rcpp::Rcout << " double = " << a << " push index = " << instruction_str.at(static_cast<Do>(calc.user_number_count));
            #endif
            
            #ifdef R_BUILD
            if(calc.user_number_count >=10)Rcpp::stop("Only ten user numbers currently permitted.");
            #endif
            
            calc.instructions.push_back(static_cast<Do>(calc.user_number_count));
            calc.numbers[calc.user_number_count] = a;
            calc.user_number_count++;
            
          } else {
            #ifdef ENABLE_DEBUG
            Rcpp::Rcout << " parameter name: " << token_as_str;
            #endif
            // interpret any other string as the name of a parameter
            calc.instructions.push_back(Do::PushParameter);
            calc.parameter_names.push_back(token_as_str);
            calc.indexes.push_back(calc.parameter_count);
            calc.parameter_count++;
            added_a_parameter = true;
            if(bracket_flag)calc.any_nonlinear = true;
          }
        }
      }
    }
  }
  
  #ifdef ENABLE_DEBUG
  Rcpp::Rcout << "\nOriginal formula: ";
  for(const auto& j: formula)Rcpp::Rcout << j;
  Rcpp::Rcout << "\ncalc:\n";
  for(const auto& i: calc.instructions)Rcpp::Rcout << instruction_str.at(i) << " ";
  Rcpp::Rcout << "\nParameter count: " << calc.parameter_count << " indexes:\n";
  glmmr::print_vec_1d<intvec>(calc.indexes);
  Rcpp::Rcout<< "\nnumbers: ";
  if(calc.user_number_count>0)for(const auto& i: calc.numbers)Rcpp::Rcout << i << " ";
  Rcpp::Rcout << "\nAny nonlinear: " << calc.any_nonlinear;
  #endif
  
  return added_a_parameter;
}

}

