#ifndef FORMULAPARSE_H
#define FORMULAPARSE_H

#include "general.h"

namespace glmmr{

inline bool parse_formula(std::vector<char>& formula,
                          glmmr::calculator& calc,
                          const ArrayXXd& data,
                          const strvec& colnames,
                          MatrixXd& Xdata){
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
    if(formula[cursor]=='(')bracket_count++;
    if(formula[cursor]==')')bracket_count--;
    if((formula[cursor]=='+' || formula[cursor]=='-') && bracket_count == 0){
      has_found_symbol = true;
      break;
    } else {
      s1.push_back(formula[cursor]);
    }
    cursor++;
  }
  if(has_found_symbol){
    // split at +/-
    if(formula[cursor]=='+'){
      calc.instructions.push_back(3);
    } else if(formula[cursor]=='-'){
      calc.instructions.push_back(4);
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
    parse_formula(s1,calc,data,colnames,Xdata);
    parse_formula(s2,calc,data,colnames,Xdata);
  } else {
    // no +/- to split at, try *//
    s1.clear();
    s2.clear();
    cursor=0;
    bracket_count = 0;
    has_found_symbol = false;
    while(!has_found_symbol && cursor < nchar){
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
      // split at *//
      if(formula[cursor]=='*'){
        calc.instructions.push_back(5);
      } else if(formula[cursor]=='/'){
        calc.instructions.push_back(6);
      } 
      
      cursor++;
      while(cursor < nchar){
        s2.push_back(formula[cursor]);
        cursor++;
      }
      
      s1_check = parse_formula(s1,calc,data,colnames,Xdata);
      s2_check = parse_formula(s2,calc,data,colnames,Xdata);
      if(s1_check && s2_check)calc.any_nonlinear = true;
    } else {
      // no * to split at, try pow
      s1.clear();
      s2.clear();
      cursor=0;
      bracket_count = 0;
      while(!has_found_symbol && cursor < nchar){
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
        if(formula[cursor]=='^'){
          calc.instructions.push_back(8);
        } 
        cursor++;
        while(cursor < nchar){
          s2.push_back(formula[cursor]);
          cursor++;
        }
        
        s1_check = parse_formula(s1,calc,data,colnames,Xdata);
        s2_check = parse_formula(s2,calc,data,colnames,Xdata);
        if(s1_check || s2_check)calc.any_nonlinear = true;
      } else {
        // no pow, try brackets
        s1.clear();
        s2.clear();
        cursor=0;
        bracket_count = 0;
        while(!has_found_symbol && cursor < nchar){
          if(cursor==0 && formula[cursor]=='(')break;
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
          cursor++;
          while(!(bracket_count == 0 && formula[cursor]==')') && cursor < nchar){
            s2.push_back(formula[cursor]);
            if(formula[cursor]=='(')bracket_count++;
            if(formula[cursor]==')')bracket_count--;
            cursor++;
          }
          // process s1 as function (if size > 0)
          if(s1.size()>0){
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
                //bool itemsinxb = calc.instructions.size() > 0;
                for(int i = 0; i < factorrange; i++){
                  if(i < (factorrange - 1))calc.instructions.push_back(3);
                  calc.instructions.push_back(5);
                  if(Xdata.cols()<=calc.data_count)Xdata.conservativeResize(NoChange,calc.data_count+1);
                  for(int j = 0; j < data.rows(); j++){
                    //calc.data[j].push_back(data(j,column_index)==unique_values[i] ? 1.0 : 0.0);
                    Xdata(j,calc.data_count) = data(j,column_index)==unique_values[i] ? 1.0 : 0.0;
                  }
                  calc.indexes.push_back(calc.data_count);
                  calc.data_count++;
                  calc.instructions.push_back(0);
                  // parameter
                  calc.instructions.push_back(2);
                  str parname = "b_" + token_as_str2 + "_" + std::to_string(unique_values[i])[0];
                  calc.parameter_names.push_back(parname);
                  calc.indexes.push_back(calc.parameter_count);
                  calc.parameter_count++;
                }
                added_a_parameter = true;
              } 
            } else {
              if(token_as_str == "exp"){
                calc.instructions.push_back(9);
              } else if(token_as_str == "log"){
                calc.instructions.push_back(16);
              } else if(token_as_str == "sqrt"){
                calc.instructions.push_back(7);
              } else if(token_as_str == "sin"){
                calc.instructions.push_back(13);
              } else if(token_as_str == "cos"){
                calc.instructions.push_back(14);
              } 
              
              s2_check = parse_formula(s2,calc,data,colnames,Xdata);
              if(s2_check)calc.any_nonlinear = true;
            }
          }
        } else {
          // no brackets - now check data
          auto colidx = std::find(colnames.begin(),colnames.end(),token_as_str);
          if(colidx != colnames.end()){
            // token is the name of a variable
            calc.instructions.push_back(0);
            int column_index = colidx - colnames.begin();
            calc.indexes.push_back(calc.data_count);
            if(Xdata.cols()<=calc.data_count)Xdata.conservativeResize(NoChange,calc.data_count+1);
            Xdata.col(calc.data_count) = data.col(column_index);
            calc.data_count++;
          } else if(glmmr::is_number(token_as_str)) {
            // add an integer to the stack
            int p = s1.size();
            int addint;
            if(p > 1){
              for(int i = 0; i < (p-1); i++){
                calc.instructions.push_back(3);
              }
            }
            for(int k = 1; k <= p; k++){
              int number = s1[p-k] - '0';
              addint = number + 20;
              if(k==1){
                calc.instructions.push_back(addint);
              } else {
                calc.instructions.push_back(5);
                calc.instructions.push_back(addint);
                calc.instructions.push_back(6);
                calc.instructions.push_back(k-1);
                calc.instructions.push_back(20);
              }
            }
          } else {
            // interpret any other string as the name of a parameter
            calc.instructions.push_back(2);
            calc.parameter_names.push_back(token_as_str);
            calc.indexes.push_back(calc.parameter_count);
            calc.parameter_count++;
            added_a_parameter = true;
          }
        }
      }
    }
  }
  return added_a_parameter;
}

}



#endif