#ifndef INTERPRETER_H
#define INTERPRETER_H

#include "general.h"
#include "calculator.hpp"

namespace glmmr {

inline intvec interpret_re(const std::string& fn,
                           const intvec& A){
  intvec B;
  switch(string_to_case.at(fn)){
  case 1:
    B = {2}; //var par here
    break;
  case 2:
    B.push_back(2);
    B.insert(B.end(), A.begin(), A.end());
    B.push_back(2);
    B.push_back(8);
    B.push_back(5);
    break;
  case 3:
     {
      const intvec C = {6,10,9};
      B.push_back(2);
      B.insert(B.end(), A.begin(), A.end());
      B.insert(B.end(), C.begin(), C.end());
      break;
     }
  case 4:
    {
       const intvec C = {6,10,9,2,5};  //var par here
       B.push_back(2);
       B.insert(B.end(), A.begin(), A.end());
       B.insert(B.end(), C.begin(), C.end());
       break;
    }
  case 5:
    {
      const intvec C1 = {2,2,5};
      const intvec C2 = {5,6,10,9};
      B.insert(B.end(), C1.begin(), C1.end());
      B.insert(B.end(), A.begin(), A.end());
      B.insert(B.end(), A.begin(), A.end());
      B.insert(B.end(), C2.begin(), C2.end());
      break;
    }
  case 6:
    {
      const intvec C1 = {2,2,5};
      const intvec C2 = {5,6,10,9,2,5};//var par here
      B.insert(B.end(), C1.begin(), C1.end());
      B.insert(B.end(), A.begin(), A.end());
      B.insert(B.end(), A.begin(), A.end());
      B.insert(B.end(), C2.begin(), C2.end());
      break;
    }
  case 7:
    {
      const intvec C = {6,11};
      B.push_back(2);
      B.insert(B.end(), A.begin(), A.end());
      B.insert(B.end(), C.begin(), C.end());
      break;
    }
  case 8:
    {
      const intvec C1 = {2,12,2,21,4,22,8,6,22,2,5,7,2};
      const intvec C2 = {6,5,2,8,5,2,22,2,5,7,2};
      const intvec C3 = {6,5,15,5};
      B.insert(B.end(), C1.begin(), C1.end());
      B.insert(B.end(), A.begin(), A.end());
      B.insert(B.end(), C2.begin(), C2.end());
      B.insert(B.end(), A.begin(), A.end());
      B.insert(B.end(), C3.begin(), C3.end());
      break;
    }
  case 9:
    {
      const intvec C = {21,4,8,5};
      B.push_back(2);
      B.push_back(2);
      B.insert(B.end(), A.begin(), A.end());
      B.insert(B.end(), C.begin(), C.end());
      break;
    }
  case 10:
    {
      const intvec C1 = {2,2,21,3};
      const intvec C2 = {5,21,3,5,2,21,3};
      const intvec C3 = {21,4,8,5};
      B.insert(B.end(), C1.begin(), C1.end());
      B.insert(B.end(), A.begin(), A.end());
      B.insert(B.end(), C2.begin(), C2.end());
      B.insert(B.end(), A.begin(), A.end());
      B.insert(B.end(), C3.begin(), C3.end());
      break;
    }
  case 11:
    {
      const intvec C1 = {2,21};
      const intvec C2 = {22,2,3,5,3,23,21,4,21,2,22,3,2,22,3,5,4,5};
      const intvec C3 = {5,5,3,5,2,22,3};
      const intvec C4 = {21,4,8,5};
      B.insert(B.end(), C1.begin(), C1.end());
      B.insert(B.end(), A.begin(), A.end());
      B.insert(B.end(), C2.begin(), C2.end());
      B.insert(B.end(), A.begin(), A.end());
      B.insert(B.end(), A.begin(), A.end());
      B.insert(B.end(), C3.begin(), C3.end());
      B.insert(B.end(), A.begin(), A.end());
      B.insert(B.end(), C4.begin(), C4.end());
      break;
    }
  case 12:
    {
      const intvec C1 = {2,2,12,2,21,4,22,8,6,5,2};
      const intvec C2 = {8,5,2};
      const intvec C3 = {15,5};
      const intvec C4 = {5,20,20,5,20,27,3,3,4,5,22,20,21,3,6};
      const intvec C5 = {5,3,21,3,5};
      const intvec C6 = {21,4,5};
      B.insert(B.end(), C1.begin(), C1.end());
      B.insert(B.end(), A.begin(), A.end());
      B.insert(B.end(), C2.begin(), C2.end());
      B.insert(B.end(), A.begin(), A.end());
      B.insert(B.end(), C3.begin(), C3.end());
      B.insert(B.end(), A.begin(), A.end());
      B.insert(B.end(), A.begin(), A.end());
      B.insert(B.end(), C4.begin(), C4.end());
      B.insert(B.end(), A.begin(), A.end());
      B.insert(B.end(), C5.begin(), C5.end());
      B.insert(B.end(), A.begin(), A.end());
      B.insert(B.end(), C6.begin(), C6.end());
      break;
    }
  case 13:
    {
      const intvec C1 = {8,21,4,23,10,8,2,5};
      const intvec C2 = {30,5,14};
      const intvec C3 = {21,4,5};
      const intvec C4 = {30,13,30,21,6,5,3,5};
      B.push_back(2);
      B.insert(B.end(), A.begin(), A.end());
      B.insert(B.end(), C1.begin(), C1.end());
      B.insert(B.end(), A.begin(), A.end());
      B.insert(B.end(), C2.begin(), C2.end());
      B.insert(B.end(), A.begin(), A.end());
      B.insert(B.end(), C3.begin(), C3.end());
      B.insert(B.end(), A.begin(), A.end());
      break;
    }
  case 14:
    {
      const intvec C1 = {8,10,9,2,22,30};
      const intvec C2 = {5,5,22,30};
      const intvec C3 = {5,5,13,6};
      const intvec C4 = {21,4,5,22,30};
      const intvec C5 = {5,5,22,30};
      const intvec C6 = {5,5,14,21,4,6,30,21,6,5,3,5};
      B.push_back(2);
      B.insert(B.end(), A.begin(), A.end());
      B.insert(B.end(), C1.begin(), C1.end());
      B.insert(B.end(), A.begin(), A.end());
      B.insert(B.end(), C2.begin(), C2.end());
      B.insert(B.end(), A.begin(), A.end());
      B.insert(B.end(), C3.begin(), C3.end());
      B.insert(B.end(), A.begin(), A.end());
      B.insert(B.end(), C4.begin(), C4.end());
      B.insert(B.end(), A.begin(), A.end());
      B.insert(B.end(), C5.begin(), C5.end());
      B.insert(B.end(), A.begin(), A.end());
      B.insert(B.end(), C6.begin(), C6.end());
      break;
    }
  case 15: case 16:
    B.insert(B.end(), A.begin(), A.end());
    B.push_back(2);
    B.push_back(8);
    break;
  }
  return B;
}

//add in the indexes for each function
inline intvec interpret_re_par(const std::string& fn,
                               const intvec& col_idx,
                               const intvec& par_idx){
  intvec B;
  
  
  auto addA = [&] (){
    for(int i = 0; i<col_idx.size();i++){
      B.push_back(col_idx[i]);
      B.push_back(col_idx[i]);
      B.push_back(col_idx[i]);
      B.push_back(col_idx[i]);
    }
  };
  
  auto addPar2 = [&] (int i){
    B.push_back(par_idx[i]);
    B.push_back(par_idx[i]);
  };
  
  
  switch(string_to_case.at(fn)){
  case 1:
    //addPar2(0);
    B.push_back(par_idx[0]);
    break;
  case 2: 
    B.push_back(par_idx[0]);
    addA();
    B.push_back(par_idx[1]);
    break;
  case 3: case 7:
    B.push_back(par_idx[0]);
    addA();
    break;
  case 4:
    B.push_back(par_idx[1]);
    addA();
    //addPar2(0);
    B.push_back(par_idx[0]);
    break;
  case 5:
    addPar2(0);
    addA();
    addA();
    break;
  case 6:
    addPar2(1);
    addA();
    addA();
    //addPar2(0);
    B.push_back(par_idx[0]);
    break;
  case 8:
    addPar2(0);
    addPar2(0);
    B.push_back(par_idx[1]);
    addA();
    addPar2(0);
    B.push_back(par_idx[0]);
    B.push_back(par_idx[1]);
    addA();
    break;
  case 9:
    B.push_back(par_idx[0]);
    B.push_back(par_idx[1]);
    addA();
    break;
  case 10:
    addPar2(0);
    B.push_back(par_idx[1]);
    addA();
    break;
  case 11:
    B.push_back(par_idx[0]);
    addA();
    addPar2(1);
    addA();
    addA();
    B.push_back(par_idx[1]);
    addA();
    break;
  case 12:
    B.push_back(par_idx[0]);
    addPar2(1);
    addPar2(0);
    addA();
    addA();
    addA();
    addA();
    addA();
    addA();
    break;
  case 13:
    B.push_back(par_idx[1]);
    addA();
    B.push_back(par_idx[0]);
    addA();
    addA();
    addA();
    break;
  case 14:
    B.push_back(par_idx[1]);
    addA();
    B.push_back(par_idx[0]);
    addA();
    addA();
    addA();
    addA();
    addA();
    break;
  case 15: case 16:
    addA();
    B.push_back(par_idx[0]);
    break;
  }
  return B;
}

inline void re_linear_predictor(glmmr::calculator& calc,
                                const int& Q){
  
  intvec re_instruct;
  intvec re_seq = {0,2,5,3};
  for(int i = 0; i < Q; i++){
    re_instruct.insert(re_instruct.end(),re_seq.begin(),re_seq.end());
    calc.parameter_names.push_back("v_"+std::to_string(i));
    calc.indexes.push_back(i+calc.data_count);
    calc.indexes.push_back(i+calc.data_count);
  }
  calc.parameter_count += Q;
  calc.instructions.insert(calc.instructions.end(),re_instruct.begin(),re_instruct.end());
  calc.data_count += Q;
}

inline void linear_predictor_to_link(glmmr::calculator& calc,
                                     const str& link){
  intvec out;
  intvec addzu = {18,3};
  calc.instructions.insert(calc.instructions.end(),addzu.begin(),addzu.end());
  const static std::unordered_map<std::string, int> link_to_case{
    {"logit",1},
    {"log",2},
    {"probit",3},
    {"identity",4},
    {"inverse",5}
  };
  switch (link_to_case.at(link)) {
  case 1:
    {
      out = calc.instructions;
      intvec logit_instruct = {10,9,21,3,21,6};
      out.insert(out.end(),logit_instruct.begin(),logit_instruct.end());
      break;
    }
  case 2:
    {
      out = calc.instructions;
      out.push_back(9);
      break;
    }
  case 3:
    {
      // probit is a pain in the ass because of the error function!
      // this uses Abramowitz and Stegun approximation.
      intvec iStar = {22,7};
      iStar.insert(iStar.end(),calc.instructions.begin(),calc.instructions.end());
      iStar.push_back(6);
      intvec M = iStar;
      intvec MStar = {31,5,21,3,21,6};
      M.insert(M.end(),MStar.begin(),MStar.end());
      intvec Ltail = {8,5,3};
      intvec L1 = {32};
      L1.insert(L1.end(),M.begin(),M.end());
      L1.push_back(5);
      intvec L2 = {33,22};
      L1.insert(L1.end(),L2.begin(),L2.end());
      L1.insert(L1.end(),M.begin(),M.end());
      L1.insert(L1.end(),Ltail.begin(),Ltail.end());
      L2 = {34,23};
      L1.insert(L1.end(),L2.begin(),L2.end());
      L1.insert(L1.end(),M.begin(),M.end());
      L1.insert(L1.end(),Ltail.begin(),Ltail.end());
      L2 = {35,24};
      L1.insert(L1.end(),L2.begin(),L2.end());
      L1.insert(L1.end(),M.begin(),M.end());
      L1.insert(L1.end(),Ltail.begin(),Ltail.end());
      L2 = {36,25};
      L1.insert(L1.end(),L2.begin(),L2.end());
      L1.insert(L1.end(),M.begin(),M.end());
      L1.push_back(8);
      L1.push_back(5);
      intvec L3 = {22};
      L3.insert(L3.end(),iStar.begin(),iStar.end());
      intvec L4 = {6,10,8};
      L3.insert(L3.end(),L4.begin(),L4.end());
      out = L1;
      out.insert(out.end(),L3.begin(),L3.end());
      out.push_back(5);
      out.push_back(21);
      out.push_back(4);
      break;
    }
  case 4:
    {
      out = calc.instructions;
      break;
    }
  case 5:
    {
      out = calc.instructions;
      intvec inverse_instruct = {21,6};
      out.insert(out.end(),inverse_instruct.begin(),inverse_instruct.end());
      break;
    }
  }
  
  calc.instructions = out;
}

inline void link_to_likelihood(glmmr::calculator& calc,
                               const str& family){
  
  intvec out;
  intvec idx;
  const static std::unordered_map<std::string, int> family_to_case{
    {"gaussian",1},
    {"bernoulli",2},
    {"poisson",3},
    {"gamma",4},
    {"beta",5},
    {"binomial",6}
  };
  
  
  switch (family_to_case.at(family)){
    case 1:
      {
        intvec gaus_instruct = {19,4,17,6,22,21,6,5,22,30,5,16,22,21,6,5,3,41,16,22,21,6,5,3,10};
        out.push_back(41);
        out.insert(out.end(),calc.instructions.begin(),calc.instructions.end());
        idx.insert(idx.end(),calc.indexes.begin(),calc.indexes.end());
        out.insert(out.end(),gaus_instruct.begin(),gaus_instruct.end());
        break;
      }
    case 2:
      {
        intvec binom_instruct = {16,5,19,21,4};
        intvec binom_instruct2 = {21,4,16,5,3};
        out.push_back(19);
        out.insert(out.end(),calc.instructions.begin(),calc.instructions.end());
        idx.insert(idx.end(),calc.indexes.begin(),calc.indexes.end());
        out.insert(out.end(),binom_instruct.begin(),binom_instruct.end());
        out.insert(out.end(),calc.instructions.begin(),calc.instructions.end());
        idx.insert(idx.end(),calc.indexes.begin(),calc.indexes.end());
        out.insert(out.end(),binom_instruct2.begin(),binom_instruct2.end());
        break;
      }
    case 3:
      {
        intvec poisson_instruct = {19,40,3,19};
        intvec poisson_instruct2 = {16,5,4};
        out.insert(out.end(),calc.instructions.begin(),calc.instructions.end());
        idx.insert(idx.end(),calc.indexes.begin(),calc.indexes.end());
        out.insert(out.end(),poisson_instruct.begin(),poisson_instruct.end());
        out.insert(out.end(),calc.instructions.begin(),calc.instructions.end());
        idx.insert(idx.end(),calc.indexes.begin(),calc.indexes.end());
        out.insert(out.end(),poisson_instruct2.begin(),poisson_instruct2.end());
        break;
      }
    case 4:
      {
        intvec gamma_instruct = {41,19,5,6};
        intvec gamma_instruct2 = {41,19,5,6,16,41,5,4,41,12,19,5,21,6,1,6,3};
        out.insert(out.end(),calc.instructions.begin(),calc.instructions.end());
        idx.insert(idx.end(),calc.indexes.begin(),calc.indexes.end());
        out.insert(out.end(),gamma_instruct.begin(),gamma_instruct.end());
        out.insert(out.end(),calc.instructions.begin(),calc.instructions.end());
        idx.insert(idx.end(),calc.indexes.begin(),calc.indexes.end());
        out.insert(out.end(),gamma_instruct2.begin(),gamma_instruct2.end());
        break;
      }
    case 5:
      {
        intvec beta_instruct = {41,4,19,16,5,21};
        intvec beta_instruct2 = {21,4,41,5,4,19,21,4,16,5,3};
        intvec beta_instruct3 = {41,5,12,16,10,3};
        intvec beta_instruct4 = {21,4,41,5,12,16,10,3,41,12,16,3};
        out.push_back(21);
        out.insert(out.end(),calc.instructions.begin(),calc.instructions.end());
        idx.insert(idx.end(),calc.indexes.begin(),calc.indexes.end());
        out.insert(out.end(),beta_instruct.begin(),beta_instruct.end());
        out.insert(out.end(),calc.instructions.begin(),calc.instructions.end());
        idx.insert(idx.end(),calc.indexes.begin(),calc.indexes.end());
        out.insert(out.end(),beta_instruct2.begin(),beta_instruct2.end());
        out.insert(out.end(),calc.instructions.begin(),calc.instructions.end());
        idx.insert(idx.end(),calc.indexes.begin(),calc.indexes.end());
        out.insert(out.end(),beta_instruct3.begin(),beta_instruct3.end());
        out.insert(out.end(),calc.instructions.begin(),calc.instructions.end());
        idx.insert(idx.end(),calc.indexes.begin(),calc.indexes.end());
        out.insert(out.end(),beta_instruct4.begin(),beta_instruct4.end());
        break;
      }
    case 6:
      {
        intvec binom_instruct = {19,40,19,41,4,3,41,40,3};
        intvec binom_instruct2 = {16,19,5,3};
        intvec binom_instruct3 = {21,4,16,19,41,4,5,3};
        out.insert(out.end(),binom_instruct.begin(),binom_instruct.end());
        out.insert(out.end(),calc.instructions.begin(),calc.instructions.end());
        idx.insert(idx.end(),calc.indexes.begin(),calc.indexes.end());
        out.insert(out.end(),binom_instruct2.begin(),binom_instruct2.end());
        out.insert(out.end(),calc.instructions.begin(),calc.instructions.end());
        idx.insert(idx.end(),calc.indexes.begin(),calc.indexes.end());
        out.insert(out.end(),binom_instruct3.begin(),binom_instruct3.end());
      }
  }
  calc.instructions = out;
  calc.indexes = idx;
}

}

#endif