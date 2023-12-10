#pragma once

#include "general.h"
#include "calculator.hpp"

namespace glmmr {

inline std::vector<Do> interpret_re(const CovFunc& fn){
  using instructs = std::vector<Do>;
  instructs B;
  switch(fn){
  case CovFunc::gr:
    B = {Do::PushParameter}; 
    break;
  case CovFunc::ar:
    B.push_back(Do::PushParameter);
    B.push_back(Do::PushCovData);
    B.push_back(Do::PushParameter);
    B.push_back(Do::Power);
    B.push_back(Do::Multiply);
    break;
  case CovFunc::fexp0:
     {
      const instructs C = {Do::Divide,Do::Negate,Do::Exp};
      B.push_back(Do::PushParameter);
      B.push_back(Do::PushCovData);
      B.insert(B.end(), C.begin(), C.end());
      break;
     }
  case CovFunc::fexp:
    {
       const instructs C = {Do::Divide,Do::Negate,Do::Exp,Do::PushParameter,Do::Multiply};  //var par here
       B.push_back(Do::PushParameter);
       B.push_back(Do::PushCovData);
       B.insert(B.end(), C.begin(), C.end());
       break;
    }
  case CovFunc::sqexp0:
    {
      const instructs C1 = {Do::PushParameter,Do::Square,Do::PushCovData,
                            Do::Square,Do::Divide,Do::Negate,Do::Exp};
      B.insert(B.end(), C1.begin(), C1.end());
      break;
    }
  case CovFunc::sqexp:
    {
      const instructs C1 = {Do::PushParameter,Do::Square};
      const instructs C2 = {Do::PushCovData,Do::Square,Do::Divide,Do::Negate,Do::Exp,
                            Do::PushParameter,Do::Multiply};
      B.insert(B.end(), C1.begin(), C1.end());
      B.insert(B.end(), C2.begin(), C2.end());
      break;
    }
  case CovFunc::bessel:
    {
      const instructs C = {Do::Divide,Do::Bessel};
      B.push_back(Do::PushParameter);
      B.push_back(Do::PushCovData);
      B.insert(B.end(), C.begin(), C.end());
      break;
    }
  case CovFunc::matern:
    {
      const instructs C1 = {Do::PushParameter,Do::Gamma,Do::PushParameter,Do::Int1,
                            Do::Subtract,Do::Int2,Do::Power,Do::Divide,
                            Do::Int2,Do::PushParameter,Do::Multiply,Do::Sqrt,Do::PushParameter};
      const instructs C2 = {Do::Divide,Do::Multiply,Do::PushParameter,Do::Power,Do::Multiply,
                            Do::PushParameter,Do::Int2,Do::PushParameter,
                            Do::Multiply,Do::Sqrt,Do::PushParameter};
      const instructs C3 = {Do::Divide,Do::Multiply,Do::BesselK,Do::Multiply};
      B.insert(B.end(), C1.begin(), C1.end());
      B.push_back(Do::PushCovData);
      B.insert(B.end(), C2.begin(), C2.end());
      B.push_back(Do::PushCovData);
      B.insert(B.end(), C3.begin(), C3.end());
      break;
    }
  case CovFunc::wend0:
    {
      const instructs C = {Do::Int1,Do::Subtract,Do::Power,Do::Multiply};
      B.push_back(Do::PushParameter);
      B.push_back(Do::PushParameter);
      B.push_back(Do::PushCovData);
      B.insert(B.end(), C.begin(), C.end());
      break;
    }
  case CovFunc::wend1:
    {
      const instructs C1 = {Do::PushParameter,Do::PushParameter,Do::Int1,Do::Add};
      const instructs C2 = {Do::Multiply,Do::Int1,Do::Add,Do::Multiply,Do::PushParameter,Do::Int1,Do::Add};
      const instructs C3 = {Do::Int1,Do::Subtract,Do::Power,Do::Multiply};
      B.insert(B.end(), C1.begin(), C1.end());
      B.push_back(Do::PushCovData);
      B.insert(B.end(), C2.begin(), C2.end());
      B.push_back(Do::PushCovData);
      B.insert(B.end(), C3.begin(), C3.end());
      break;
    }
  case CovFunc::wend2:
    {
      const instructs C1 = {Do::PushParameter,Do::Int1};
      const instructs C2 = {Do::Int2,Do::PushParameter,Do::Add,Do::Multiply,Do::Add,Do::Int3,
                            Do::Int1,Do::Subtract,Do::Int1,
                            Do::PushParameter,Do::Int2,Do::Add,Do::PushParameter,Do::Int2,
                            Do::Add,Do::Multiply,Do::Subtract,Do::Multiply};
      const instructs C3 = {Do::Multiply,Do::Multiply,Do::Add,Do::Multiply,Do::PushParameter,
                            Do::Int2,Do::Add};
      const instructs C4 = {Do::Int1,Do::Subtract,Do::Power,Do::Multiply};
      B.insert(B.end(), C1.begin(), C1.end());
      B.push_back(Do::PushCovData);
      B.insert(B.end(), C2.begin(), C2.end());
      B.push_back(Do::PushCovData);
      B.push_back(Do::PushCovData);
      B.insert(B.end(), C3.begin(), C3.end());
      B.push_back(Do::PushCovData);
      B.insert(B.end(), C4.begin(), C4.end());
      break;
    }
  case CovFunc::prodwm:
    {
      const instructs C1 = {Do::PushParameter,Do::PushParameter,Do::Gamma,Do::PushParameter,
                            Do::Int1,Do::Subtract,Do::Int2,Do::Power,Do::Divide,
                            Do::Multiply,Do::PushParameter};
      const instructs C2 = {Do::Power,Do::Multiply,Do::PushParameter};
      const instructs C3 = {Do::BesselK,Do::Multiply};
      const instructs C4 = {Do::Multiply,Do::Int10,Do::Int10,Do::Multiply,Do::Int10,Do::Int7,
                            Do::Add,Do::Add,Do::Subtract,
                            Do::Multiply,Do::Int2,Do::Int10,Do::Int1,Do::Add,Do::Divide};
      const instructs C5 = {Do::Multiply,Do::Add,Do::Int1,Do::Add,Do::Multiply};
      const instructs C6 = {Do::Int1,Do::Subtract,Do::Multiply};
      B.insert(B.end(), C1.begin(), C1.end());
      B.push_back(Do::PushCovData);
      B.insert(B.end(), C2.begin(), C2.end());
      B.push_back(Do::PushCovData);
      B.insert(B.end(), C3.begin(), C3.end());
      B.push_back(Do::PushCovData);
      B.push_back(Do::PushCovData);
      B.insert(B.end(), C4.begin(), C4.end());
      B.push_back(Do::PushCovData);
      B.insert(B.end(), C5.begin(), C5.end());
      B.push_back(Do::PushCovData);
      B.insert(B.end(), C6.begin(), C6.end());
      break;
    }
  case CovFunc::prodcb:
    {
      const instructs C1 = {Do::Power,Do::Int1,Do::Subtract,Do::Int3,Do::Negate,Do::Power,
                            Do::PushParameter,Do::Multiply};
      const instructs C2 = {Do::Pi,Do::Multiply,Do::Cos};
      const instructs C3 = {Do::Int1,Do::Subtract,Do::Multiply};
      const instructs C4 = {Do::Pi,Do::Sin,Do::Pi,Do::Int1,Do::Divide,Do::Multiply,Do::Add,
                            Do::Multiply};
      B.push_back(Do::PushParameter);
      B.push_back(Do::PushCovData);
      B.insert(B.end(), C1.begin(), C1.end());
      B.push_back(Do::PushCovData);
      B.insert(B.end(), C2.begin(), C2.end());
      B.push_back(Do::PushCovData);
      B.insert(B.end(), C3.begin(), C3.end());
      B.push_back(Do::PushCovData);
      break;
    }
  case CovFunc::prodek:
    {
      const instructs C1 = {Do::Power,Do::Negate,Do::Exp,Do::PushParameter,Do::Int2,Do::Pi};
      const instructs C2 = {Do::Multiply,Do::Multiply,Do::Int2,Do::Pi};
      const instructs C3 = {Do::Multiply,Do::Multiply,Do::Sin,Do::Divide};
      const instructs C4 = {Do::Int1,Do::Subtract,Do::Multiply,Do::Int2,Do::Pi};
      const instructs C6 = {Do::Multiply,Do::Multiply,Do::Cos,Do::Int1,Do::Subtract,
                            Do::Divide,Do::Pi,Do::Int1,Do::Divide,Do::Multiply,Do::Add,
                            Do::Multiply};
      B.push_back(Do::PushParameter);
      B.push_back(Do::PushCovData);
      B.insert(B.end(), C1.begin(), C1.end());
      B.push_back(Do::PushCovData);
      B.insert(B.end(), C2.begin(), C2.end());
      B.push_back(Do::PushCovData);
      B.insert(B.end(), C3.begin(), C3.end());
      B.push_back(Do::PushCovData);
      B.insert(B.end(), C4.begin(), C4.end());
      B.push_back(Do::PushCovData);
      B.insert(B.end(), C2.begin(), C2.end());
      B.push_back(Do::PushCovData);
      B.insert(B.end(), C6.begin(), C6.end());
      break;
    }
  case CovFunc::ar1: case CovFunc::ar0:
    B.push_back(Do::PushCovData);
    B.push_back(Do::PushParameter);
    B.push_back(Do::Power);
    break;
  case CovFunc::dist:
    B.push_back(Do::PushCovData);
    break;
  }
  return B;
}

//add in the indexes for each function
inline intvec interpret_re_par(const CovFunc& fn,
                               const int col_idx,
                               const intvec& par_idx){
  intvec B;

  auto addA = [&] (){
    B.push_back(col_idx);
  };
  
  auto addPar2 = [&] (int i){
    B.push_back(par_idx[i]);
    B.push_back(par_idx[i]);
  };
  
  
  switch(fn){
  case CovFunc::gr:
    B.push_back(par_idx[0]);
    break;
  case CovFunc::ar: 
    B.push_back(par_idx[0]);
    addA();
    B.push_back(par_idx[1]);
    break;
  case CovFunc::fexp0: case CovFunc::bessel:
    B.push_back(par_idx[0]);
    addA();
    break;
  case CovFunc::fexp:
    B.push_back(par_idx[1]);
    addA();
    B.push_back(par_idx[0]);
    break;
  case CovFunc::sqexp0:
    addPar2(0);
    addA();
    addA();
    break;
  case CovFunc::sqexp:
    addPar2(1);
    addA();
    addA();
    B.push_back(par_idx[0]);
    break;
  case CovFunc::matern:
    addPar2(0);
    addPar2(0);
    B.push_back(par_idx[1]);
    addA();
    addPar2(0);
    B.push_back(par_idx[0]);
    B.push_back(par_idx[1]);
    addA();
    break;
  case CovFunc::wend0:
    B.push_back(par_idx[0]);
    B.push_back(par_idx[1]);
    addA();
    break;
  case CovFunc::wend1:
    addPar2(0);
    B.push_back(par_idx[1]);
    addA();
    break;
  case CovFunc::wend2:
    B.push_back(par_idx[0]);
    addA();
    addPar2(1);
    addA();
    addA();
    B.push_back(par_idx[1]);
    addA();
    break;
  case CovFunc::prodwm:
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
  case CovFunc::prodcb:
    B.push_back(par_idx[1]);
    addA();
    B.push_back(par_idx[0]);
    addA();
    addA();
    addA();
    break;
  case CovFunc::prodek:
    B.push_back(par_idx[1]);
    addA();
    B.push_back(par_idx[0]);
    addA();
    addA();
    addA();
    addA();
    addA();
    break;
  case CovFunc::ar1: case CovFunc::ar0:
    addA();
    B.push_back(par_idx[0]);
    break;
  case CovFunc::dist:
    addA();
    break;
  }
  return B;
}

inline void re_linear_predictor(glmmr::calculator& calc,
                                const int Q){
  using instructs = std::vector<Do>;
  
  instructs re_seq = {Do::PushData,Do::PushParameter,Do::Multiply,Do::Add};
  for(int i = 0; i < Q; i++){
    calc.instructions.insert(calc.instructions.end(),re_seq.begin(),re_seq.end());
    calc.parameter_names.push_back("v_"+std::to_string(i));
    calc.data_names.push_back("z_"+std::to_string(i));
    calc.indexes.push_back(calc.data_count);
    calc.indexes.push_back(calc.parameter_count);
    calc.parameter_count++;
    calc.data_count++;
  }
}

inline void linear_predictor_to_link(glmmr::calculator& calc,
                                     const Link link){
  using instructs = std::vector<Do>;
  instructs out;
  instructs addzu = {Do::PushExtraData,Do::Add};
  calc.instructions.insert(calc.instructions.end(),addzu.begin(),addzu.end());
  
  switch (link) {
  case Link::logit:
    {
      out = calc.instructions;
      instructs logit_instruct = {Do::Negate,Do::Exp,Do::Int1,Do::Add,Do::Int1,Do::Divide};
      out.insert(out.end(),logit_instruct.begin(),logit_instruct.end());
      break;
    }
  case Link::loglink:
    {
      out = calc.instructions;
      out.push_back(Do::Exp);
      break;
    }
  case Link::probit:
    {
      // probit is a pain because of the error function!
      // this uses Abramowitz and Stegun approximation.
      instructs iStar = {Do::Int2,Do::Sqrt};
      iStar.insert(iStar.end(),calc.instructions.begin(),calc.instructions.end());
      iStar.push_back(Do::Divide);
      instructs M = iStar;
      instructs MStar = {Do::Constant1,Do::Multiply,Do::Int1,Do::Add,Do::Int1,Do::Divide};
      M.insert(M.end(),MStar.begin(),MStar.end());
      instructs Ltail = {Do::Power,Do::Multiply,Do::Add};
      instructs L1 = {Do::Constant2};
      L1.insert(L1.end(),M.begin(),M.end());
      L1.push_back(Do::Multiply);
      instructs L2 = {Do::Constant3,Do::Int2};
      L1.insert(L1.end(),L2.begin(),L2.end());
      L1.insert(L1.end(),M.begin(),M.end());
      L1.insert(L1.end(),Ltail.begin(),Ltail.end());
      L2 = {Do::Constant4,Do::Int3};
      L1.insert(L1.end(),L2.begin(),L2.end());
      L1.insert(L1.end(),M.begin(),M.end());
      L1.insert(L1.end(),Ltail.begin(),Ltail.end());
      L2 = {Do::Constant5,Do::Int4};
      L1.insert(L1.end(),L2.begin(),L2.end());
      L1.insert(L1.end(),M.begin(),M.end());
      L1.insert(L1.end(),Ltail.begin(),Ltail.end());
      L2 = {Do::Constant6,Do::Int5};
      L1.insert(L1.end(),L2.begin(),L2.end());
      L1.insert(L1.end(),M.begin(),M.end());
      L1.push_back(Do::Power);
      L1.push_back(Do::Multiply);
      instructs L3 = {Do::Int2};
      L3.insert(L3.end(),iStar.begin(),iStar.end());
      instructs L4 = {Do::Divide,Do::Negate,Do::Power};
      L3.insert(L3.end(),L4.begin(),L4.end());
      out = L1;
      out.insert(out.end(),L3.begin(),L3.end());
      out.push_back(Do::Multiply);
      out.push_back(Do::Int1);
      out.push_back(Do::Subtract);
      break;
    }
  case Link::identity:
    {
      out = calc.instructions;
      break;
    }
  case Link::inverse:
    {
      out = calc.instructions;
      instructs inverse_instruct = {Do::Int1,Do::Divide};
      out.insert(out.end(),inverse_instruct.begin(),inverse_instruct.end());
      break;
    }
  }
  
  calc.instructions = out;
}

inline void link_to_likelihood(glmmr::calculator& calc,
                               const Fam family){
  using instructs = std::vector<Do>;
  instructs out;
  intvec idx;
  
  switch (family){
    case Fam::gaussian:
      {
        instructs gaus_instruct = {Do::PushY,Do::Subtract,Do::Square,Do::Divide,Do::Int2,Do::Int1,
                                   Do::Divide,Do::Multiply,Do::Int2,Do::Pi,Do::Multiply,
                                   Do::Log,Do::Int2,Do::Int1,Do::Divide,Do::Multiply,Do::Add,
                                   Do::PushVariance,Do::Log,Do::Int2,Do::Int1,Do::Divide,
                                   Do::Multiply,Do::Add,Do::Negate};
        out.push_back(Do::PushVariance);
        out.insert(out.end(),calc.instructions.begin(),calc.instructions.end());
        idx.insert(idx.end(),calc.indexes.begin(),calc.indexes.end());
        out.insert(out.end(),gaus_instruct.begin(),gaus_instruct.end());
        break;
      }
    case Fam::bernoulli:
      {
        instructs binom_instruct = {Do::Log,Do::Multiply,Do::PushY,Do::Int1,Do::Subtract};
        instructs binom_instruct2 = {Do::Int1,Do::Subtract,Do::Log,Do::Multiply,Do::Add};
        out.push_back(Do::PushY);
        out.insert(out.end(),calc.instructions.begin(),calc.instructions.end());
        idx.insert(idx.end(),calc.indexes.begin(),calc.indexes.end());
        out.insert(out.end(),binom_instruct.begin(),binom_instruct.end());
        out.insert(out.end(),calc.instructions.begin(),calc.instructions.end());
        idx.insert(idx.end(),calc.indexes.begin(),calc.indexes.end());
        out.insert(out.end(),binom_instruct2.begin(),binom_instruct2.end());
        break;
      }
    case Fam::poisson:
      {
        instructs poisson_instruct = {Do::PushY,Do::LogFactorialApprox,Do::Add,Do::PushY};
        instructs poisson_instruct2 = {Do::Log,Do::Multiply,Do::Subtract};
        out.insert(out.end(),calc.instructions.begin(),calc.instructions.end());
        idx.insert(idx.end(),calc.indexes.begin(),calc.indexes.end());
        out.insert(out.end(),poisson_instruct.begin(),poisson_instruct.end());
        out.insert(out.end(),calc.instructions.begin(),calc.instructions.end());
        idx.insert(idx.end(),calc.indexes.begin(),calc.indexes.end());
        out.insert(out.end(),poisson_instruct2.begin(),poisson_instruct2.end());
        break;
      }
    case Fam::gamma:
      {
        instructs gamma_instruct = {Do::PushVariance,Do::PushY,Do::Multiply,Do::Divide};
        instructs gamma_instruct2 = {Do::Log,Do::PushVariance,Do::Log,Do::Subtract,
                                     Do::PushVariance,Do::Multiply,Do::PushY,Do::Log,
                                     Do::Int1,Do::PushVariance,
                                     Do::Subtract,Do::Multiply,Do::Add,Do::Subtract};
        out.insert(out.end(),calc.instructions.begin(),calc.instructions.end());
        idx.insert(idx.end(),calc.indexes.begin(),calc.indexes.end());
        out.insert(out.end(),gamma_instruct.begin(),gamma_instruct.end());
        out.insert(out.end(),calc.instructions.begin(),calc.instructions.end());
        idx.insert(idx.end(),calc.indexes.begin(),calc.indexes.end());
        out.insert(out.end(),gamma_instruct2.begin(),gamma_instruct2.end());
        break;
      }
    case Fam::beta:
      {
        instructs beta_instruct = {Do::PushVariance,Do::Subtract,Do::PushY,Do::Log,Do::Multiply,
                                   Do::Int1};
        instructs beta_instruct2 = {Do::Int1,Do::Subtract,Do::PushVariance,Do::Multiply,Do::Subtract,
                                    Do::PushY,Do::Int1,Do::Subtract,Do::Log,
                                    Do::Multiply,Do::Add};
        instructs beta_instruct3 = {Do::PushVariance,Do::Multiply,Do::Gamma,Do::Log,Do::Negate,Do::Add};
        instructs beta_instruct4 = {Do::Int1,Do::Subtract,Do::PushVariance,Do::Multiply,Do::Gamma,Do::Log,
                                    Do::Negate,Do::Add,Do::PushVariance,
                                    Do::Gamma,Do::Log,Do::Add};
        out.push_back(Do::Int1);
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
    case Fam::binomial:
      {
        instructs binom_instruct = {Do::PushY,Do::LogFactorialApprox,Do::PushY,Do::PushVariance,
                                    Do::Subtract,Do::Add,Do::PushVariance,
                                    Do::LogFactorialApprox,Do::Add};
        instructs binom_instruct2 = {Do::Log,Do::PushY,Do::Multiply,Do::Add};
        instructs binom_instruct3 = {Do::Int1,Do::Subtract,Do::Log,Do::PushY,Do::PushVariance,
                                     Do::Subtract,Do::Multiply,Do::Add};
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
