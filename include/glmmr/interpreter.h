#pragma once

#include "general.h"
#include "calculator.hpp"

namespace glmmr {


enum class CovFunc {
  gr = 0,
    ar = 1,
    fexp0 = 2,
    fexp = 3,
    sqexp0 = 4,
    sqexp = 5,
    bessel = 6,
    matern = 7,
    truncpow2 = 8,
    truncpow3 = 9,
    truncpow4 = 10,
    cauchy = 11,
    cauchy3 = 12,
    truncpow20 = 13,
    truncpow30 = 14,
    truncpow40 = 15,
    cauchy0 = 16,
    cauchy30 = 17,
    ar0 = 18,
    ar1 = 19,
    dist = 20
};

const std::map<str, CovFunc> str_to_covfunc = {
  {"gr", CovFunc::gr},
  {"ar", CovFunc::ar},
  {"fexp0", CovFunc::fexp0},
  {"fexp", CovFunc::fexp},
  {"sqexp0",CovFunc::sqexp0},
  {"sqexp",CovFunc::sqexp},
  {"bessel",CovFunc::bessel},
  {"matern",CovFunc::matern},
  {"truncpow2",CovFunc::truncpow2},
  {"truncpow3",CovFunc::truncpow3},
  {"truncpow4",CovFunc::truncpow4},
  {"cauchy",CovFunc::cauchy},
  {"cauchy3",CovFunc::cauchy3},
  {"truncpow20",CovFunc::truncpow20},
  {"truncpow30",CovFunc::truncpow30},
  {"truncpow40",CovFunc::truncpow40},
  {"cauchy0",CovFunc::cauchy0},
  {"cauchy30",CovFunc::cauchy30},
  {"ar0", CovFunc::ar0},
  {"ar1", CovFunc::ar1},
  {"dist",CovFunc::dist}
};

// unfortunately need bidirectional map so need to duplicate this unless there's
// a better way??
const std::map<CovFunc, str> covfunc_to_str = {
  {CovFunc::gr, "gr"},
  {CovFunc::ar, "ar"},
  {CovFunc::fexp0, "fexp0"},
  {CovFunc::fexp, "fexp"},
  {CovFunc::sqexp0, "sqexp0"},
  {CovFunc::sqexp, "sqexp"},
  {CovFunc::bessel, "bessel"},
  {CovFunc::matern, "matern"},
  {CovFunc::truncpow2, "truncpow2"},
  {CovFunc::truncpow3, "truncpow3"},
  {CovFunc::truncpow4, "truncpow4"},
  {CovFunc::cauchy, "cauchy"},
  {CovFunc::cauchy3, "cauchy3"},
  {CovFunc::truncpow20, "truncpow20"},
  {CovFunc::truncpow30, "truncpow30"},
  {CovFunc::truncpow40, "truncpow40"},
  {CovFunc::cauchy0, "cauchy0"},
  {CovFunc::cauchy30, "cauchy30"},
  {CovFunc::ar0, "ar0"},
  {CovFunc::ar1, "ar1"},
  {CovFunc::dist, "dist"}
};

const std::map<CovFunc, int> covfunc_to_nvar = {
  {CovFunc::gr, 1},
  {CovFunc::ar, 2},
  {CovFunc::fexp0, 1},
  {CovFunc::fexp, 2},
  {CovFunc::sqexp0, 1},
  {CovFunc::sqexp, 2},
  {CovFunc::bessel, 1},
  {CovFunc::matern, 2},
  {CovFunc::truncpow2, 2},
  {CovFunc::truncpow3, 2},
  {CovFunc::truncpow4, 2},
  {CovFunc::cauchy, 3},
  {CovFunc::cauchy3, 2},
  {CovFunc::truncpow20, 1},
  {CovFunc::truncpow30, 1},
  {CovFunc::truncpow40, 1},
  {CovFunc::cauchy0, 2},
  {CovFunc::cauchy30, 1},
  {CovFunc::ar0, 1},
  {CovFunc::ar1, 1},
  {CovFunc::dist, 0}
};

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
  case CovFunc::truncpow2:
    {
      const instructs C = {Do::PushParameter,Do::Int2,Do::PushParameter,Do::PushCovData,Do::Power,Do::Int1,Do::Subtract,Do::Power,Do::Multiply};
      B.insert(B.end(), C.begin(), C.end());
      break;
    }
  case CovFunc::truncpow3:
  {
    const instructs C = {Do::PushParameter,Do::Int3,Do::PushParameter,Do::PushCovData,Do::Power,Do::Int1,Do::Subtract,Do::Power,Do::Multiply};
    B.insert(B.end(), C.begin(), C.end());
    break;
  }
  case CovFunc::truncpow4:
  {
    const instructs C = {Do::PushParameter,Do::Int4,Do::PushParameter,Do::PushCovData,Do::Power,Do::Int1,Do::Subtract,Do::Power,Do::Multiply};
    B.insert(B.end(), C.begin(), C.end());
    break;
  }
  case CovFunc::cauchy3:
  {
    const instructs C = {Do::PushParameter,Do::Int3,Do::Negate,Do::PushParameter,Do::PushCovData,Do::Power,Do::Int1,Do::Add,Do::Power,Do::Multiply};
    B.insert(B.end(), C.begin(), C.end());
    break;
  }
  case CovFunc::cauchy:
  {
    const instructs C = {Do::PushParameter,Do::PushParameter,Do::PushParameter,Do::Divide,Do::Negate,Do::PushParameter,Do::PushCovData,Do::Power,
                         Do::Int1,Do::Add,Do::Power,Do::Multiply};
    B.insert(B.end(), C.begin(), C.end());
    break;
  }
  case CovFunc::truncpow20:
  {
    const instructs C = {Do::Int2,Do::PushParameter,Do::PushCovData,Do::Power,Do::Int1,Do::Subtract,Do::Power};
    B.insert(B.end(), C.begin(), C.end());
    break;
  }
  case CovFunc::truncpow30:
  {
    const instructs C = {Do::Int3,Do::PushParameter,Do::PushCovData,Do::Power,Do::Int1,Do::Subtract,Do::Power};
    B.insert(B.end(), C.begin(), C.end());
    break;
  }
  case CovFunc::truncpow40:
  {
    const instructs C = {Do::Int4,Do::PushParameter,Do::PushCovData,Do::Power,Do::Int1,Do::Subtract,Do::Power};
    B.insert(B.end(), C.begin(), C.end());
    break;
  }
  case CovFunc::cauchy30:
  {
    const instructs C = {Do::Int3,Do::Negate,Do::PushParameter,Do::PushCovData,Do::Power,Do::Int1,Do::Add,Do::Power};
    B.insert(B.end(), C.begin(), C.end());
    break;
  }
  case CovFunc::cauchy0:
  {
    const instructs C = {Do::PushParameter,Do::PushParameter,Do::Divide,Do::Negate,Do::PushParameter,Do::PushCovData,Do::Power,
                         Do::Int1,Do::Add,Do::Power};
    B.insert(B.end(), C.begin(), C.end());
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
  case CovFunc::fexp0: case CovFunc::bessel: case CovFunc::sqexp0:
    B.push_back(par_idx[0]);
    addA();
    break;
  case CovFunc::fexp: case CovFunc::sqexp:
    B.push_back(par_idx[1]);
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
  case CovFunc::truncpow2: case CovFunc::truncpow3: case CovFunc::truncpow4: case CovFunc::cauchy3:
    B.push_back(par_idx[0]);
    B.push_back(par_idx[1]);
    addA();
    break;
  case CovFunc::cauchy:
    B.push_back(par_idx[0]);
    B.push_back(par_idx[1]);
    B.push_back(par_idx[2]);
    B.push_back(par_idx[1]);
    addA();
    break;
  case CovFunc::truncpow20: case CovFunc::truncpow30: case CovFunc::truncpow40: case CovFunc::cauchy30:
    B.push_back(par_idx[0]);
    addA();
    break;
  case CovFunc::cauchy0:
    B.push_back(par_idx[0]);
    B.push_back(par_idx[1]);
    B.push_back(par_idx[0]);
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

inline void re_log_likelihood(glmmr::calculator& calc,
                                const int Q){
  using instructs = std::vector<Do>;
  instructs re_seq = {Do::PushParameter,Do::Square,Do::Add};
  for(int i = 0; i < Q; i++){
    calc.instructions.insert(calc.instructions.end(),re_seq.begin(),re_seq.end());
    auto uidx = std::find(calc.parameter_names.begin(),calc.parameter_names.end(),"v_"+std::to_string(i));
    if(uidx == calc.parameter_names.end()){
      throw std::runtime_error("Error finding name of random effect in calculator");
    } else {
      int idx_to_add = uidx - calc.parameter_names.begin();
      calc.indexes.push_back(idx_to_add);
    }
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
      out.push_back(Do::SqrtTwo);
       out = calc.instructions;
       instructs probit_instruct = {Do::Divide, Do::ErrorFunc, Do::Int1, Do::Add, Do::Half, Do::Multiply};
      out.insert(out.end(),probit_instruct.begin(),probit_instruct.end());
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

// // many of these could be optimised better!
// inline void link_to_likelihood(glmmr::calculator& calc,
//                                const Fam family){
//   using instructs = std::vector<Do>;
//   instructs out;
//   intvec idx;
//   
//   switch (family){
//     case Fam::gaussian:
//       {
//         instructs gaus_instruct = {Do::PushY,Do::Subtract,Do::Square,Do::Divide,Do::Half,
//                                    Do::Multiply,Do::HalfLog2Pi,Do::Add,
//                                    Do::PushVariance,Do::Log,Do::Half,
//                                    Do::Multiply,Do::Add,Do::Negate};
//         out.push_back(Do::PushVariance);
//         out.insert(out.end(),calc.instructions.begin(),calc.instructions.end());
//         idx.insert(idx.end(),calc.indexes.begin(),calc.indexes.end());
//         out.insert(out.end(),gaus_instruct.begin(),gaus_instruct.end());
//         break;
//       }
//     case Fam::bernoulli:
//       {
//         instructs binom_instruct = {Do::Log,Do::Multiply,Do::PushY,Do::Int1,Do::Subtract};
//         instructs binom_instruct2 = {Do::Int1,Do::Subtract,Do::Log,Do::Multiply,Do::Add};
//         out.push_back(Do::PushY);
//         out.insert(out.end(),calc.instructions.begin(),calc.instructions.end());
//         idx.insert(idx.end(),calc.indexes.begin(),calc.indexes.end());
//         out.insert(out.end(),binom_instruct.begin(),binom_instruct.end());
//         out.insert(out.end(),calc.instructions.begin(),calc.instructions.end());
//         idx.insert(idx.end(),calc.indexes.begin(),calc.indexes.end());
//         out.insert(out.end(),binom_instruct2.begin(),binom_instruct2.end());
//         break;
//       }
//     case Fam::poisson:
//       {
//         out.push_back(Do::PushY);
//         out.push_back(Do::LogFactorialApprox);
//         out.insert(out.end(),calc.instructions.begin(),calc.instructions.end());
//         idx.insert(idx.end(),calc.indexes.begin(),calc.indexes.end());
//         out.insert(out.end(),calc.instructions.begin(),calc.instructions.end());
//         idx.insert(idx.end(),calc.indexes.begin(),calc.indexes.end());
//         out.push_back(Do::Log);
//         out.push_back(Do::PushY);
//         out.push_back(Do::Multiply);
//         out.push_back(Do::Subtract);
//         out.push_back(Do::Subtract);
//         break;
//       }
//     case Fam::gamma:
//       {
//         instructs gamma_instruct = {Do::PushVariance,Do::PushY,Do::Multiply,Do::Divide};
//         instructs gamma_instruct2 = {Do::Log,Do::PushVariance,Do::Log,Do::Subtract,
//                                      Do::PushVariance,Do::Multiply,Do::PushY,Do::Log,
//                                      Do::Int1,Do::PushVariance,
//                                      Do::Subtract,Do::Multiply,Do::Add,Do::Subtract};
//         out.insert(out.end(),calc.instructions.begin(),calc.instructions.end());
//         idx.insert(idx.end(),calc.indexes.begin(),calc.indexes.end());
//         out.insert(out.end(),gamma_instruct.begin(),gamma_instruct.end());
//         out.insert(out.end(),calc.instructions.begin(),calc.instructions.end());
//         idx.insert(idx.end(),calc.indexes.begin(),calc.indexes.end());
//         out.insert(out.end(),gamma_instruct2.begin(),gamma_instruct2.end());
//         break;
//       }
//     case Fam::beta:
//       {
//         instructs beta_instruct = {Do::PushVariance,Do::Subtract,Do::PushY,Do::Log,Do::Multiply,
//                                    Do::Int1};
//         instructs beta_instruct2 = {Do::Int1,Do::Subtract,Do::PushVariance,Do::Multiply,Do::Subtract,
//                                     Do::PushY,Do::Int1,Do::Subtract,Do::Log,
//                                     Do::Multiply,Do::Add};
//         instructs beta_instruct3 = {Do::PushVariance,Do::Multiply,Do::Gamma,Do::Log,Do::Negate,Do::Add};
//         instructs beta_instruct4 = {Do::Int1,Do::Subtract,Do::PushVariance,Do::Multiply,Do::Gamma,Do::Log,
//                                     Do::Negate,Do::Add,Do::PushVariance,
//                                     Do::Gamma,Do::Log,Do::Add};
//         out.push_back(Do::Int1);
//         out.insert(out.end(),calc.instructions.begin(),calc.instructions.end());
//         idx.insert(idx.end(),calc.indexes.begin(),calc.indexes.end());
//         out.insert(out.end(),beta_instruct.begin(),beta_instruct.end());
//         out.insert(out.end(),calc.instructions.begin(),calc.instructions.end());
//         idx.insert(idx.end(),calc.indexes.begin(),calc.indexes.end());
//         out.insert(out.end(),beta_instruct2.begin(),beta_instruct2.end());
//         out.insert(out.end(),calc.instructions.begin(),calc.instructions.end());
//         idx.insert(idx.end(),calc.indexes.begin(),calc.indexes.end());
//         out.insert(out.end(),beta_instruct3.begin(),beta_instruct3.end());
//         out.insert(out.end(),calc.instructions.begin(),calc.instructions.end());
//         idx.insert(idx.end(),calc.indexes.begin(),calc.indexes.end());
//         out.insert(out.end(),beta_instruct4.begin(),beta_instruct4.end());
//         break;
//       }
//     case Fam::binomial:
//       {
//         instructs binom_instruct = {Do::PushY,Do::LogFactorialApprox,Do::PushY,Do::PushVariance,
//                                     Do::Subtract,Do::Add,Do::PushVariance,
//                                     Do::LogFactorialApprox,Do::Add};
//         instructs binom_instruct2 = {Do::Log,Do::PushY,Do::Multiply,Do::Add};
//         instructs binom_instruct3 = {Do::Int1,Do::Subtract,Do::Log,Do::PushY,Do::PushVariance,
//                                      Do::Subtract,Do::Multiply,Do::Add};
//         out.insert(out.end(),binom_instruct.begin(),binom_instruct.end());
//         out.insert(out.end(),calc.instructions.begin(),calc.instructions.end());
//         idx.insert(idx.end(),calc.indexes.begin(),calc.indexes.end());
//         out.insert(out.end(),binom_instruct2.begin(),binom_instruct2.end());
//         out.insert(out.end(),calc.instructions.begin(),calc.instructions.end());
//         idx.insert(idx.end(),calc.indexes.begin(),calc.indexes.end());
//         out.insert(out.end(),binom_instruct3.begin(),binom_instruct3.end());
//       }
//   }
//   calc.instructions = out;
//   calc.indexes = idx;
// }

}
