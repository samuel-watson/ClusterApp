#pragma once

#include "general.h"

namespace glmmr {
const static std::unordered_map<std::string,int> model_to_int{
  {"poissonlog",1},
  {"poissonidentity",2},
  {"bernoullilogit",3},
  {"bernoullilog",4},
  {"bernoulliidentity",5},
  {"bernoulliprobit",6},
  {"gaussianidentity",7},
  {"gaussianlog",8},
  {"gammalog",9},
  {"gammainverse",10},
  {"gammaidentity",11},
  {"betalogit",12},
  {"binomiallogit",13},
  {"binomiallog",14},
  {"binomialidentity",15},
  {"binomialprobit",16}
};

class Family{
public:
  std::string family;
  std::string link;
  int flink;
  Family(std::string family_, std::string link_): family(family_), link(link_) {flink = glmmr::model_to_int.at(family_+link_);};
  Family(){};
};
}
