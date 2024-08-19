#pragma once

#include "general.h"

namespace glmmr {

class Family{
public:
  Fam     family;
  Link    link;
  double  quantile = 0.5;
  Family(std::string family_, std::string link_): family(str_to_family.at(family_)), link(str_to_link.at(link_)) {};
  Family(){};
  Family(const glmmr::Family& fam) : family(fam.family), link(fam.link) {};
  
  void set_quantile(const double& q){
    if(q <= 0 || q >= 1) throw std::runtime_error("q !in [0,1]");
    if(! (family == Fam::quantile || family == Fam::quantile_scaled)) throw std::runtime_error("Quantile only relevant for quantile family");
    quantile = q;
  }
};
}
