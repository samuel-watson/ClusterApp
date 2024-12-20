#pragma once

#include "general.h"

namespace glmmr {

enum class Fam {
  gaussian = 0,
    bernoulli = 1,
    poisson = 2,
    gamma = 3,
    beta = 4,
    binomial = 5,
    quantile = 6, // quantile is the asymmetric Laplacian distribution
    quantile_scaled = 7
};

enum class Link {
  logit = 0,
    loglink = 1, // to avoid conflicting with log() function
    probit = 2,
    identity = 3,
    inverse = 4
};

const std::map<str, Fam> str_to_family = {
  {"gaussian",Fam::gaussian},
  {"bernoulli",Fam::bernoulli},
  {"poisson",Fam::poisson},
  {"gamma",Fam::gamma},
  {"Gamma",Fam::gamma},
  {"beta",Fam::beta},
  {"binomial",Fam::binomial},
  {"quantile",Fam::quantile},
  {"quantile_scaled",Fam::quantile_scaled}
};

const std::map<str, Link> str_to_link = {
  {"logit",Link::logit},
  {"log",Link::loglink},
  {"probit",Link::probit},
  {"identity",Link::identity},
  {"inverse",Link::inverse}
};

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
