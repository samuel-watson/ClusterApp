#pragma once

#include "general.h"

namespace glmmr {

class Family{
public:
  Fam     family;
  Link    link;
  Family(std::string family_, std::string link_): family(str_to_family.at(family_)), link(str_to_link.at(link_)) {};
  Family(){};
  Family(const glmmr::Family& fam) : family(fam.family), link(fam.link) {};
};
}
