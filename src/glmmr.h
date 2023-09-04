#pragma once

#include "glmmr/general.h"
#include "glmmr/maths.h"
#include "glmmr/formula.hpp"
#include "glmmr/covariance.hpp"
#include "glmmr/linearpredictor.hpp"
#include "glmmr/model.hpp"
#include "glmmr/modelbits.hpp"
#include "glmmr/openmpheader.h"

typedef glmmr::ModelBits<glmmr::Covariance, glmmr::LinearPredictor> bits;
typedef glmmr::Model<glmmr::ModelBits<glmmr::Covariance, glmmr::LinearPredictor> > glmm;
