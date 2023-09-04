#pragma once

#include "general.h"
#include "modelbits.hpp"
#include "randomeffects.hpp"
#include "modelmatrix.hpp"
#include "modelmcmc.hpp"
#include "modeloptim.hpp"

namespace glmmr {

    using namespace Eigen;

    template<typename modeltype>
    class Model {
    public:
        modeltype model;
        glmmr::RandomEffects<modeltype> re;
        glmmr::ModelMatrix<modeltype> matrix;
        glmmr::ModelOptim<modeltype> optim;
        glmmr::ModelMCMC<modeltype> mcmc;

        Model(const std::string& formula_,
            const ArrayXXd& data_,
            const strvec& colnames_,
            std::string family_,
            std::string link_) : model(formula_, data_, colnames_, family_, link_), re(model), matrix(model, re), optim(model, matrix, re), mcmc(model, matrix, re) {};

        virtual void set_offset(const VectorXd& offset_);
        virtual void set_weights(const ArrayXd& weights_);
        virtual void set_y(const VectorXd& y_);
        virtual void update_beta(const dblvec& beta_);
        virtual void update_theta(const dblvec& theta_);
        virtual void update_u(const MatrixXd& u_);
        virtual void set_trace(int trace_);
        // add in functions to run the whole fitting algorithm here
    };

}

template<typename modeltype>
inline void glmmr::Model<modeltype>::set_offset(const VectorXd& offset_) {
    model.data.set_offset(offset_);
}

template<typename modeltype>
inline void glmmr::Model<modeltype>::set_weights(const ArrayXd& weights_) {
    model.data.set_weights(weights_);
    if ((weights_ != 1.0).any()) {
        model.weighted = true;
    }
}

template<typename modeltype>
inline void glmmr::Model<modeltype>::set_y(const VectorXd& y_) {
    model.data.update_y(y_);
}

template<typename modeltype>
inline void glmmr::Model<modeltype>::update_beta(const dblvec& beta_) {
    model.linear_predictor.update_parameters(beta_);
}

template<typename modeltype>
inline void glmmr::Model<modeltype>::update_theta(const dblvec& theta_) {
    model.covariance.update_parameters(theta_);
    re.ZL = model.covariance.ZL_sparse();
    re.zu_ = re.ZL * re.u_;
}

template<typename modeltype>
inline void glmmr::Model<modeltype>::update_u(const MatrixXd& u_) {
    if (u_.cols() != re.u(false).cols()) {
        re.u_.conservativeResize(model.covariance.Q(), u_.cols());
        re.zu_.conservativeResize(model.covariance.Q(), u_.cols());
    }
    re.u_ = u_;
    re.zu_ = re.ZL * re.u_;
}

template<typename modeltype>
inline void glmmr::Model<modeltype>::set_trace(int trace_) {
    optim.trace = trace_;
    mcmc.trace = trace_;
    if (trace_ > 0) {
        mcmc.verbose = true;
    }
    else {
        mcmc.verbose = false;
    }
}
