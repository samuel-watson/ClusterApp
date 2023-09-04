#pragma once

#include "general.h"
#include "covariance.hpp"
#include "linearpredictor.hpp"
#include "family.hpp"
#include "modelextradata.hpp"
#include "calculator.hpp"
#include "formula.hpp"

namespace glmmr {

    using namespace Eigen;

    template<typename cov, typename linpred>
    class ModelBits {
    public:
        glmmr::Formula formula;
        cov covariance;
        linpred linear_predictor;
        glmmr::ModelExtraData data;
        glmmr::Family family;
        glmmr::calculator calc;
        glmmr::calculator vcalc;
        bool weighted = false;
        ModelBits(const std::string& formula_,
            const ArrayXXd& data_,
            const strvec& colnames_,
            std::string family_,
            std::string link_) :
            formula(formula_),
            covariance(formula, data_, colnames_),
            linear_predictor(formula, data_, colnames_),
            data(data_.rows()),
            family(family_, link_) {
            setup_calculator();
        };
        int n() { return linear_predictor.n(); };
        ArrayXd xb() { return linear_predictor.xb() + data.offset; };
        virtual void make_covariance_sparse();
        virtual void make_covariance_dense();

    private:
        void setup_calculator();
    };

}

template<typename cov, typename linpred>
inline void glmmr::ModelBits<cov, linpred>::setup_calculator() {
    dblvec yvec(n(), 0.0);
    calc = linear_predictor.calc;
    glmmr::linear_predictor_to_link(calc, family.link);
    glmmr::link_to_likelihood(calc, family.family);
    calc.y = yvec;
    calc.variance.conservativeResize(yvec.size());
    calc.variance = data.variance;
    vcalc = linear_predictor.calc;
    glmmr::re_linear_predictor(vcalc, covariance.Q());
    glmmr::linear_predictor_to_link(vcalc, family.link);
    glmmr::link_to_likelihood(vcalc, family.family);
    vcalc.y = yvec;
    vcalc.variance.conservativeResize(yvec.size());
    vcalc.variance = data.variance;
}

template<typename cov, typename linpred>
void glmmr::ModelBits<cov, linpred>::make_covariance_sparse() {
    covariance.set_sparse(true);
}

template<typename cov, typename linpred>
void glmmr::ModelBits<cov, linpred>::make_covariance_dense() {
    covariance.set_sparse(false);
}
