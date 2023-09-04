#pragma once

#include "general.h"
#include "modelbits.hpp"
#include "randomeffects.hpp"
#include "modelmatrix.hpp"
#include "openmpheader.h"
#include "maths.h"
#include "algo.h"
#include "sparseext.h"
#include "calculator.hpp"

namespace glmmr {

    using namespace rminqa;
    using namespace Eigen;

    template<typename modeltype>
    class ModelOptim {
    public:
        modeltype& model;
        glmmr::ModelMatrix<modeltype>& matrix;
        glmmr::RandomEffects<modeltype>& re;
        int trace = 0;

        ModelOptim(modeltype& model_,
            glmmr::ModelMatrix<modeltype>& matrix_,
            glmmr::RandomEffects<modeltype>& re_) : model(model_), matrix(matrix_), re(re_) {};

        virtual void update_beta(const dblvec& beta);
        virtual void update_theta(const dblvec& theta);
        virtual void update_beta(const VectorXd& beta);
        virtual void update_theta(const VectorXd& theta);
        virtual void update_u(const MatrixXd& u_);
        virtual double log_likelihood();
        virtual double full_log_likelihood();
        virtual void nr_beta();
        virtual void laplace_nr_beta_u();
        virtual void update_var_par(const double& v);
        virtual void update_var_par(const ArrayXd& v);
        virtual void ml_beta();
        virtual void ml_theta();
        virtual void ml_all();
        virtual void laplace_ml_beta_u();
        virtual void laplace_ml_theta();
        virtual void laplace_ml_beta_theta();
        virtual double aic();
        virtual ArrayXd optimum_weights(double N, VectorXd C, double tol = 1e-5, int max_iter = 501);

    protected:
        void calculate_var_par();
        dblvec get_start_values(bool beta, bool theta, bool var = true);
        dblvec get_lower_values(bool beta, bool theta, bool var = true);
        dblvec get_upper_values(bool beta, bool theta, bool var = true);

        class L_likelihood : public Functor<dblvec> {
            ModelOptim<modeltype>& M;
            double ll;
        public:
            L_likelihood(ModelOptim<modeltype>& M_) :
                M(M_), ll(0.0) {};
            double operator()(const dblvec& par);
        };

        class D_likelihood : public Functor<dblvec> {
            ModelOptim<modeltype>& M;
            const MatrixXd& Lu;
            double logl;
        public:
            D_likelihood(ModelOptim<modeltype>& M_,
                const MatrixXd& Lu_) :
                M(M_),
                Lu(Lu_),
                logl(0.0) {};
            double operator()(const dblvec& par);
        };

        class F_likelihood : public Functor<dblvec> {
            ModelOptim<modeltype>& M;
            int G;
            bool importance;
            double ll;
            double denomD;
        public:
            F_likelihood(ModelOptim<modeltype>& M_,
                double denomD_ = 0,
                bool importance_ = false) :
                M(M_),
                G(M_.model.covariance.npar()),
                importance(importance_),
                ll(0.0),
                denomD(denomD_) {};
            double operator()(const dblvec& par);
        };

        class LA_likelihood : public Functor<dblvec> {
            ModelOptim<modeltype>& M;
            MatrixXd v;
            MatrixXd LZWZL;
            double LZWdet;
            double logl;
            double ll;
        public:
            LA_likelihood(ModelOptim<modeltype>& M_) :
                M(M_),
                v(M.model.covariance.Q(), 1),
                LZWZL(MatrixXd::Zero(M.model.covariance.Q(), M.model.covariance.Q())),
                LZWdet(0.0),
                logl(0.0), ll(0.0) {
                M.matrix.W.update();
                LZWZL = M.model.covariance.LZWZL(M.matrix.W.W());
                LZWdet = glmmr::maths::logdet(LZWZL);
            };
            double operator()(const dblvec& par);
        };

        class LA_likelihood_cov : public Functor<dblvec> {
            ModelOptim<modeltype>& M;
            MatrixXd LZWZL;
            double LZWdet;
            double logl;
            double ll;
        public:
            LA_likelihood_cov(ModelOptim<modeltype>& M_) :
                M(M_),
                LZWZL(MatrixXd::Zero(M.model.covariance.Q(), M.model.covariance.Q())),
                LZWdet(0.0), logl(0.0), ll(0.0) {};
            double operator()(const dblvec& par);
        };

        class LA_likelihood_btheta : public Functor<dblvec> {
            ModelOptim<modeltype>& M;
            MatrixXd LZWZL;
            double LZWdet;
            double logl;
            double ll;
        public:
            LA_likelihood_btheta(ModelOptim<modeltype>& M_) :
                M(M_),
                LZWZL(MatrixXd::Zero(M.model.covariance.Q(), M.model.covariance.Q())),
                LZWdet(0.0), logl(0.0), ll(0.0) {};
            double operator()(const dblvec& par);
        };

    };

}

template<typename modeltype>
inline void glmmr::ModelOptim<modeltype>::update_beta(const dblvec& beta) {
    model.linear_predictor.update_parameters(beta);
}

template<typename modeltype>
inline void glmmr::ModelOptim<modeltype>::update_beta(const VectorXd& beta) {
    model.linear_predictor.update_parameters(beta.array());
}

template<typename modeltype>
inline void glmmr::ModelOptim<modeltype>::update_theta(const dblvec& theta) {
    model.covariance.update_parameters(theta);
    re.ZL = model.covariance.ZL_sparse();
    re.zu_ = re.ZL * re.u_;
}

template<typename modeltype>
inline void glmmr::ModelOptim<modeltype>::update_theta(const VectorXd& theta) {
    model.covariance.update_parameters(theta.array());
    re.ZL = model.covariance.ZL_sparse();
    re.zu_ = re.ZL * re.u_;
}

template<typename modeltype>
inline void glmmr::ModelOptim<modeltype>::update_u(const MatrixXd& u_) {
    if (u_.cols() != re.u(false).cols()) {
        re.u_.conservativeResize(model.covariance.Q(), u_.cols());
        re.zu_.resize(model.covariance.Q(), u_.cols());
    }
    re.u_ = u_;
    re.zu_ = re.ZL * re.u_;
}

template<typename modeltype>
inline double glmmr::ModelOptim<modeltype>::log_likelihood() {
    double ll = 0;
    ArrayXd xb(model.xb());
    if (model.weighted) {
        if (model.family.family == "gaussian") {
#pragma omp parallel for reduction (+:ll) collapse(2)
            for (int j = 0; j < re.Zu().cols(); j++) {
                for (int i = 0; i < model.n(); i++) {
                    ll += glmmr::maths::log_likelihood(model.data.y(i), xb(i) + re.zu_(i, j), model.data.variance(i) / model.data.weights(i), model.family.flink);
                }
            }
        }
        else {
#pragma omp parallel for reduction (+:ll) collapse(2)
            for (int j = 0; j < re.Zu().cols(); j++) {
                for (int i = 0; i < model.n(); i++) {
                    ll += model.data.weights(i) * glmmr::maths::log_likelihood(model.data.y(i), xb(i) + re.zu_(i, j), model.data.variance(i), model.family.flink);
                }
            }
            ll *= model.data.weights.sum() / model.n();
        }
    }
    else {
#pragma omp parallel for reduction (+:ll) collapse(2)
        for (int j = 0; j < re.Zu().cols(); j++) {
            for (int i = 0; i < model.n(); i++) {
                ll += glmmr::maths::log_likelihood(model.data.y(i), xb(i) + re.zu_(i, j), model.data.variance(i), model.family.flink);
            }
        }
    }

    // to use the calculator object instead... seems to be generally slower so have opted 
    // for specific formulae above. Will try to optimise this in future versions
    // #pragma omp parallel for reduction (+:ll) collapse(2) 
    //  for(int j=0; j<zu_.cols() ; j++){
    //    for(int i = 0; i<n_; i++){
    //      double ozu = offset_(i)+zu_(i,j);
    //      ll += calc_.calculate(i,linpred_.parameters_,linpred_.Xdata_,0,0,ozu)[0];
    //    }
    //  }

    return ll / re.Zu().cols();
}

template<typename modeltype>
inline double glmmr::ModelOptim<modeltype>::full_log_likelihood() {
    double ll = log_likelihood();
    double logl = 0;
    MatrixXd Lu = model.covariance.Lu(re.u(false));
#pragma omp parallel for reduction (+:logl)
    for (int i = 0; i < Lu.cols(); i++) {
        logl += model.covariance.log_likelihood(Lu.col(i));
    }
    logl *= 1 / Lu.cols();
    return ll + logl;
}

template<typename modeltype>
inline dblvec glmmr::ModelOptim<modeltype>::get_start_values(bool beta, bool theta, bool var) {
    dblvec start;
    if (beta) {
        for (int i = 0; i < model.linear_predictor.P(); i++)start.push_back(model.linear_predictor.parameters[i]);
        if (theta) {
            for (int i = 0; i < model.covariance.npar(); i++) {
                start.push_back(model.covariance.parameters_[i]);
            }
        }
    }
    else {
        start = model.covariance.parameters_;
    }
    if (var && (model.family.family == "gaussian" || model.family.family == "Gamma" || model.family.family == "beta")) {
        start.push_back(model.data.var_par);
    }
    return start;
}

template<typename modeltype>
inline dblvec glmmr::ModelOptim<modeltype>::get_lower_values(bool beta, bool theta, bool var) {
    dblvec lower;
    if (beta) {
        for (int i = 0; i < model.linear_predictor.P(); i++) {
            lower.push_back(-100000000);
        }
    }
    if (theta) {
        for (int i = 0; i < model.covariance.npar(); i++) {
            lower.push_back(1e-6);
        }
    }
    if (var && (model.family.family == "gaussian" || model.family.family == "Gamma" || model.family.family == "beta")) {
        lower.push_back(0.0);
    }
    return lower;
}

template<typename modeltype>
inline dblvec glmmr::ModelOptim<modeltype>::get_upper_values(bool beta, bool theta, bool var) {
    dblvec upper;
    if (beta) {
        for (int i = 0; i < model.linear_predictor.P(); i++) {
            upper.push_back(100000000);
        }
    }
    if (theta) {
        for (int i = 0; i < model.covariance.npar(); i++) {
            upper.push_back(100000000);
        }
    }
    if (var && (model.family.family == "gaussian" || model.family.family == "Gamma" || model.family.family == "beta")) {
        upper.push_back(100000000);
    }
    return upper;
}

template<typename modeltype>
inline void glmmr::ModelOptim<modeltype>::nr_beta() {
    int niter = re.u(false).cols();
    MatrixXd zd = matrix.linpred();
    ArrayXd sigmas(niter);
    if (model.linear_predictor.any_nonlinear()) {
        vector_matrix score = matrix.b_score();
        MatrixXd infomat = score.mat.llt().solve(MatrixXd::Identity(model.linear_predictor.P(), model.linear_predictor.P()));
        VectorXd bplus = infomat * score.vec;
        for (int i = 0; i < bplus.size(); i++)model.linear_predictor.parameters[i] += bplus(i);
    }
    else {
        MatrixXd XtXW = MatrixXd::Zero(model.linear_predictor.P() * niter, model.linear_predictor.P());
        MatrixXd Wu = MatrixXd::Zero(model.n(), niter);
        ArrayXd nvar_par(model.n());
        if (model.family.family == "gaussian") {
            nvar_par = model.data.variance;
        }
        else if (model.family.family == "Gamma") {
            nvar_par = model.data.variance.inverse();
        }
        else if (model.family.family == "beta") {
            nvar_par = (1 + model.data.variance);
        }
        else if (model.family.family == "binomial") {
            nvar_par = model.data.variance.inverse();
        }
        else {
            nvar_par.setConstant(1.0);
        }
#pragma omp parallel for
        for (int i = 0; i < niter; ++i) {
            VectorXd w = glmmr::maths::dhdmu(zd.col(i), model.family);
            w = ((w.array() * nvar_par).inverse() * model.data.weights).matrix();
            VectorXd zdu = glmmr::maths::mod_inv_func(zd.col(i), model.family.link);
            VectorXd dmu = glmmr::maths::detadmu(zd.col(i), model.family.link);
            if (model.family.family == "binomial") {
                zdu = zdu.cwiseProduct(model.data.variance.matrix());
                dmu = dmu.cwiseProduct(model.data.variance.inverse().matrix());
            }
            ArrayXd resid = (model.data.y - zdu);
            XtXW.block(model.linear_predictor.P() * i, 0, model.linear_predictor.P(), model.linear_predictor.P()) = model.linear_predictor.X().transpose() * w.asDiagonal() * model.linear_predictor.X();
            w = w.cwiseProduct(dmu);
            w = w.cwiseProduct(resid.matrix());
            Wu.col(i) = w;
        }
        XtXW *= (double)1 / niter;
        MatrixXd XtWXm = XtXW.block(0, 0, model.linear_predictor.P(), model.linear_predictor.P());
        for (int i = 1; i < niter; i++) XtWXm += XtXW.block(model.linear_predictor.P() * i, 0, model.linear_predictor.P(), model.linear_predictor.P());
        XtWXm = XtWXm.inverse();
        VectorXd Wum = Wu.rowwise().mean();
        VectorXd bincr = XtWXm * (model.linear_predictor.X().transpose()) * Wum;
        update_beta(model.linear_predictor.parameter_vector() + bincr);
    }
    calculate_var_par();
}

template<typename modeltype>
inline void glmmr::ModelOptim<modeltype>::laplace_nr_beta_u() {
    matrix.W.update();
    VectorXd zd = (matrix.linpred()).col(0);
    VectorXd dmu = glmmr::maths::detadmu(zd, model.family.link);
    MatrixXd infomat = matrix.observed_information_matrix();
    infomat = infomat.llt().solve(MatrixXd::Identity(model.linear_predictor.P() + model.covariance.Q(), model.linear_predictor.P() + model.covariance.Q()));
    VectorXd zdu = glmmr::maths::mod_inv_func(zd, model.family.link);
    if (model.family.family == "binomial") {
        zdu = zdu.cwiseProduct(model.data.variance.matrix());
        dmu = dmu.cwiseProduct(model.data.variance.inverse().matrix());
    }
    ArrayXd resid = (model.data.y - zdu).array();
    VectorXd w = matrix.W.W();
    w = w.cwiseProduct(dmu);
    w = w.cwiseProduct(resid.matrix());
    VectorXd params(model.linear_predictor.P() + model.covariance.Q());
    params.head(model.linear_predictor.P()) = model.linear_predictor.parameter_vector();
    params.tail(model.covariance.Q()) = re.u_.col(0);
    VectorXd pderiv(model.linear_predictor.P() + model.covariance.Q());
    pderiv.head(model.linear_predictor.P()) = (model.linear_predictor.X()).transpose() * w;
    pderiv.tail(model.covariance.Q()) = matrix.log_gradient(re.u_.col(0));
    params += infomat * pderiv;
    update_beta(params.head(model.linear_predictor.P()));
    update_u(params.tail(model.covariance.Q()));
    calculate_var_par();
}

template<typename modeltype>
inline void glmmr::ModelOptim<modeltype>::update_var_par(const double& v) {
    model.data.var_par = v;
    model.data.variance.setConstant(v);
    model.calc.variance = model.data.variance;
}

template<typename modeltype>
inline void glmmr::ModelOptim<modeltype>::update_var_par(const ArrayXd& v) {
    model.data.variance = v;
    model.calc.variance = model.data.variance;
}

template<typename modeltype>
inline void glmmr::ModelOptim<modeltype>::calculate_var_par() {
    if (model.family.family == "gaussian") {
        // revise this for beta and Gamma re residuals
        int niter = re.u(false).cols();
        ArrayXd sigmas(niter);
        MatrixXd zd = matrix.linpred();
#pragma omp parallel for
        for (int i = 0; i < niter; ++i) {
            VectorXd zdu = glmmr::maths::mod_inv_func(zd.col(i), model.family.link);
            ArrayXd resid = (model.data.y - zdu);
            resid *= model.data.weights.sqrt();
            sigmas(i) = (resid - resid.mean()).square().sum() / (resid.size() - 1);
        }
        update_var_par(sigmas.mean());
    }
}

template<typename modeltype>
inline void glmmr::ModelOptim<modeltype>::ml_beta() {
    L_likelihood ldl(*this);
    Rbobyqa<L_likelihood, dblvec> opt;
    opt.control.iprint = trace;
    dblvec start = get_start_values(true, false, false);
    dblvec lower = get_lower_values(true, false, false);
    opt.set_lower(lower);
    opt.control.iprint = trace;
    opt.minimize(ldl, start);
    calculate_var_par();
}

template<typename modeltype>
inline void glmmr::ModelOptim<modeltype>::ml_theta() {
    MatrixXd Lu = model.covariance.Lu(re.u(false));
    D_likelihood ddl(*this, Lu);
    Rbobyqa<D_likelihood, dblvec> opt;
    dblvec lower = get_lower_values(false, true, false);
    opt.set_lower(lower);
    opt.control.iprint = trace;
    dblvec start_t = get_start_values(false, true, false);
    opt.minimize(ddl, start_t);
}

template<typename modeltype>
inline void glmmr::ModelOptim<modeltype>::ml_all() {
    MatrixXd Lu = model.covariance.Lu(re.u(false));
    double denomD = 0;
    for (int i = 0; i < Lu.cols(); i++) {
        denomD += model.covariance.log_likelihood(Lu.col(i));
    }
    denomD *= 1 / Lu.cols();
    F_likelihood dl(*this, denomD, true);
    Rbobyqa<F_likelihood, dblvec> opt;
    dblvec start = get_start_values(true, true, false);
    dblvec lower = get_lower_values(true, true, false);
    opt.set_lower(lower);
    opt.control.iprint = trace;
    opt.minimize(dl, start);
    calculate_var_par();
}

template<typename modeltype>
inline void glmmr::ModelOptim<modeltype>::laplace_ml_beta_u() {
    LA_likelihood ldl(*this);
    Rbobyqa<LA_likelihood, dblvec> opt;
    opt.control.iprint = trace;
    dblvec start = get_start_values(true, false, false);
    for (int i = 0; i < model.covariance.Q(); i++)start.push_back(re.u_(i, 0));
    opt.control.iprint = trace;
    opt.minimize(ldl, start);
    calculate_var_par();
}

template<typename modeltype>
inline void glmmr::ModelOptim<modeltype>::laplace_ml_theta() {
    LA_likelihood_cov ldl(*this);
    Rbobyqa<LA_likelihood_cov, dblvec> opt;
    dblvec lower = get_lower_values(false, true, false);
    dblvec start = get_start_values(false, true, false);
    opt.control.iprint = trace;
    opt.set_lower(lower);
    opt.minimize(ldl, start);
}

template<typename modeltype>
inline void glmmr::ModelOptim<modeltype>::laplace_ml_beta_theta() {
    LA_likelihood_btheta ldl(*this);
    Rbobyqa<LA_likelihood_btheta, dblvec> opt;
    dblvec lower = get_lower_values(true, true, false);
    dblvec start = get_start_values(true, true, false);
    opt.set_lower(lower);
    opt.control.iprint = trace;
    opt.minimize(ldl, start);
    calculate_var_par();
}

template<typename modeltype>
inline double glmmr::ModelOptim<modeltype>::L_likelihood::operator()(const dblvec& par) {
    M.update_beta(par);
    ll = M.log_likelihood();
    return -1 * ll;
}

template<typename modeltype>
inline double glmmr::ModelOptim<modeltype>::D_likelihood::operator()(const dblvec& par) {
    M.update_theta(par);
    logl = 0;
#pragma omp parallel for reduction (+:logl)
    for (int i = 0; i < M.re.u(false).cols(); i++) {
        logl += M.model.covariance.log_likelihood(Lu.col(i));
    }
    return -1 * logl / Lu.cols();
}

template<typename modeltype>
inline double glmmr::ModelOptim<modeltype>::F_likelihood::operator()(const dblvec& par) {
    auto first = par.begin();
    auto last1 = par.begin() + M.model.linear_predictor.P();
    auto last2 = par.begin() + M.model.linear_predictor.P() + G;
    dblvec beta(first, last1);
    dblvec theta(last1, last2);
    M.update_beta(beta);
    M.update_theta(theta);
    if (M.model.family.family == "gaussian" || M.model.family.family == "Gamma" || M.model.family.family == "beta")M.update_var_par(par[M.model.linear_predictor.P() + G]);
    ll = M.full_log_likelihood();
    if (importance) {
        return -1.0 * log(exp(ll) / exp(denomD));
    }
    else {
        return -1.0 * ll;
    }
}

template<typename modeltype>
inline double glmmr::ModelOptim<modeltype>::LA_likelihood::operator()(const dblvec& par) {
    logl = 0;
    auto start = par.begin();
    auto end = par.begin() + M.model.linear_predictor.P();
    dblvec beta(start, end);
    for (int i = 0; i < M.model.covariance.Q(); i++)v(i, 0) = par[M.model.linear_predictor.P() + i];
    M.update_beta(beta);
    M.update_u(v);
    logl = v.col(0).transpose() * v.col(0);
    ll = M.log_likelihood();
    if (M.model.family.family != "gaussian") {
        M.matrix.W.update();
        LZWZL = M.model.covariance.LZWZL(M.matrix.W.W());
        LZWdet = glmmr::maths::logdet(LZWZL);
    }
    return -1.0 * (ll - 0.5 * logl - 0.5 * LZWdet);
}

template<typename modeltype>
inline double glmmr::ModelOptim<modeltype>::LA_likelihood_cov::operator()(const dblvec& par) {
    M.update_theta(par);
    M.matrix.W.update();
    logl = M.re.u_.col(0).transpose() * M.re.u_.col(0);
    ll = M.log_likelihood();
    LZWZL = M.model.covariance.LZWZL(M.matrix.W.W());
    LZWdet = glmmr::maths::logdet(LZWZL);
    return -1 * (ll - 0.5 * logl - 0.5 * LZWdet);
}

template<typename modeltype>
inline double glmmr::ModelOptim<modeltype>::LA_likelihood_btheta::operator()(const dblvec& par) {
    auto start = par.begin();
    auto end1 = par.begin() + M.model.linear_predictor.P();
    auto end2 = par.begin() + M.model.linear_predictor.P() + M.model.covariance.npar();
    dblvec beta(start, end1);
    dblvec theta(end1, end2);
    M.update_beta(beta);
    M.update_theta(theta);
    ll = M.log_likelihood();
    logl = M.re.u_.col(0).transpose() * M.re.u_.col(0);
    M.matrix.W.update();
    LZWZL = M.model.covariance.LZWZL(M.matrix.W.W());
    LZWdet = glmmr::maths::logdet(LZWZL);
    return -1 * (ll - 0.5 * logl - 0.5 * LZWdet);
}

template<typename modeltype>
inline double glmmr::ModelOptim<modeltype>::aic() {
    MatrixXd Lu = model.covariance.Lu(re.u(false));
    int dof = model.linear_predictor.P() + model.covariance.npar();
    double logl = 0;
#pragma omp parallel for reduction (+:logl)
    for (int i = 0; i < Lu.cols(); i++) {
        logl += model.covariance.log_likelihood(Lu.col(i));
    }
    double ll = log_likelihood();

    return (-2 * (ll + logl) + 2 * dof);
}

template<typename modeltype>
inline ArrayXd glmmr::ModelOptim<modeltype>::optimum_weights(double N,
    VectorXd C,
    double tol,
    int max_iter) {
    VectorXd Cvec(C);
    ArrayXd weights = ArrayXd::Constant(model.n(), 1.0 * model.n());
    VectorXd holder(model.n());
    weights = weights.inverse();
    ArrayXd weightsnew(weights);
    ArrayXd w = (matrix.W.W()).array().inverse();
    std::vector<MatrixXd> ZDZ;
    std::vector<MatrixXd> Sigmas;
    std::vector<MatrixXd> Xs;
    std::vector<glmmr::SigmaBlock> SB(matrix.get_sigma_blocks());
    int maxprint = model.n() < 10 ? model.n() : 10;
    for (unsigned int i = 0; i < SB.size(); i++) {
        sparse ZLs = submat_sparse(model.covariance.ZL_sparse(), SB[i].RowIndexes);
        MatrixXd ZL = sparse_to_dense(ZLs, false);
        MatrixXd S = ZL * ZL.transpose();
        ZDZ.push_back(S);
        Sigmas.push_back(S);
        ArrayXi rows = Map<ArrayXi, Unaligned>(SB[i].RowIndexes.data(), SB[i].RowIndexes.size());
        MatrixXd X = glmmr::Eigen_ext::submat(model.linear_predictor.X(), rows, ArrayXi::LinSpaced(model.linear_predictor.P(), 0, model.linear_predictor.P() - 1));
        Xs.push_back(X);
    }

    double diff = 1;
    int block_size;
    MatrixXd M(model.linear_predictor.P(), model.linear_predictor.P());
    int iter = 0;
    while (diff > tol && iter < max_iter) {
        iter++;
        //add check to remove weights that are below a certain threshold
        if ((weights < 1e-8).any()) {
            for (unsigned int i = 0; i < SB.size(); i++) {
                auto it = SB[i].RowIndexes.begin();
                while (it != SB[i].RowIndexes.end()) {
                    if (weights(*it) < 1e-8) {
                        weights(*it) = 0;
                        int idx = it - SB[i].RowIndexes.begin();
                        glmmr::Eigen_ext::removeRow(Xs[i], idx);
                        glmmr::Eigen_ext::removeRow(ZDZ[i], idx);
                        glmmr::Eigen_ext::removeColumn(ZDZ[i], idx);
                        Sigmas[i].conservativeResize(ZDZ[i].rows(), ZDZ[i].cols());
                        it = SB[i].RowIndexes.erase(it);
                    }
                    else {
                        it++;
                    }
                }
            }
        }

        M.setZero();
        for (unsigned int i = 0; i < SB.size(); i++) {
            Sigmas[i] = ZDZ[i];
            for (int j = 0; j < Sigmas[i].rows(); j++) {
                // sigma_sq
                Sigmas[i](j, j) += w(SB[i].RowIndexes[j]) / (N * weights(SB[i].RowIndexes[j]));
            }
            Sigmas[i] = Sigmas[i].llt().solve(MatrixXd::Identity(Sigmas[i].rows(), Sigmas[i].cols()));
            M += Xs[i].transpose() * Sigmas[i] * Xs[i];
        }

        //check if positive definite, if not remove the offending column(s)
        bool isspd = glmmr::Eigen_ext::issympd(M);
        if (isspd) {
            ArrayXd M_row_sums = M.rowwise().sum();
            int fake_it = 0;
            int countZero = 0;
            for (int j = 0; j < M_row_sums.size(); j++) {
                if (M_row_sums(j) == 0) {
                    for (unsigned int k = 0; k < Xs.size(); k++) {
                        glmmr::Eigen_ext::removeColumn(Xs[k], fake_it);
                    }
                    glmmr::Eigen_ext::removeElement(Cvec, fake_it);
                    countZero++;
                }
                else {
                    fake_it++;
                }
            }
            M.conservativeResize(M.rows() - countZero, M.cols() - countZero);
            M.setZero();
            for (unsigned int k = 0; k < SB.size(); k++) {
                M += Xs[k].transpose() * Sigmas[k] * Xs[k];
            }
        }
        M = M.llt().solve(MatrixXd::Identity(M.rows(), M.cols()));
        VectorXd Mc = M * Cvec;
        weightsnew.setZero();
        for (unsigned int i = 0; i < SB.size(); i++) {
            block_size = SB[i].RowIndexes.size();
            holder.segment(0, block_size) = Sigmas[i] * Xs[i] * Mc;
            for (int j = 0; j < block_size; j++) {
                weightsnew(SB[i].RowIndexes[j]) = holder(j);
            }
        }
        weightsnew = weightsnew.abs();
        weightsnew *= 1 / weightsnew.sum();
        diff = ((weights - weightsnew).abs()).maxCoeff();
        weights = weightsnew;
    }
    return weights;
}
