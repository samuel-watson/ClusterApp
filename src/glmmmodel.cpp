#include "clusterclasses.h"
#include <boost/math/distributions/students_t.hpp>

void ClusterApp::glmmModel::update_formula() {
    std::string new_formula = "int";
    if (option.two_treatments) {
        new_formula += "+int2+int12";
    }
    if (statmodel.beta_pars.size() > 1) {
        if (statmodel.linearpredictor == ClusterApp::LinearPredictor::time_fixed_effects) {
            new_formula += "+factor(t)";
        }
        else {
            new_formula += "+factor(t)";
            //need to update this so that it creates new data with the cluster linear trends; add these by default to the data?
        }
    }
    if (statmodel.include_intercept == 0)new_formula += "-1";
    // covariance
    switch (statmodel.covariance) {
    case ClusterApp::Covariance::exchangeable:
        new_formula += "+(1|gr(cl))";
        break;
    case ClusterApp::Covariance::nested_exchangeable:
        new_formula += "+(1|gr(cl))+(1|gr(cl,t))";
        break;
    case ClusterApp::Covariance::autoregressive:
        new_formula += "+(1|gr(cl)*ar0(t))";
        break;
    case ClusterApp::Covariance::exponential:
        new_formula += "+(1|gr(cl)*fexp0(t))";
        break;
    case ClusterApp::Covariance::squared_exponential:
        new_formula += "+(1|gr(cl)*sqexp0(t))";
        break;
    }
    if (statmodel.sampling == ClusterApp::Sampling::cohort) {
        switch (statmodel.ind_covariance) {
        case ClusterApp::IndividualCovariance::exchangeable:
            new_formula += "+(1|gr(cl))";
            break;
        case ClusterApp::IndividualCovariance::autoregressive:
            new_formula += "+(1|gr(cl)*ar0(t))";
            break;
        case ClusterApp::IndividualCovariance::exponential:
            new_formula += "+(1|gr(cl)*fexp0(t))";
            break;
        case ClusterApp::IndividualCovariance::squared_exponential:
            new_formula += "+(1|gr(cl)*sqexp0(t))";
            break;
        }
    }
    formula = new_formula;
    switch (statmodel.family) {
    case ClusterApp::Family::gaussian:
        family = "gaussian";
        break;
    case ClusterApp::Family::binomial:
        family = "binomial";
        break;
    case ClusterApp::Family::poisson:
        family = "poisson";
        break;
    case ClusterApp::Family::gamma:
        family = "gamma";
        break;
    case ClusterApp::Family::beta:
        family = "beta";
        break;
    }
    switch (statmodel.link) {
    case ClusterApp::Link::identity:
        link = "identity";
        break;
    case ClusterApp::Link::log:
        link = "log";
        break;
    case ClusterApp::Link::logit:
        link = "logit";
        break;
    case ClusterApp::Link::probit:
        link = "probit";
        break;
    case ClusterApp::Link::inverse:
        link = "inverse";
        break;
    }
}
void ClusterApp::glmmModel::update_parameters() {
    std::vector<double> beta;
    double mean_n = designs.mean_n();
    beta.push_back(statmodel.te_pars[0]);
    if (option.two_treatments) {
        beta.push_back(statmodel.te_pars[1]);
        beta.push_back(statmodel.te_pars[2]);
    }
    beta.insert(beta.end(), statmodel.beta_pars.begin(), statmodel.beta_pars.end());
    std::vector<double> theta;
    if (statmodel.family == ClusterApp::Family::gaussian) {
        if (statmodel.covariance == ClusterApp::Covariance::exchangeable) {
            theta.push_back(statmodel.ixx_pars[0]);
        }       
        else if (statmodel.covariance == ClusterApp::Covariance::nested_exchangeable) {
            double tau1 = statmodel.ixx_pars[0] * statmodel.ixx_pars[1];
            double tau2 = statmodel.ixx_pars[0] * (1 - statmodel.ixx_pars[1]);
            theta.push_back(tau1);
            theta.push_back(tau2);
        }
        else if (statmodel.covariance == ClusterApp::Covariance::autoregressive || statmodel.covariance == ClusterApp::Covariance::exponential || statmodel.covariance == ClusterApp::Covariance::squared_exponential) {
            theta.push_back(statmodel.ixx_pars[0]);
            theta.push_back(statmodel.cov_pars[1]);
        }

        if (statmodel.sampling == ClusterApp::Sampling::cohort) {
            double tau3 = statmodel.ixx_pars[2] * (1 - statmodel.ixx_pars[0]);
            tau3 = tau3 / mean_n;
            theta.push_back(tau3);
            if (statmodel.ind_covariance != ClusterApp::IndividualCovariance::exchangeable) {
                theta.push_back(statmodel.cov_pars[4]);
            }
        }

        if (statmodel.link == ClusterApp::Link::identity) {
            if (statmodel.sampling == ClusterApp::Sampling::cohort) {
                double tau4 = (1 - statmodel.ixx_pars[0]) * (1 - statmodel.ixx_pars[2]);
                (*model).model.data.set_var_par(tau4);
            }
            else {
                (*model).model.data.set_var_par(1 - statmodel.ixx_pars[0]);
            }           
        }
        else {
            (*model).model.data.set_var_par(1 - statmodel.cov_pars[2]);
        }

    }
    else {
        theta.push_back(statmodel.cov_pars[0]);
        if (statmodel.covariance != ClusterApp::Covariance::exchangeable) {
            theta.push_back(statmodel.cov_pars[1]);
        }
        if (statmodel.sampling == ClusterApp::Sampling::cohort) {
            double tau3 = statmodel.cov_pars[3];
            tau3 = tau3 / mean_n;
            theta.push_back(tau3);
            if (statmodel.ind_covariance != ClusterApp::IndividualCovariance::exchangeable) {
                theta.push_back(statmodel.cov_pars[4]);
            }
        }
        if (statmodel.family == ClusterApp::Family::beta || statmodel.family == ClusterApp::Family::gamma) {
            (*model).model.data.set_var_par(1 - statmodel.cov_pars[2]);
        }
    }

    (*model).update_beta(beta);
    (*model).update_theta(theta);
    (*model).matrix.W.update();
}

void ClusterApp::glmmModel::update_model_data(const Eigen::ArrayXXd& data) {
    if (modelbits)modelbits.release();
    if (model)model.release();
    modelbits = std::unique_ptr<glmmr::ModelBits>(new glmmr::ModelBits(formula,
        data, colnames, family, link));
    model = std::unique_ptr<glmmr::Model>(new glmmr::Model(*modelbits));
    // update weights
    Eigen::ArrayXd weights(data.rows());
    for (int i = 0; i < data.rows(); i++) {
        if (statmodel.family == ClusterApp::Family::poisson) {
            weights(i) = log(data(i, 2));
        }
        else {
            weights(i) = data(i, 2);
        }
    }
    if (statmodel.family == ClusterApp::Family::gaussian) {
        (*model).set_weights(weights);
    }
    else if (statmodel.family == ClusterApp::Family::binomial) {
        (*model).model.data.set_variance(weights);
    }
    else if (statmodel.family == ClusterApp::Family::poisson) {
        (*model).set_offset(weights);
    }
    update_parameters();
};

void ClusterApp::glmmModel::power(ClusterApp::modelSummary& summary) {
    if (model) {
        Eigen::MatrixXd M = (*model).matrix.information_matrix();
        M = M.llt().solve(Eigen::MatrixXd::Identity(M.rows(), M.cols()));
        int idx = statmodel.include_intercept == 1 ? 1 : 0;
        double zval;
        double bvar = M(idx, idx);
        if (!isnan(bvar) && bvar > 0) {
            zval = abs(statmodel.te_pars[0] / sqrt(bvar));
            summary.power = boost::math::cdf(norm, zval - zcutoff) * 100;
            summary.dof = (*model).model.n();
            summary.se = sqrt(M(idx, idx));
            summary.ci_width = zcutoff * summary.se;
        }
        else {
            summary.power = 909;
            summary.dof = 0;
            summary.se = 0;
            summary.ci_width = 0;
        }
        if (option.two_treatments) {
            bvar = M(idx + 1, idx + 1);
            if (!isnan(bvar) && bvar > 0) {
                zval = abs(statmodel.te_pars[1] / sqrt(bvar));
                summary.power_2 = boost::math::cdf(norm, zval - zcutoff) * 100;
                summary.dof_2 = (*model).model.n();
                summary.se_2 = sqrt(M(idx + 1, idx + 1));
                summary.ci_width_2 = zcutoff * summary.se_2;
            }
            else {
                summary.power_2 = 909;
                summary.dof_2 = 0;
                summary.se_2 = 0;
                summary.ci_width_2 = 0;
            }
            bvar = M(idx + 2, idx + 2);
            if (!isnan(bvar) && bvar > 0) {
                zval = abs(statmodel.te_pars[2] / sqrt(bvar));
                summary.power_12 = boost::math::cdf(norm, zval - zcutoff) * 100;
                summary.dof_12 = (*model).model.n();
                summary.se_12 = sqrt(M(idx + 2, idx + 2));
                summary.ci_width_12 = zcutoff * summary.se_12;
            }
            else {
                summary.power_12 = 909;
                summary.dof_12 = 0;
                summary.se_12 = 0;
                summary.ci_width_12 = 0;
            }
            
        }
    }
}

void ClusterApp::glmmModel::power_kr(ClusterApp::modelSummary& summary) {
    if (model) {
        kenward_data res = (*model).matrix.kenward_roger();
        int idx = statmodel.include_intercept == 1 ? 1 : 0;
        double dofkr = res.dof(idx) > 1 ? res.dof(idx) : 1.0;
        boost::math::students_t dist(dofkr);
        double bvar = res.vcov_beta(idx, idx);
        double tval, tcutoff;
        if (!isnan(bvar) && bvar > 0) {
            tval = abs(statmodel.te_pars[0] / sqrt(bvar));
            tcutoff = boost::math::quantile(dist, 0.975);
            summary.power_kr = res.dof(idx) > 1 ? boost::math::cdf(dist, tval - tcutoff) * 100 : 0.0;
            summary.dof_kr = res.dof(idx);
            summary.se_kr = sqrt(bvar);
            summary.ci_width_kr = tcutoff * summary.se_kr;            
        }
        else {
            summary.power_kr = 909;
            summary.dof_kr = 0;
            summary.se_kr = 0;
            summary.ci_width_kr = 0;
        }

        if (option.two_treatments) {
            bvar = res.vcov_beta(idx + 1, idx + 1);
            if (!isnan(bvar) && bvar >= 0) {
                tval = abs(statmodel.te_pars[1] / sqrt(bvar));
                summary.power_kr_2 = res.dof(idx + 1) > 1 ? boost::math::cdf(dist, tval - tcutoff) * 100 : 0.0;
                summary.dof_kr_2 = res.dof(idx + 1);
                summary.se_kr_2 = sqrt(bvar);
                summary.ci_width_kr_2 = tcutoff * summary.se_2;
            }
            else {
                summary.power_kr_2 = 909;
                summary.dof_kr_2 = 0;
                summary.se_kr_2 = 0;
                summary.ci_width_kr_2 = 0;
            }
            
            bvar = res.vcov_beta(idx + 2, idx + 2);
            if (!isnan(bvar) && bvar >= 0) {
                tval = abs(statmodel.te_pars[2] / sqrt(bvar));
                summary.power_kr_12 = res.dof(idx + 2) > 1 ? boost::math::cdf(dist, tval - tcutoff) * 100 : 0.0;
                summary.dof_kr_12 = res.dof(idx + 2);
                summary.se_kr_12 = sqrt(bvar);
                summary.ci_width_kr_12 = tcutoff * summary.se_12;
            }
            else {
                summary.power_kr_12 = 909;
                summary.dof_kr_12 = 0;
                summary.se_kr_12 = 0;
                summary.ci_width_kr_12 = 0;
            }
            
        }
    }
}

void ClusterApp::glmmModel::power_bw(ClusterApp::modelSummary& summary) {
    if (model) {
        Eigen::MatrixXd M = (*model).matrix.information_matrix();
        M = M.llt().solve(Eigen::MatrixXd::Identity(M.rows(), M.cols()));
        int idx = statmodel.include_intercept == 1 ? 1 : 0;
        double zval;
        double bvar = M(idx, idx);        
        double dofbw = dof - (*model).model.linear_predictor.P();
        dofbw = dofbw > 1 ? dofbw : 1.0;
        boost::math::students_t dist(dofbw);
        double tcutoff = boost::math::quantile(dist, 0.975);
        if (!isnan(bvar) && bvar > 0) {
            zval = abs(statmodel.te_pars[0] / sqrt(bvar));                         
            summary.power_bw = dofbw == 1 ? 0.0 : boost::math::cdf(dist, zval - tcutoff) * 100;
            summary.dof_bw = dof - (*model).model.linear_predictor.P();
            summary.ci_width_bw = tcutoff * summary.se;
        }
        else {
            summary.power_bw = 0;
            summary.dof_bw = 0;
            summary.ci_width_bw = 0;
        }
        if (option.two_treatments) {
            bvar = M(idx + 1, idx + 1);
            if (!isnan(bvar) && bvar > 0) {
                zval = abs(statmodel.te_pars[1] / sqrt(bvar));
                summary.power_bw_2 = boost::math::cdf(dist, zval - tcutoff) * 100;
                summary.dof_bw_2 = dofbw;
                summary.ci_width_bw_2 = tcutoff * summary.se_2;
            }
            else {
                summary.power_bw_2 = 909;
                summary.dof_bw_2 = 0;
                summary.ci_width_bw_2 = 0;
            }            
            bvar = M(idx + 2, idx + 2);
            if (!isnan(bvar) && bvar > 0) {
                zval = abs(statmodel.te_pars[2] / sqrt(M(idx + 2, idx + 2)));
                summary.power_bw_12 = boost::math::cdf(dist, zval - tcutoff) * 100;
                summary.dof_bw_12 = dofbw;
                summary.ci_width_bw_12 = tcutoff * summary.se_12;
            }
            else {
                summary.power_bw_12 = 909;
                summary.dof_bw_12 = 0;
                summary.ci_width_bw_12 = 0;
            }
            
        }
    }
}

void ClusterApp::glmmModel::optimum(int N) {
    if (model) {
        Eigen::VectorXd C = Eigen::VectorXd::Zero((*model).model.linear_predictor.P());
        int idx = statmodel.include_intercept == 1 ? 1 : 0;
        C(idx) = statmodel.c_vals[0];
        if (option.two_treatments) {
            C(idx+1) = statmodel.c_vals[1];
            C(idx + 2) = statmodel.c_vals[2];
        }
        Eigen::ArrayXd weights = (*model).optim.optimum_weights(N,C);
        if (weights.size() != optimal_weights.size()) optimal_weights.resize(weights.size());
        for (int i = 0; i < weights.size(); i++) optimal_weights[i] = weights(i);
    }
}
