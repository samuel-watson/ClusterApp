#include "clusterclasses.h"
#include "glmmr/maths.h"


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
    else if (statmodel.sampling == ClusterApp::Sampling::open_cohort) {
        new_formula += "+(1|gr(cl)*ar0(t))";
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

double ClusterApp::glmmModel::mean_individual_variance(bool weighted) {
    Eigen::VectorXd w = glmmr::maths::dhdmu((*model).model.xb(),(*model).model.family);
    // weighted mean
    double var = 0;
    if (weighted) {
        int count = 0;
        for (int i = 0; i < designs.sequences; i++) {
            for (int t = 0; t < designs.time; t++) {
                if (*designs.active(i, t)) {
                    var += w(count) * (*designs.n(i,t));
                }
            }
        }
        var *= (1.0 / designs.total_n());
    }
    else {
        var = w.array().mean();
    }
    if(option.log)logger.AddLog("[%05d] [%s] Calculate mean individual variance %.3f \n", ImGui::GetFrameCount(), logger.cat[0], var);
    return var;
}

std::pair<double, double> ClusterApp::glmmModel::mean_outcome() {   
    Eigen::Vector2d xb;
    xb(0) = statmodel.beta_pars[0];
    xb(1) = statmodel.beta_pars[0] + statmodel.te_pars[0];
    Eigen::Vector2d ymean = glmmr::maths::mod_inv_func(xb, (*model).model.family.link);
    std::pair<double, double> result = { ymean(0), ymean(1) };
    if (option.log)logger.AddLog("[%05d] [%s] Calculate group means [%.3f, %.3f] \n", ImGui::GetFrameCount(), logger.cat[0], result.first, result.second);
    return result;
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

    // if non-gaussian then convert parameters both ways

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
        } else if (statmodel.sampling == ClusterApp::Sampling::open_cohort) {
            double tau3 = statmodel.ixx_pars[2] * (1 - statmodel.ixx_pars[0]);
            tau3 = tau3 / mean_n;
            theta.push_back(tau3);
            theta.push_back(1-statmodel.cov_pars[4]);
        }

        if (statmodel.link == ClusterApp::Link::identity) {
            if (statmodel.sampling == ClusterApp::Sampling::cohort || statmodel.sampling == ClusterApp::Sampling::open_cohort) {
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

        // switch between ICC and parameter specifications
        statmodel.cov_pars[2] = mean_individual_variance(false);
        double ind_level_error = statmodel.cov_pars[2];
        double cl_level_error;
        if (option.use_icc_for_non_gaussian) {
            if (statmodel.sampling == ClusterApp::Sampling::cohort) {
                statmodel.cov_pars[3] = statmodel.ixx_pars[2] * statmodel.cov_pars[2] / (1 - statmodel.ixx_pars[2]);
                ind_level_error += statmodel.cov_pars[3];
            }
            cl_level_error = statmodel.ixx_pars[0] * ind_level_error / (1 - statmodel.ixx_pars[0]);
            if (statmodel.covariance == ClusterApp::Covariance::nested_exchangeable) {
                statmodel.cov_pars[0] = statmodel.ixx_pars[1] * cl_level_error;
                statmodel.cov_pars[1] = ( 1- statmodel.ixx_pars[1])  * cl_level_error;
            }
            else  {
                statmodel.cov_pars[0] = cl_level_error;
            }
        }
        else {
            if (statmodel.sampling == ClusterApp::Sampling::cohort) {                
                ind_level_error += statmodel.cov_pars[3];
            }
            cl_level_error = statmodel.cov_pars[0];
            if (statmodel.covariance == ClusterApp::Covariance::nested_exchangeable) {
                cl_level_error += statmodel.cov_pars[1];
            }
            statmodel.ixx_pars[0] = cl_level_error / (cl_level_error + ind_level_error);
            if (statmodel.covariance == ClusterApp::Covariance::nested_exchangeable) {
                statmodel.ixx_pars[1] = statmodel.cov_pars[1] / cl_level_error;
            }
            if (statmodel.sampling == ClusterApp::Sampling::cohort) {
                statmodel.ixx_pars[2] = statmodel.cov_pars[3] / ind_level_error;
            }
        }

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
        } else if (statmodel.sampling == ClusterApp::Sampling::open_cohort) {
            double tau3 = statmodel.cov_pars[3];
            tau3 = tau3 / mean_n;
            theta.push_back(tau3);
            theta.push_back(1 - statmodel.cov_pars[4]);
        }
        if (statmodel.family == ClusterApp::Family::beta || statmodel.family == ClusterApp::Family::gamma) {
            (*model).model.data.set_var_par(1 - statmodel.cov_pars[2]);
        }
    }

    (*model).update_beta(beta);
    (*model).update_theta(theta);
    if (option.log) {
        logger.AddLog("[%05d] [%s] Update beta: \n", ImGui::GetFrameCount(), logger.cat[0]);
        ClusterApp::AddVectorToLog(beta, logger);
        logger.AddLog("[%05d] [%s] Update theta: \n", ImGui::GetFrameCount(), logger.cat[0]);
        ClusterApp::AddVectorToLog(theta, logger);
    }
    (*model).matrix.W.update();
}

void ClusterApp::glmmModel::update_model_data(const Eigen::ArrayXXd& data) {
    if(option.log)logger.AddLog("[%05d] [%s] Update GLMM model data \n", ImGui::GetFrameCount(), logger.cat[0]);
    if (model)model.release();
    model = std::unique_ptr<glmm>(new glmm(formula, data, colnames, family, link));
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
    if (option.log)logger.AddLog("[%05d] [%s] Calculate GLS power \n", ImGui::GetFrameCount(), logger.cat[0]);
    if (model) {
        zcutoff = boost::math::quantile(norm, 1-statmodel.alpha/2);
        Eigen::MatrixXd M = (*model).matrix.information_matrix();
        if (option.log && logger.ShowMatrix) {
            logger.AddLog("[%05d] [%s] Information matrix: \n", ImGui::GetFrameCount(), logger.cat[0]);
            ClusterApp::AddMatrixToLog(M, logger);
        }
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
            if (option.log) {
                logger.AddLog("[%05d] [%s] GLS information matrix not positive definite \n", ImGui::GetFrameCount(), logger.cat[2]);
            }
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
                if (option.log) {
                    logger.AddLog("[%05d] [%s] GLS information matrix (2nd treatment) not positive definite \n", ImGui::GetFrameCount(), logger.cat[2]);
                }
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
                if (option.log) {
                    logger.AddLog("[%05d] [%s] GLS information matrix (interaction) not positive definite \n", ImGui::GetFrameCount(), logger.cat[2]);
                }
                summary.power_12 = 909;
                summary.dof_12 = 0;
                summary.se_12 = 0;
                summary.ci_width_12 = 0;
            }
            
        }
    }
}

void ClusterApp::glmmModel::power_box(ClusterApp::modelSummary& summary) {
    if (option.log)logger.AddLog("[%05d] [%s] Calculate Box correction power \n", ImGui::GetFrameCount(), logger.cat[0]);
    if (statmodel.family == Family::gaussian && statmodel.link == Link::identity) {
        BoxResults res = model->matrix.box();
        int idx = statmodel.include_intercept == 1 ? 1 : 0;
        boost::math::fisher_f fdist(1.0, res.dof[idx]);
        double fcutoff = boost::math::quantile(fdist, 1 - statmodel.alpha);
        double fval = res.test_stat[idx] / res.scale[idx];
        summary.power_box = res.dof[idx] > 1 ? boost::math::cdf(fdist, fval - fcutoff) * 100 : 0.0;
        summary.dof_box = res.dof[idx];

        if (option.two_treatments) {
            boost::math::fisher_f fdist2(1.0, res.dof[idx+1]);
            fcutoff = boost::math::quantile(fdist2, 1 - statmodel.alpha);
            fval = res.test_stat[idx+1] / res.scale[idx+1];
            summary.power_box_2 = res.dof[idx+1] > 1 ? boost::math::cdf(fdist2, fval - fcutoff) * 100 : 0.0;
            summary.dof_box_2 = res.dof[idx+1];

            boost::math::fisher_f fdist3(1.0, res.dof[idx + 2]);
            fcutoff = boost::math::quantile(fdist3, 1 - statmodel.alpha);
            fval = res.test_stat[idx + 2] / res.scale[idx + 2];
            summary.power_box_12 = res.dof[idx + 2] > 1 ? boost::math::cdf(fdist3, fval - fcutoff) * 100 : 0.0;
            summary.dof_box_12 = res.dof[idx + 2];
        }
    }
}

void ClusterApp::glmmModel::power_kr(ClusterApp::modelSummary& summary) {
    if (option.log)logger.AddLog("[%05d] [%s] Calculate Kenward-Roger power \n", ImGui::GetFrameCount(), logger.cat[0]);
    if (model) {
        zcutoff = boost::math::quantile(norm, 1 - statmodel.alpha / 2);
        int valsize = option.two_treatments ? 3 : 1;
        dblvec dofkr(valsize);
        dblvec bvar(valsize);
        dblvec bvar2(valsize);
        int idx = statmodel.include_intercept == 1 ? 1 : 0;
        bool with_improved = false;
        if (statmodel.covariance != Covariance::exchangeable && statmodel.covariance != Covariance::nested_exchangeable) {
            with_improved = true;
            CorrectionData<glmmr::SE::KRBoth> res = (*model).matrix.small_sample_correction<glmmr::SE::KRBoth>();
            dofkr[0] = res.dof(idx) > 1 ? res.dof(idx) : 1.0;
            bvar[0] = res.vcov_beta(idx, idx);
            bvar2[0] = res.vcov_beta_second(idx, idx);
            if (option.log && logger.ShowMatrix) {
                logger.AddLog("[%05d] [%s] Kenward-Roger Matrix \n", ImGui::GetFrameCount(), logger.cat[0]);
                ClusterApp::AddMatrixToLog(res.vcov_beta, logger);
                logger.AddLog("[%05d] [%s] Improved Kenward-Roger Matrix \n", ImGui::GetFrameCount(), logger.cat[0]);
                ClusterApp::AddMatrixToLog(res.vcov_beta_second, logger);
            }
            if (option.two_treatments) {
                dofkr[1] = res.dof(idx+1) > 1 ? res.dof(idx + 1) : 1.0;
                bvar[1] = res.vcov_beta(idx + 1, idx + 1);
                bvar2[1] = res.vcov_beta_second(idx + 1, idx + 1);
                dofkr[2] = res.dof(idx + 2) > 1 ? res.dof(idx + 2) : 1.0;
                bvar[2] = res.vcov_beta(idx + 2, idx + 2);
                bvar2[2] = res.vcov_beta_second(idx + 2, idx + 2);
            }
        }
        else {
            CorrectionData<glmmr::SE::KR> res = (*model).matrix.small_sample_correction<glmmr::SE::KR>();
            dofkr[0] = res.dof(idx) > 1 ? res.dof(idx) : 1.0;
            bvar[0] = res.vcov_beta(idx, idx);
            bvar2[0] = 0;
            if (option.log && logger.ShowMatrix) {
                logger.AddLog("[%05d] [%s] Kenward-Roger Matrix \n", ImGui::GetFrameCount(), logger.cat[0]);
                ClusterApp::AddMatrixToLog(res.vcov_beta, logger);
            }
            
            if (option.two_treatments) {
                dofkr[1] = res.dof(idx + 1) > 1 ? res.dof(idx + 1) : 1.0;
                bvar[1] = res.vcov_beta(idx + 1, idx + 1);
                bvar2[1] = 0;
                dofkr[2] = res.dof(idx + 2) > 1 ? res.dof(idx + 2) : 1.0;
                bvar[2] = res.vcov_beta(idx + 2, idx + 2);
                bvar2[2] = 0;
            }
        }        
        boost::math::students_t dist(dofkr[0]);
        double tval, tcutoff, tval_sat, tval2;
        tval_sat = abs(statmodel.te_pars[0] / summary.se);
        tcutoff = boost::math::quantile(dist, 1 - statmodel.alpha / 2);
        summary.dof_kr = dofkr[0];
        summary.power_sat = dofkr[0] > 1 ? boost::math::cdf(dist, tval_sat - tcutoff) * 100 : 0.0;
        summary.ci_width_sat = tcutoff * summary.se;
        if (option.log)logger.AddLog("[%05d] [%s] Calculate Satterthwaite power SE: %.3f te: %.3f \n", ImGui::GetFrameCount(), logger.cat[0], summary.se, statmodel.te_pars[0]);
        if (!isnan(bvar[0]) && bvar[0] > 0 && bvar[0] < 1e5) {
            tval = abs(statmodel.te_pars[0] / sqrt(bvar[0]));
            summary.power_kr = dofkr[0] > 1 ? boost::math::cdf(dist, tval - tcutoff) * 100 : 0.0;
            summary.se_kr = sqrt(bvar[0]);
            summary.ci_width_kr = tcutoff * summary.se_kr;  
        }
        else {
            if (option.log) {
                logger.AddLog("[%05d] [%s] Kenward-Roger Matrix not positive definite \n", ImGui::GetFrameCount(), logger.cat[2]);
            }
            summary.power_kr = 909;
            summary.dof_kr = 0;
            summary.se_kr = 0;
            summary.ci_width_kr = 0;
        }
        if (with_improved) {
            if (!isnan(bvar2[0]) && bvar2[0] > 0 && bvar2[0] < 1e5) {
                tval2 = abs(statmodel.te_pars[0] / sqrt(bvar2[0]));
                summary.power_kr2 = dofkr[0] > 1 ? boost::math::cdf(dist, tval2 - tcutoff) * 100 : 0.0;
                summary.se_kr2 = sqrt(bvar2[0]);
                summary.ci_width_kr2 = tcutoff * summary.se_kr2;
            }
            else {
                if (option.log) {
                    logger.AddLog("[%05d] [%s] Improved Kenward-Roger Matrix not positive definite or other error \n", ImGui::GetFrameCount(), logger.cat[2]);
                }
                summary.power_kr2 = 909;
                summary.se_kr2 = 0;
                summary.ci_width_kr2 = 0;
            }
        }

        if (option.two_treatments) {
            boost::math::students_t dist2(dofkr[1]);
            summary.dof_kr_2 = dofkr[1];
            tval_sat = abs(statmodel.te_pars[0] / summary.se_2);
            summary.power_sat_2 = dofkr[1] > 1 ? boost::math::cdf(dist2, tval_sat - tcutoff) * 100 : 0.0;
            summary.ci_width_sat_2 = tcutoff * summary.se_2;
            if (!isnan(bvar[1]) && bvar[1] >= 0) {
                tval = abs(statmodel.te_pars[1] / sqrt(bvar[1]));
                summary.power_kr_2 = dofkr[1] > 1 ? boost::math::cdf(dist2, tval - tcutoff) * 100 : 0.0;
                summary.se_kr_2 = sqrt(bvar[1]);
                summary.ci_width_kr_2 = tcutoff * summary.se_kr_2;
            }
            else {
                if (option.log) {
                    logger.AddLog("[%05d] [%s] Kenward-Roger Matrix (2nd treatment) not positive definite \n", ImGui::GetFrameCount(), logger.cat[2]);
                }
                summary.power_kr_2 = 909;
                summary.dof_kr_2 = 0;
                summary.se_kr_2 = 0;
                summary.ci_width_kr_2 = 0;
            }
            if (with_improved) {
                if (!isnan(bvar2[1]) && bvar2[1] > 0) {
                    tval2 = abs(statmodel.te_pars[1] / sqrt(bvar2[1]));
                    summary.power_kr2_2 = dofkr[1] > 1 ? boost::math::cdf(dist2, tval2 - tcutoff) * 100 : 0.0;
                    summary.se_kr2_2 = sqrt(bvar2[1]);
                    summary.ci_width_kr2_2 = tcutoff * summary.se_kr2_2;
                }
                else {
                    if (option.log) {
                        logger.AddLog("[%05d] [%s] Improved Kenward-Roger Matrix (2nd treatment) not positive definite or other error \n", ImGui::GetFrameCount(), logger.cat[2]);
                    }
                    summary.power_kr2_2 = 909;
                    summary.se_kr2_2 = 0;
                    summary.ci_width_kr2_2 = 0;
                }
            }

            boost::math::students_t dist3(dofkr[2]);
            tval_sat = abs(statmodel.te_pars[2] / summary.se_12);
            summary.dof_kr_12 = dofkr[2];
            summary.power_sat_12 = dofkr[2] > 1 ? boost::math::cdf(dist3, tval_sat - tcutoff) * 100 : 0.0;
            summary.ci_width_sat_12 = tcutoff * summary.se_12;
            if (!isnan(bvar[2]) && bvar[2] >= 0) {
                tval = abs(statmodel.te_pars[2] / sqrt(bvar[2]));
                summary.power_kr_12 = dofkr[2] > 1 ? boost::math::cdf(dist3, tval - tcutoff) * 100 : 0.0;
                summary.se_kr_12 = sqrt(bvar[2]);
                summary.ci_width_kr_12 = tcutoff * summary.se_kr_12;
            }
            else {
                if (option.log) {
                    logger.AddLog("[%05d] [%s] Kenward-Roger Matrix (interaction) not positive definite or other error \n", ImGui::GetFrameCount(), logger.cat[2]);
                }
                summary.power_kr_12 = 909;
                summary.dof_kr_12 = 0;
                summary.se_kr_12 = 0;
                summary.ci_width_kr_12 = 0;
            }
            if (with_improved) {
                if (!isnan(bvar2[2]) && bvar2[2] > 0) {
                    tval2 = abs(statmodel.te_pars[2] / sqrt(bvar2[2]));
                    summary.power_kr2_12 = dofkr[2] > 1 ? boost::math::cdf(dist3, tval2 - tcutoff) * 100 : 0.0;
                    summary.se_kr2_12 = sqrt(bvar2[2]);
                    summary.ci_width_kr2_12 = tcutoff * summary.se_kr2_12;
                }
                else {
                    if (option.log) {
                        logger.AddLog("[%05d] [%s] Improved Kenward-Roger Matrix (interaction) not positive definite or other error \n", ImGui::GetFrameCount(), logger.cat[2]);
                    }
                    summary.power_kr2_12 = 909;
                    summary.se_kr2_12 = 0;
                    summary.ci_width_kr2_12 = 0;
                }
            }
            
        }
    }
}

void ClusterApp::glmmModel::power_bw(ClusterApp::modelSummary& summary) {
    if (model) {
        zcutoff = boost::math::quantile(norm, 1 - statmodel.alpha / 2);
        Eigen::MatrixXd M = (*model).matrix.information_matrix();
        M = M.llt().solve(Eigen::MatrixXd::Identity(M.rows(), M.cols()));
        int idx = statmodel.include_intercept == 1 ? 1 : 0;
        double zval;
        double bvar = M(idx, idx);        
        double dofbw = dof - (*model).model.linear_predictor.P();
        dofbw = dofbw > 1 ? dofbw : 1.0;
        boost::math::students_t dist(dofbw);
        double tcutoff = boost::math::quantile(dist, 0.975);
        if (option.log)logger.AddLog("[%05d] [%s] Calculate GLS-BW power, SE %.3f DoF %.1f te %.3f \n", ImGui::GetFrameCount(), logger.cat[0], (float)sqrt(bvar), (float)dofbw, statmodel.te_pars[0]);
        if (!isnan(bvar) && bvar > 0) {
            zval = abs(statmodel.te_pars[0] / sqrt(bvar));                         
            summary.power_bw = dofbw == 1 ? 0.0 : boost::math::cdf(dist, zval - tcutoff) * 100;
            summary.dof_bw = dof - (*model).model.linear_predictor.P();
            summary.ci_width_bw = tcutoff * summary.se;
        }
        else {
            if (option.log) {
                logger.AddLog("[%05d] [%s] (B-W) GLS Matrix not positive definite \n", ImGui::GetFrameCount(), logger.cat[2]);
            }
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
    if (option.log)logger.AddLog("[%05d] [%s] Optimum design algorithm \n", ImGui::GetFrameCount(), logger.cat[0]);
    if (model) {
        Eigen::VectorXd C = Eigen::VectorXd::Zero((*model).model.linear_predictor.P());
        int idx = statmodel.include_intercept == 1 ? 1 : 0;
        C(idx) = 1; //statmodel.c_vals[0];
        if (option.two_treatments) {
            C(idx+1) = statmodel.c_vals[1];
            C(idx + 2) = statmodel.c_vals[2];
        }
        Eigen::ArrayXd fix_weights = model->model.data.weights;
        model->model.data.set_weights(fix_weights.size()); // set weights to 1
        model->matrix.W.update();
        Eigen::ArrayXd weights = (*model).optim.optimum_weights(N,C,1e-6);
        model->model.data.set_weights(fix_weights);
        model->matrix.W.update();
        if (weights.size() != optimal_weights.size()) optimal_weights.resize(weights.size());
        for (int i = 0; i < weights.size(); i++) optimal_weights[i] = weights(i);
        if (option.log) {
            logger.AddLog("[%05d] [%s] Optimum design weights: \n", ImGui::GetFrameCount(), logger.cat[0]);
            ClusterApp::AddVectorToLog(optimal_weights, logger, 10);
        }
    }
}

float ClusterApp::glmmModel::individual_n() {
    int n0 = 0;
    int n1 = 0;
    for (int i = 0; i < designs.sequences; i++) {
        for (int t = 0; t < designs.time; t++) {
            if (*designs.active(i, t)) {
                if (*designs.intervention(i, t)) {
                    n1 += *designs.n(i, t) * (*designs.n_clusters(i));
                }
                else {
                    n0 += *designs.n(i, t) * (*designs.n_clusters(i));
                }
            }
        }
    }
    float r = (float)n1 / (float)n0;
    float m = 2 * r * n0 / (1 + r);
    if (option.log)logger.AddLog("[%05d] [%s] Calculate equivalent sample size, ratio %.2f, m %.1f \n", ImGui::GetFrameCount(), logger.cat[0], r, m);
    return m;
}


double ClusterApp::glmmModel::design_effect() {
    std::vector<int> colsums(designs.time);
    std::vector<int> rowsums(designs.sequences);
    
    double deffr = 1;
    double r = 1;
    double m = designs.mean_n();
    double L = (double)designs.sequences;   
    double deffc = 1 + (m - 1) * statmodel.ixx_pars[0];   

    if (designs.time > 1) {
        double ixx_par = statmodel.covariance == ClusterApp::Covariance::nested_exchangeable ? statmodel.ixx_pars[0] * statmodel.ixx_pars[1] : statmodel.ixx_pars[0];
        if (statmodel.sampling == ClusterApp::Sampling::cohort) {
            r = m * ixx_par + (1 - statmodel.ixx_pars[0]) * statmodel.ixx_pars[2];
        }
        else {
            r = m * ixx_par;
        }
        r *= (1 / deffc);

        float B, C, D;
        for (int i = 0; i < designs.sequences; i++) {
            for (int t = 0; t < designs.time; t++) {
                if (*designs.active(i, t) && *designs.intervention(i, t)) {
                    colsums[t]++;
                    rowsums[i]++;
                }
            }
        }

        B = 0;
        C = 0;
        D = 0;
        for (int i = 0; i < designs.sequences; i++) {
            B += rowsums[i];
            C += rowsums[i] * rowsums[i];
        }
        for (int t = 0; t < designs.time; t++) {
            D += colsums[t] * colsums[t];
        }

        deffr = L * L * (1 - r) * (1 + designs.time * r);
        double deffr_denom = 4 * (L * B - D + r * (B * B + L * designs.time * B - designs.time * D - L * C));
        deffr *= 1 / deffr_denom;
    }
    float design_effect = deffc * deffr;
    if (option.log)logger.AddLog("[%05d] [%s] Calculate design effect: C %.3f R %.3f DE %.3f \n", ImGui::GetFrameCount(), logger.cat[0], deffc, deffr, design_effect);
    return design_effect;
}

void ClusterApp::glmmModel::power_de(ClusterApp::modelSummary& summary, int type) {
    //TYPES: 1 = glm model variance
    //TYPES: 2 = chi-squared (applicable to binary outcome)
    //TYPES: 3 = Fisher's exact test
    //TYPES: 4 = t-test difference in means - not needed as equivalent to linear model
    // also investigate count data and others.
    summary.individual_n = individual_n();
    summary.individual_n *= (1.0 / designs.time);
    summary.design_effect = design_effect();
    if (option.log)logger.AddLog("[%05d] [%s] Calculate design effect power \n", ImGui::GetFrameCount(), logger.cat[0]);

    if (type == 0) {
        double individual_var = mean_individual_variance(false);
        if (!isnan(individual_var) && individual_var > 0) {
            zcutoff = boost::math::quantile(norm, 1 - statmodel.alpha / 2);
            summary.individual_var = mean_individual_variance(false);
            summary.individual_se = sqrt(2 * summary.individual_var / summary.individual_n);
            summary.se_de = sqrt(2 * summary.individual_var * summary.design_effect / summary.individual_n);
            double zval = abs(statmodel.te_pars[0] / summary.se_de);
            summary.power_de = boost::math::cdf(norm, zval - zcutoff) * 100;
            summary.ci_width_de = zcutoff * summary.se_de;
        }
        else {
            if (option.log)logger.AddLog("[%05d] [%s] Design effect, individual power is nan or zero: %.3f \n", ImGui::GetFrameCount(), logger.cat[2], (float)individual_var);
            summary.individual_var = 0;
            summary.individual_n = 0;
            summary.design_effect = design_effect();
            summary.individual_se = 0;
            summary.se_de = 0;
            summary.power_de = 909;
            summary.ci_width_de = 0;
        }
    }
    else if(type == 1) {
        // chi-squared
        std::pair<double, double> ymean = mean_outcome();
        float ratio = designs.randomisation_ratio(1, false);
        if (ymean.first <= 0 || ymean.first >= 1 || ymean.second <= 0 || ymean.second >= 1 || ratio == 0 || ratio == 1) {
            if (option.log)logger.AddLog("[%05d] [%s] Design effect Chi-sq, complete case error [%.2f %.2f] ratio: %.2f  \n", ImGui::GetFrameCount(), logger.cat[2], (float)ymean.first, (float)ymean.second, (float)ratio);
            summary.individual_var = 0;
            summary.individual_n = 0;
            summary.design_effect = design_effect();
            summary.individual_se = ratio;
            summary.se_de = ymean.first;
            summary.power_de = 909;
            summary.ci_width_de = ymean.second;
        }
        else {
            double n = summary.individual_n / summary.design_effect;
            double del = ymean.second - ymean.first;
            //double non_central_param = del * del * n * ratio * (1-ratio) * (1 / (ymean.first * (1-ymean.first)));
            // for test of independence above
            double non_central_param = del * del * n / (2 * ymean.first);
            boost::math::non_central_chi_squared dist(1, non_central_param);
            boost::math::chi_squared chi_dist(1);
            double k = boost::math::quantile(chi_dist, 1 - statmodel.alpha / 2);                        
            summary.power_de = (1 - boost::math::cdf(dist, k)) * 100;
        }       
        
    }

    
    //case 3:
    //{
    //    // fishers exact test to add another time - 
    //    // proposed strategy is to define the hypergeometric distribution under the null,
    //    // and then calculate probabilities for possible y counts under alternative and 
    //    // take weighted average using binomial probablities, but likely very slow for 
    //    // larger sample sizes. Put a cap on sample sizes.
    //}
    

    // get sample sizes and design effects - also update summary to include new variables
    

    
}

std::vector<int> ClusterApp::glmmModel::round_weights(std::vector<float> w, int n) {
    if (option.log)logger.AddLog("[%05d] [%s] Rounding weights \n", ImGui::GetFrameCount(), logger.cat[0]);
    int total = w.size();
    std::vector<double> totals(total);
    std::vector<double> rem(total);
    int remn = n;
    for (int i = 0; i < total; i++) {
        double rn = n * w[i];
        totals[i] = std::floor(rn);
        rem[i] = n * w[i] - totals[i];
        remn -= totals[i];
    }
    while (remn > 0) {
        auto largest_remainder_iter = std::max_element(rem.begin(), rem.end());
        auto largest_remainder_idx = std::distance(rem.begin(), largest_remainder_iter);
        totals[largest_remainder_idx] += 1;
        rem[largest_remainder_idx] = 0;
        remn -= 1;
    }
    std::vector<int> result(total);
    for (int i = 0; i < total; i++)result[i] = (int)totals[i];
    return result;
}

std::vector<double> ClusterApp::glmmModel::sim_data() {
    if (option.log)logger.AddLog("[%05d] [%s] Simulating data \n", ImGui::GetFrameCount(), logger.cat[0]);
    Eigen::VectorXd re = model->model.covariance.sim_re();
    Eigen::ArrayXd xb = model->model.xb();
    re += xb.matrix();
    Eigen::VectorXd mu = glmmr::maths::mod_inv_func(re, model->model.family.link);
    std::vector<double> sim_data(mu.data(), mu.data()+mu.size());
    return sim_data;
}