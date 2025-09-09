#include "clusterclasses.h"

std::pair<bool, bool> ClusterApp::statisticalModel::check() {
    CRC crc;
    CRC crc2;
    crc(static_cast<int>(family));
    crc(static_cast<int>(link));
    crc(static_cast<int>(covariance));
    crc(static_cast<int>(linearpredictor));
    crc(static_cast<int>(sampling));
    crc(static_cast<int>(powertype));
    crc(mean_size);
    crc(cv_size);
    crc(cv_size_within);
    crc(option.use_exact_cell_sizes);
    crc(option.heterogeneous_te);
    crc(option.one_arm_clustering);
    crc(include_intercept);    
    crc(target_power);
    for (int i = 0; i < te_pars.size(); i++)crc2(te_pars[i]);
    for (int i = 0; i < treat_control_means.size(); i++)crc2(treat_control_means[i]);
    for (int i = 0; i < ixx_pars.size(); i++)crc2(ixx_pars[i]);
    for (int i = 0; i < cov_pars.size(); i++)crc2(cov_pars[i]);
    for (int i = 0; i < beta_pars.size(); i++)crc2(beta_pars[i]);
    crc2(alpha);
    // crc2(quantile);
    int new_val = crc.get();
    int new_val2 = crc2.get();
    std::pair<bool, bool> has_changed;
    has_changed.first = crc_val != new_val;
    has_changed.second = crc_val_pars != new_val2;
    if (has_changed.first)crc_val = crc.get();
    if (has_changed.second)crc_val_pars = crc2.get();
    return has_changed;
}

void ClusterApp::statisticalModel::update_beta(ClusterApp::design& design) {
    if (linearpredictor == ClusterApp::LinearPredictor::time_fixed_effects) {
        int T = design.active_time_periods();
        if (beta_pars.size() != T) {
            int difft = T - beta_pars.size();
            if (difft < 0) {
                beta_pars.resize(T);
            }
            else {
                for (int i = 0; i < difft; i++)beta_pars.push_back(0);
            }
        }
    }
    else {
        int totalcl = design.total_clusters();
        if (include_intercept == 1)totalcl++;
        if (beta_pars.size() != totalcl) {
            int difft = totalcl - beta_pars.size();
            if (difft < 0) {
                beta_pars.resize(totalcl);
            }
            else {
                for (int i = 0; i < difft; i++)beta_pars.push_back(0);
            }
        }
    }
}

void ClusterApp::statisticalModel::set_beta_random(const double m, const double s) {
    std::default_random_engine generator; 
    std::mt19937 rng(generator());
    rng.seed(time(0));
    std::normal_distribution<double> distribution(m , s);
    for (int i = 0; i < beta_pars.size(); i++) {
        beta_pars[i] = distribution(rng);
    }
}
