#include "clusterclasses.h"

ClusterApp::modelUpdater::modelUpdater(ClusterApp::design& designs_,
    ClusterApp::statisticalModel& model_,
    ClusterApp::modelSummary& summary_,
    ClusterApp::glmmModel& glmm_) : designs(designs_), model(model_),
    summary(summary_), glmm(glmm_) {
    update_data();
};

Eigen::ArrayXXd ClusterApp::modelUpdater::generate_data() {
    int n_rows = 0;
    for (int i = 0; i < designs.sequences; i++) {
        for (int t = 0; t < designs.time; t++) {
            if (*designs.active(i, t)) {
                n_rows += *designs.n_clusters(i);
            }
        }
    }
    Eigen::ArrayXXd newdata(n_rows, data.cols());
    int cl_number = 1;
    int row_number = 0;
    for (int i = 0; i < designs.sequences; i++) {
        for (int j = 0; j < *designs.n_clusters(i); j++) {
            for (int t = 0; t < designs.time; t++) {
                if (*designs.active(i, t)) {
                    newdata(row_number, 0) = cl_number;
                    newdata(row_number, 1) = t + 1;
                    newdata(row_number, 2) = *designs.n(i, t);
                    newdata(row_number, 3) = *designs.intervention(i, t);
                    newdata(row_number, 4) = *designs.intervention_2(i, t);
                    newdata(row_number, 5) = newdata(row_number, 3) * newdata(row_number, 4);
                    row_number++;
                }
            }
            cl_number++;
        }
    }
    return newdata;
}

void ClusterApp::modelUpdater::update_data() {
    // iterate over the cluster-periods to get the number of rows, then re-iterate to fill in the information
    model.update_beta(designs);
    int n_rows = 0;
    for (int i = 0; i < designs.sequences; i++) {
        for (int t = 0; t < designs.time; t++) {
            if (*designs.active(i, t)) {
                n_rows += *designs.n_clusters(i);
            }
        }
    }
    data.conservativeResize(n_rows, Eigen::NoChange);
    int cl_number = 1;
    int row_number = 0;
    for (int i = 0; i < designs.sequences; i++) {
        for (int j = 0; j < *designs.n_clusters(i); j++) {
            for (int t = 0; t < designs.time; t++) {
                if (*designs.active(i, t)) {
                    data(row_number, 0) = cl_number;
                    data(row_number, 1) = t + 1;
                    data(row_number, 2) = *designs.n(i, t);
                    data(row_number, 3) = *designs.intervention(i, t);
                    data(row_number, 4) = *designs.intervention_2(i, t);
                    data(row_number, 5) = data(row_number, 3) * data(row_number, 4);
                    row_number++;
                }
            }
            cl_number++;
        }
    }
    glmm.update_formula();
    glmm.update_model_data(data);
    glmm.dof = designs.total_cluster_periods();
    if (glmm.option.results)update_summary_statistics();
    if (glmm.option.optimiser) update_optimum();
};

void ClusterApp::modelUpdater::update_formula() {
    glmm.update_formula();
    if (glmm.option.results)update_summary_statistics();
    if(glmm.option.optimiser) update_optimum();
}

void ClusterApp::modelUpdater::update_parameters() {
    model.update_beta(designs);
    glmm.update_parameters();
    if (glmm.option.results) update_summary_statistics();
    if (glmm.option.optimiser) update_optimum();
}

void ClusterApp::modelUpdater::update_summary_statistics() {
    glmm.power(summary);
    glmm.power_kr(summary);
    glmm.power_bw(summary);
    glmm.power_de(summary,de_mode);
    summary.dof = (double)designs.total_n();
    if(!manual_n_optim) summary.total_n = designs.total_n();
}

void ClusterApp::modelUpdater::update_optimum() {
    glmm.optimum(summary.total_n);
    int N = designs.total_clusters();
    if (N != optimum_data.size()) {
        optimum_data.resize(N);
        optimum_n.resize(N);
    }
    int counter = 0;
    int total = glmm.optimal_weights.size();
    std::vector<double> totals(total);
    std::vector<double> rem(total);
    int n = summary.total_n;
    int remn = n;
    for (int i = 0; i < total; i++) {
        double rn = n * glmm.optimal_weights[i];
        totals[i] = std::floor(rn);
        rem[i] = n * glmm.optimal_weights[i] - totals[i];
        remn -= totals[i];
    }
    while (remn > 0) {
        auto largest_remainder_iter = std::max_element(rem.begin(), rem.end());
        auto largest_remainder_idx = std::distance(rem.begin(), largest_remainder_iter);
        totals[largest_remainder_idx]+=1;
        rem[largest_remainder_idx] = 0;
        remn -= 1;
    }
    for (int i = 0; i < N; i++) {
        optimum_data[i].resize(designs.time);
        optimum_n[i].resize(designs.time);
        int seq = designs.seq_by_cluster(i);
        for (int t = 0; t < designs.time; t++) {
            if (*designs.active(seq, t)) {
                if (counter > total-1) {
                    optimum_data[i][t] = total;
                    optimum_n[i][t] = n;
                }
                else {
                    optimum_data[i][t] = glmm.optimal_weights[counter];
                    optimum_n[i][t] =  totals[counter];
                }                
                counter++;
            }
            else {
                optimum_data[i][t] = 0.0;
                optimum_n[i][t] = 0;
            }
        }
    }
}

