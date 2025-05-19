#include "clusterclasses.h"

ClusterApp::modelUpdater::modelUpdater(ClusterApp::design& designs_,
    ClusterApp::statisticalModel& model_,
    ClusterApp::modelSummary& summary_,
    ClusterApp::glmmModel& glmm_, ClusterApp::AppLog& log_) : designs(designs_), model(model_),
    summary(summary_), glmm(glmm_), log(log_) {
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
                    newdata(row_number, 3) = glmm.option.dose_effect ? (* designs.intervention(i, t)) * (*designs.dose(i,t)) : *designs.intervention(i,t);
                    newdata(row_number, 4) = *designs.intervention_2(i, t);
                    newdata(row_number, 5) = newdata(row_number, 3) * newdata(row_number, 4);
                    row_number++;
                }
            }
            cl_number++;
        }
    }
    newdata.col(6) = 1 - newdata.col(3);
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
                    data(row_number, 3) = glmm.option.dose_effect ? (*designs.intervention(i, t)) * (*designs.dose(i, t)) : *designs.intervention(i, t);
                    data(row_number, 4) = *designs.intervention_2(i, t);
                    data(row_number, 5) = data(row_number, 3) * data(row_number, 4);
                    row_number++;
                }
            }
            cl_number++;
        }
    }
    data.col(6) = 1 - data.col(3);
    log.AddLog("[%05d] [%s] %s \n", ImGui::GetFrameCount(), log.cat[0], ("Data update: rows = " + std::to_string(data.rows())).c_str());
    log.AddLog("[%05d] [%s] %s \n", ImGui::GetFrameCount(), log.cat[0], ("Data update: time periods = " + std::to_string(designs.time)).c_str());
    log.AddLog("[%05d] [%s] %s \n", ImGui::GetFrameCount(), log.cat[0], ("Data update: sequences = " + std::to_string(designs.sequences)).c_str());
    log.AddLog("[%05d] [%s] %s \n", ImGui::GetFrameCount(), log.cat[0], ("Data update: dose effect = " + std::to_string(glmm.option.dose_effect)).c_str());
    glmm.update_formula();
    glmm.update_model_data(data);
    glmm.dof = designs.total_cluster_periods();
    if (glmm.option.results)update_summary_statistics();
    if (glmm.option.optimiser) update_optimum();
};

void ClusterApp::modelUpdater::update_formula() {
    glmm.update_formula();
    log.AddLog("[%05d] [%s] %s \n", ImGui::GetFrameCount(), log.cat[0], ("Formula update: "+glmm.formula).c_str());
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
    glmm.power(summary, model.powertype);
    //glmm.power_kr(summary);
    //glmm.power_bw(summary);
    // if(glmm.option.show_box) glmm.power_box(summary);
    //if (model.covariance == Covariance::exchangeable || model.covariance == Covariance::nested_exchangeable) {
    //    glmm.power_de(summary, de_mode);
    //}    
    // summary.dof = (double)designs.total_n();
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

int ClusterApp::modelUpdater::sample_size_search(const bool& clusters, const ClusterApp::PowerType& powertype){
    // save original size
    ClusterApp::modelSummary summary(designs);
    int m = glmm.statmodel.mean_size;
    int J = designs.total_clusters();
    const int lower = clusters ? designs.sequences : 2;
    const int upper = 200;
    int current = clusters ? J : m; 
    int test = current;
    glmm.power(summary, powertype);
    float current_power = summary.power;
    float target = glmm.statmodel.target_power / 100.0;
    int diff = clusters ? (int)J*2 : 100;
    diff = diff < 20 ? 20 : diff;
    bool decrease_n = true;

    std::vector<float> n_weights;
    n_weights.resize(designs.sequences);
    for (int i = 0; i < designs.sequences; i++) {
        n_weights[i] = (float)(*designs.n_clusters(i)) / (float)J;
    }

    Eigen::ArrayXd weights(glmm.model->model.n());
    int iter = 0;
    
    while(diff > 1){
        // direction
        decrease_n = current_power > glmm.statmodel.target_power;

        if(decrease_n){
            test = current - diff < lower ? lower : current - diff;
        } else {
            test = current + diff > upper ? upper : current + diff;
        }

        if(!clusters){
            weights.setConstant((double)test);

            if (glmm.statmodel.family == ClusterApp::Family::gaussian) {
                glmm.model->set_weights(weights);
            }
            else if (glmm.statmodel.family == ClusterApp::Family::binomial) {
                glmm.model->model.data.set_variance(weights);
            }
            else if (glmm.statmodel.family == ClusterApp::Family::poisson) {
                glmm.model->set_offset(weights);
            }
            glmm.model->matrix.W.update();
        } else {
            std::vector<int> new_n_cl = glmm.round_weights(n_weights, test);
            for (int k = 0; k < designs.sequences; k++) {
				*designs.n_clusters(k) = new_n_cl[k];
			}
            glmm.update_model_data(generate_data());
        }
        
        glmm.power(summary, powertype);

        if(abs(summary.power - glmm.statmodel.target_power) < abs(current_power - glmm.statmodel.target_power)){
            current = test;
            current_power = summary.power;
        }
        if (glmm.option.log)log.AddLog("[%05d] [%s] Linesearch (iter %d): current val %d diff %d test %d current power %.2f \n", ImGui::GetFrameCount(), log.cat[0], iter, current, diff, test, current_power);
        diff = (int)diff*0.5;
        iter++;
    }

    if (glmm.option.log && abs(current_power - glmm.statmodel.target_power) > 0.1 )log.AddLog("[%05d] [%s] Linesearch terminated without finding target power, starting difference likely too large \n", ImGui::GetFrameCount(), log.cat[1]);
    
    if(!clusters){
        weights.setConstant((double)m);
        if (glmm.statmodel.family == ClusterApp::Family::gaussian) {
            glmm.model->set_weights(weights);
        }
        else if (glmm.statmodel.family == ClusterApp::Family::binomial) {
            glmm.model->model.data.set_variance(weights);
        }
        else if (glmm.statmodel.family == ClusterApp::Family::poisson) {
            glmm.model->set_offset(weights);
        }
        glmm.model->matrix.W.update();
    } else {
        std::vector<int> new_n_cl = glmm.round_weights(n_weights, J);
        for (int k = 0; k < designs.sequences; k++) {
            *designs.n_clusters(k) = new_n_cl[k];
        }
        glmm.update_model_data(generate_data());
    }
    

    return current;
    
}