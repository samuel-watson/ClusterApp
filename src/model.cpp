#include "clusterapp.h"

namespace ClusterApp {


    void RenderModel(ClusterApp::design& design, ClusterApp::statisticalModel& model, ClusterApp::options& option) {
        ImGui::Begin("Statistical Model");

        ImGui::TextWrapped("Set the statistical model using the two option trees below. The first sets the type of mixed model, including covariance function, the second specifies the parameter values, including \
the treatment effect.");

        const char* family_items[] = { "Gaussian", "Binomial", "Poisson", "Beta", "Gamma" };
        const char* gaussian_link_items[] = { "Identity", "Log" };
        const char* binomial_link_items[] = { "Identity", "Log", "Logit", "Probit" };
        const char* poisson_link_items[] = { "Identity", "Log" };
        const char* beta_link_items[] = { "Identity", "Log", "Logit", "Probit" };
        const char* gamma_link_items[] = { "Identity", "Log", "Inverse" };
        static int family_item_current = 0;
        static int link_item_current = 0;
        static int structure_sampling = 0;
        const char* cluster_covariance_items[] = { "Exchangeable", "Nested exchangeable", "Autoregressive", "Exponential", "Squared exponential" };
        const char* individual_covariance_items[] = { "Exchangeable", "Autoregressive", "Exponential", "Squared exponential" };
        static int cl_cov_item_current = 0;
        static int ind_cov_item_current = 0;
        const char* linpred_items[] = { "Time period fixed effects", "Cluster-specific linear trends" };
        static int linpred_item_current = 0;
        static float control_mean = 0.5;
        static float treatment_mean = 0.5;

        if (ImGui::TreeNode("Statistical Model")) {
            ImGui::SetNextItemWidth(200);
            ImGui::Combo("Family", &family_item_current, family_items, IM_ARRAYSIZE(family_items));

            switch (family_item_current) {
            case 0:
            {
                model.family = ClusterApp::Family::gaussian;
                ImGui::SetNextItemWidth(200);
                ImGui::Combo("Link", &link_item_current, gaussian_link_items, IM_ARRAYSIZE(gaussian_link_items));
                switch (link_item_current) {
                case 0:
                    model.link = ClusterApp::Link::identity;
                    break;
                case 1:
                    model.link = ClusterApp::Link::log;
                    break;
                }
                break;
            }
            case 1:
            {
                model.family = ClusterApp::Family::binomial;
                ImGui::SetNextItemWidth(200);
                ImGui::Combo("Link", &link_item_current, binomial_link_items, IM_ARRAYSIZE(binomial_link_items));
                switch (link_item_current) {
                case 0:
                    model.link = ClusterApp::Link::identity;
                    break;
                case 1:
                    model.link = ClusterApp::Link::log;
                    break;
                case 2:
                    model.link = ClusterApp::Link::logit;
                    break;
                case 3:
                    model.link = ClusterApp::Link::probit;
                    break;
                }
                break;
            }
            case 2:
            {
                model.family = ClusterApp::Family::poisson;
                ImGui::SetNextItemWidth(200);
                ImGui::Combo("Link", &link_item_current, poisson_link_items, IM_ARRAYSIZE(poisson_link_items));
                switch (link_item_current) {
                case 0:
                    model.link = ClusterApp::Link::identity;
                    break;
                case 1:
                    model.link = ClusterApp::Link::log;
                    break;
                }
                break;
            }
            case 3:
            {
                model.family = ClusterApp::Family::beta;
                ImGui::SetNextItemWidth(200);
                ImGui::Combo("Link", &link_item_current, beta_link_items, IM_ARRAYSIZE(beta_link_items));
                switch (link_item_current) {
                case 0:
                    model.link = ClusterApp::Link::identity;
                    break;
                case 1:
                    model.link = ClusterApp::Link::log;
                    break;
                case 2:
                    model.link = ClusterApp::Link::logit;
                    break;
                case 3:
                    model.link = ClusterApp::Link::probit;
                    break;
                }
                break;
            }
            case 4:
            {
                model.family = ClusterApp::Family::gamma;
                ImGui::SetNextItemWidth(200);
                ImGui::Combo("Link", &link_item_current, gamma_link_items, IM_ARRAYSIZE(gamma_link_items));
                switch (link_item_current) {
                case 0:
                    model.link = ClusterApp::Link::identity;
                    break;
                case 1:
                    model.link = ClusterApp::Link::log;
                    break;
                case 2:
                    model.link = ClusterApp::Link::inverse;
                    break;
                }
                break;
            }
            }


            ImGui::Text("Sampling");
            ImGui::RadioButton("Cross-section", &structure_sampling, 0); ImGui::SameLine();
            ImGui::RadioButton("Closed cohort", &structure_sampling, 1); ImGui::SameLine();
            ImGui::RadioButton("Open cohort", &structure_sampling, 2); 
            model.sampling = static_cast<ClusterApp::Sampling>(structure_sampling + 1);
            ImGui::Text("Covariance");

            int options_size = 1;
            if (design.time > 1)options_size = IM_ARRAYSIZE(cluster_covariance_items);

            ImGui::SetNextItemWidth(200);
            ImGui::Combo("Cluster-level", &cl_cov_item_current, cluster_covariance_items, options_size); ImGui::SameLine(); HelpMarker(
                "Exhangeable means a single cluster-level random effect, nested exchangeable also includes a cluster-period random effect. The other functions describe the within-cluster temporal covariance.");
            if (structure_sampling == 1) {
                ImGui::SetNextItemWidth(200);
                ImGui::Combo("Individual-level", &ind_cov_item_current, individual_covariance_items, IM_ARRAYSIZE(individual_covariance_items));
            }
            else if (structure_sampling == 2) {
                ind_cov_item_current = 0;
                ImGui::SetNextItemWidth(200);
                ImGui::DragFloat("Individual replacement rate", &model.cov_pars[4], 0.01f, 0.0f, 1.0f, "%.3f", ImGuiSliderFlags_None); ImGui::SameLine(); HelpMarker(
                    "This is the proportion of individuals that dropout each time period and are replaced in the cohort. A value of 1 is a closed cohort and a value of 0 is a cross-sectional study.");
            }

            model.covariance = static_cast<ClusterApp::Covariance>(cl_cov_item_current + 1);
            model.ind_covariance = static_cast<ClusterApp::IndividualCovariance>(ind_cov_item_current + 1);

            ImGui::Text("Linear predictor");
            ImGui::SetNextItemWidth(200);
            ImGui::Combo("Fixed effects", &linpred_item_current, linpred_items, IM_ARRAYSIZE(linpred_items));
            ImGui::Text("Include intercept?"); ImGui::SameLine();
            ImGui::RadioButton("No", &model.include_intercept, 0); ImGui::SameLine();
            ImGui::RadioButton("Yes", &model.include_intercept, 1); ImGui::SameLine(); HelpMarker(
                "Excluding the intercept will include all time period fixed effects. It will affect the parameter values that should be entered in non-linear model. This choice will be reflected in the parameter values section below.");

            model.linearpredictor = static_cast<ClusterApp::LinearPredictor>(linpred_item_current + 1);
            ImGui::TreePop();
        }

        model.update_beta(design);

        if (ImGui::TreeNode("Parameters")) {
            if (ImGui::TreeNode("Treatment effects")) {
                ImGui::SetNextItemWidth(200);
                ImGui::DragFloat("Treatment effect parameter", &model.te_pars[0], 0.005f, -FLT_MAX, +FLT_MAX, "%.3f", ImGuiSliderFlags_None); ImGui::SameLine();  HelpMarker(
                    "Use CTRL+Click to directly input the value.");
                if (option.two_treatments) {
                    ImGui::SetNextItemWidth(200);
                    ImGui::DragFloat("Treatment 2 effect parameter", &model.te_pars[1], 0.005f, -FLT_MAX, +FLT_MAX, "%.3f", ImGuiSliderFlags_None);
                    ImGui::SetNextItemWidth(200);
                    ImGui::DragFloat("Interaction effect parameter", &model.te_pars[2], 0.005f, -FLT_MAX, +FLT_MAX, "%.3f", ImGuiSliderFlags_None);
                }
                ImGui::TreePop();
            }

            if (ImGui::TreeNode("Covariance parameters")) {

                if (!(family_item_current == 0 && link_item_current == 0)) {
                    ImGui::TextWrapped("You can set the variance terms for non-Gaussian models using the ICC and related statistics, or by setting the parameter values directly."); ImGui::SameLine();
                    HelpMarker("Use of ICC, CAC, and IAC values for non-Gaussian-identity models uses the mean individual-level variance, which is approximated using the GLM weights, to convert to covariance parameter values.");
                    ImGui::Checkbox("Use ICC?", &option.use_icc_for_non_gaussian);
                }

                if ((family_item_current == 0 && link_item_current == 0) || option.use_icc_for_non_gaussian) {
                    ImGui::SetNextItemWidth(200);
                    ImGui::DragFloat("ICC", &model.ixx_pars[0], 0.001f, 0.0f, 1.0f, "%.3f", ImGuiSliderFlags_None); ImGui::SameLine();  HelpMarker(
                        "Intra-class correlation coefficient.");
                    if (cl_cov_item_current == 1) {
                        ImGui::SetNextItemWidth(200);
                        ImGui::DragFloat("CAC", &model.ixx_pars[1], 0.001f, 0.0f, 1.0f, "%.3f", ImGuiSliderFlags_None); ImGui::SameLine();  HelpMarker(
                            "Cluster autocorrelation coefficient.");
                    }
                    if (cl_cov_item_current == 2) {
                        ImGui::SetNextItemWidth(200);
                        //ImGui::PushFont(unifont);ImGui::PopFont();
                        ImGui::DragFloat("Autoregressive", &model.cov_pars[1], 0.01f, 0.0f, 1.0f, "%.3f", ImGuiSliderFlags_None); ImGui::SameLine(); HelpMarker(
                            "Cluster-period autoregressive parameter.");
                    }
                    if (cl_cov_item_current == 3 || cl_cov_item_current == 4) {
                        ImGui::SetNextItemWidth(200);
                        ImGui::DragFloat("Denominator", &model.cov_pars[1], 0.01f, 0.0f, +FLT_MAX, "%.3f", ImGuiSliderFlags_None); ImGui::SameLine(); HelpMarker(
                            "Cluster-period denominator parameter, e.g. exp(|t-t'|/parameter) for exponential covariance.");
                    }
                    if (structure_sampling == 1 || structure_sampling == 2) {
                        ImGui::SetNextItemWidth(200);
                        ImGui::DragFloat("IAC", &model.ixx_pars[2], 0.001f, 0.0f, 1.0f, "%.3f", ImGuiSliderFlags_None); ImGui::SameLine();  HelpMarker(
                            "Individual autocorrelation coefficient. For open cohorts, set this parameter as if it were a closed cohort.");
                        if (ind_cov_item_current == 1 && structure_sampling != 2) {
                            ImGui::SetNextItemWidth(200);
                            ImGui::DragFloat("Autoregressive (individual)", &model.cov_pars[4], 0.01f, 0.0f, 1.0f, "%.3f", ImGuiSliderFlags_None); ImGui::SameLine();  HelpMarker(
                                "Individual autoregressive parameter.");
                        }
                        if ((ind_cov_item_current == 2 || ind_cov_item_current == 3) && structure_sampling != 2) {
                            ImGui::SetNextItemWidth(200);
                            ImGui::DragFloat("Denominator (individual)", &model.cov_pars[4], 0.01f, 0.0f, +FLT_MAX, "%.3f", ImGuiSliderFlags_None); ImGui::SameLine();  HelpMarker(
                                "Individual denominator parameter, e.g. exp(|t-t'|/parameter) for exponential covariance.");
                        }
                    }
                }
                else {
                    ImGui::SetNextItemWidth(200);
                    ImGui::DragFloat("Cluster-level variance", &model.cov_pars[0], 0.01f, 0.0f, +FLT_MAX, "%.3f", ImGuiSliderFlags_None); ImGui::SameLine();  HelpMarker(
                        "Cluster-level variance parameter.");
                    if (cl_cov_item_current == 1) {
                        ImGui::SetNextItemWidth(200);
                        ImGui::DragFloat("Cluster-period level variance", &model.cov_pars[1], 0.01f, 0.0f, +FLT_MAX, "%.3f", ImGuiSliderFlags_None); ImGui::SameLine(); HelpMarker(
                            "Cluster-period level variance parameter.");
                    }
                    if (cl_cov_item_current == 2) {
                        ImGui::SetNextItemWidth(200);
                        ImGui::DragFloat("Cluster-period autoregressive", &model.cov_pars[1], 0.01f, 0.0f, 1.0f, "%.3f", ImGuiSliderFlags_None); ImGui::SameLine(); HelpMarker(
                            "Cluster-period autoregressive parameter.");
                    }
                    if (cl_cov_item_current == 3 || cl_cov_item_current == 4) {
                        ImGui::SetNextItemWidth(200);
                        ImGui::DragFloat("Cluster-period denominator", &model.cov_pars[1], 0.01f, 0.0f, +FLT_MAX, "%.3f", ImGuiSliderFlags_None); ImGui::SameLine(); HelpMarker(
                            "Cluster-period denominator parameter, e.g. exp(|t-t'|/parameter) for exponential covariance.");
                    }
                    if (structure_sampling == 1 || structure_sampling == 2) {
                        ImGui::SetNextItemWidth(200);
                        ImGui::DragFloat("Individual-level variance", &model.cov_pars[3], 0.01f, 0.0f, +FLT_MAX, "%.3f", ImGuiSliderFlags_None); ImGui::SameLine(); HelpMarker(
                            "Individual-level variance term");
                        if (ind_cov_item_current == 2 && structure_sampling != 2) {
                            ImGui::SetNextItemWidth(200);
                            ImGui::DragFloat("Autoregressive (individual)", &model.cov_pars[4], 0.01f, 0.0f, 1.0f, "%.3f", ImGuiSliderFlags_None); ImGui::SameLine(); HelpMarker(
                                "Individual autoregressive parameter.");
                        }
                        if ((ind_cov_item_current == 3 || ind_cov_item_current == 4) && structure_sampling != 2) {
                            ImGui::SetNextItemWidth(200);
                            ImGui::DragFloat("Denominator (individual)", &model.cov_pars[4], 0.01f, 0.0f, +FLT_MAX, "%.3f", ImGuiSliderFlags_None); ImGui::SameLine();  HelpMarker(
                                "Individual denominator parameter, e.g. exp(|t-t'|/parameter) for exponential covariance.");
                        }
                    }
                    if (family_item_current == 3 || family_item_current == 4) {
                        ImGui::SetNextItemWidth(200);
                        ImGui::DragFloat("Observation-level variance", &model.cov_pars[2], 0.01f, 0.0f, +FLT_MAX, "%.3f", ImGuiSliderFlags_None); ImGui::SameLine();  HelpMarker(
                            "Observation-level variance term.");
                    }

                }
                ImGui::TreePop();
            }

            if (ImGui::TreeNode("Fixed effect parameters")) {
                // add set all zero, random, constant
                ImGui::Text("Set value defaults");
                if (ImGui::SmallButton("Set all zero")) {
                    for (int i = 0; i < model.beta_pars.size(); i++) {
                        model.beta_pars[i] = 0;
                    }
                }

                static float beta_m = 0;
                static float beta_s = 1;
                static float beta_c = 0;
                if (ImGui::SmallButton("Set random normal")) {
                    model.set_beta_random(beta_m, beta_s);
                }
                ImGui::SameLine();
                ImGui::SetNextItemWidth(100);
                ImGui::DragFloat("Mean", &beta_m, 0.01f, -FLT_MAX, +FLT_MAX, "%.2f", ImGuiSliderFlags_None); ImGui::SameLine();
                ImGui::SetNextItemWidth(100);
                ImGui::DragFloat("sd", &beta_s, 0.01f, -FLT_MAX, +FLT_MAX, "%.2f", ImGuiSliderFlags_None);

                if (ImGui::SmallButton("Set constant")) {
                    for (int i = 0; i < model.beta_pars.size(); i++) {
                        model.beta_pars[i] = beta_c;
                    }
                }
                ImGui::SameLine();
                ImGui::SetNextItemWidth(100);
                ImGui::DragFloat("sd", &beta_c, 0.01f, -FLT_MAX, +FLT_MAX, "%.2f", ImGuiSliderFlags_None);

                ImGui::Dummy(ImVec2(20, 20));

                if (model.include_intercept == 1) {
                    ImGui::SetNextItemWidth(200);
                    ImGui::DragFloat("Intercept", &(model.beta_pars[0]), 0.01f, -FLT_MAX, +FLT_MAX, "%.2f", ImGuiSliderFlags_None);
                }
                if (linpred_item_current == 0) {
                    for (int t = 0; t < (design.time - 1); t++) {
                        ImGui::Text("Time period effects:");
                        ImGui::SetNextItemWidth(200);
                        int shift_idx = model.include_intercept == 0 ? 0 : 1;
                        ImGui::DragFloat(int_to_char(t + 1), &(model.beta_pars[t + shift_idx]), 0.01f, -FLT_MAX, +FLT_MAX, "%.2f", ImGuiSliderFlags_None);
                    }
                    if (model.include_intercept == 0) {
                        ImGui::SetNextItemWidth(200);
                        ImGui::DragFloat(int_to_char(design.time), &(model.beta_pars[design.time - 1]), 0.01f, -FLT_MAX, +FLT_MAX, "%.2f", ImGuiSliderFlags_None);
                    }
                }
                else {
                    ImGui::Text("Cluster specific time trend parameters:");
                    for (int l = 0; l < model.beta_pars.size(); l++) {
                        ImGui::SetNextItemWidth(200);
                        ImGui::DragFloat(int_to_char(l + 1), &(model.beta_pars[l]), 0.01f, -FLT_MAX, +FLT_MAX, "%.2f", ImGuiSliderFlags_None);
                    }
                }

                ImGui::TreePop();
            }

            if (!option.two_treatments) {
                if (ImGui::TreeNode("Set using group means")) {
                    ImGui::TextWrapped("Automatically set fixed effect parameters by specifying the mean outcomes in treatment and control groups. Temporal variation is assumed to be zero. For non-linear models, a correction is applied for the random effect. Parameter values can be found above.");


                    switch (model.family) {
                    case ClusterApp::Family::gaussian:
                        ImGui::SetNextItemWidth(150);
                        ImGui::DragFloat("Control group mean", &control_mean, 0.01f, -FLT_MAX, +FLT_MAX, "%.2f", ImGuiSliderFlags_None);
                        ImGui::SetNextItemWidth(150);
                        ImGui::DragFloat("Treatment group mean", &treatment_mean, 0.01f, -FLT_MAX, +FLT_MAX, "%.2f", ImGuiSliderFlags_None);
                        break;
                    case ClusterApp::Family::binomial: case ClusterApp::Family::beta:
                        ImGui::SetNextItemWidth(150);
                        ImGui::DragFloat("Control group mean", &control_mean, 0.01f, 0.0f, 1.0f, "%.2f", ImGuiSliderFlags_None);
                        ImGui::SetNextItemWidth(150);
                        ImGui::DragFloat("Treatment group mean", &treatment_mean, 0.01f, 0.0f, 1.0f, "%.2f", ImGuiSliderFlags_None);
                        break;
                    case ClusterApp::Family::poisson: case ClusterApp::Family::gamma:
                        ImGui::SetNextItemWidth(150);
                        ImGui::DragFloat("Control group mean", &control_mean, 0.01f, 0.0f, +FLT_MAX, "%.2f", ImGuiSliderFlags_None);
                        ImGui::SetNextItemWidth(150);
                        ImGui::DragFloat("Treatment group mean", &treatment_mean, 0.01f, 0.0f, +FLT_MAX, "%.2f", ImGuiSliderFlags_None);
                        break;
                    }


                    if (ImGui::Button("Update")) {
                        float re_adj = model.cov_pars[0];
                        if (model.covariance == ClusterApp::Covariance::nested_exchangeable)re_adj += model.cov_pars[1];
                        if (model.sampling == ClusterApp::Sampling::cohort) re_adj += model.cov_pars[3] / design.mean_n();
                        re_adj *= 0.5;
                        switch (model.link) {
                        case ClusterApp::Link::identity:
                        {
                            if (model.include_intercept == 1) {
                                model.beta_pars[0] = control_mean;
                            }
                            else {
                                for (int l = 0; l < model.beta_pars.size(); l++) {
                                    model.beta_pars[l] = control_mean;
                                }
                            }
                            model.te_pars[0] = treatment_mean - control_mean;
                            break;
                        }
                        case ClusterApp::Link::log:
                        {
                            if (model.include_intercept == 1) {
                                model.beta_pars[0] = log(control_mean) - re_adj;
                            }
                            else {
                                for (int l = 0; l < model.beta_pars.size(); l++) {
                                    model.beta_pars[l] = log(control_mean) - re_adj;
                                }
                            }
                            model.te_pars[0] = log(treatment_mean) - log(control_mean);
                            break;
                        }
                        case ClusterApp::Link::logit:
                        {
                            if (model.include_intercept == 1) {
                                model.beta_pars[0] = log(control_mean / (1 - control_mean)) - re_adj;
                            }
                            else {
                                for (int l = 0; l < model.beta_pars.size(); l++) {
                                    model.beta_pars[l] = log(control_mean / (1 - control_mean)) - re_adj;
                                }
                            }
                            model.te_pars[0] = log(treatment_mean / (1 - treatment_mean)) - log(control_mean / (1 - control_mean));
                            break;
                        }
                        case ClusterApp::Link::probit:
                        {
                            boost::math::normal norm = boost::math::normal(0.0, 1.0);
                            if (model.include_intercept == 1) {
                                model.beta_pars[0] = boost::math::quantile(norm, control_mean) - re_adj;
                            }
                            else {
                                for (int l = 0; l < model.beta_pars.size(); l++) {
                                    model.beta_pars[l] = boost::math::quantile(norm, control_mean) - re_adj;
                                }
                            }
                            model.te_pars[0] = boost::math::quantile(norm, treatment_mean) - model.beta_pars[0];
                            break;
                        }
                        case ClusterApp::Link::inverse:
                        {
                            if (model.include_intercept == 1) {
                                model.beta_pars[0] = 1/control_mean;
                            }
                            else {
                                for (int l = 0; l < model.beta_pars.size(); l++) {
                                    model.beta_pars[l] = 1/control_mean;
                                }
                            }
                            model.te_pars[0] = 1/(treatment_mean - control_mean);
                            break;
                        }

                        }
                    }


                    ImGui::TreePop();
                }
            }


            ImGui::TreePop();
        }
        /*if (option.debug_info) {
            if (ImGui::TreeNode("Model info")) {
                ImGui::Text("Model checksum"); ImGui::SameLine();
                ImGui::Text("%d", model.crc_val);
                ImGui::Text("Parameters checksum"); ImGui::SameLine();
                ImGui::Text("%d", model.crc_val_pars);
                ImGui::TreePop();
            }
        }*/


        ImGui::End();
    };

    
}

