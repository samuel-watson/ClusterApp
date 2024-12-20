#include "clusterapp.h"

namespace ClusterApp {


    void RenderMenuBar(ClusterApp::design& design, ClusterApp::statisticalModel& model, ClusterApp::modelChecker& checker, ClusterApp::options& option) {
        ImGui::Begin("Menu Bar");
            static colourPicker colours;
            const char* family_items[] = { "Gaussian", "Binomial", "Poisson", "Beta", "Gamma"};
            const char* gaussian_link_items[] = { "Identity", "Log" };
            const char* binomial_link_items[] = { "Identity", "Log", "Logit", "Probit" };
            const char* poisson_link_items[] = { "Identity", "Log" };
            const char* beta_link_items[] = { "Identity", "Log", "Logit", "Probit" };
            const char* gamma_link_items[] = { "Identity", "Log", "Inverse" };
            static int family_item_current = 0;
            static int link_item_current = 0;
           

           if (ImGui::CollapsingHeader("Parallel")){
                    static int parallel_t = 1;
                    static int parallel_n = 10;
                    static int parallel_J = 1;
                    ImGui::SetNextItemWidth(100);
                    ImGui::InputInt("Time periods", &parallel_t, 1, 10, 0);
                    ImGui::SetNextItemWidth(100);
                    ImGui::InputInt("Observations per cluster-period", &parallel_n, 1, 10, 0);
                    ImGui::SetNextItemWidth(100);
                    ImGui::InputInt("Clusters per arm", &parallel_J, 1, 10, 0);
                    if (ImGui::Button("Set")) {
                        parallel_t = parallel_t < 1 ? 1 : parallel_t;
                        parallel_n = parallel_n < 1 ? 1 : parallel_n;
                        parallel_J = parallel_J < 1 ? 1 : parallel_J;
                        design.set_parallel(parallel_t, parallel_n, parallel_J);
                        checker.updater.update_data();
                    }
           }

           if (ImGui::CollapsingHeader("Parallel with baseline")){
                    static int parallelb_t = 1;
                    static int parallelb_n = 10;
                    static int parallelb_J = 1;
                    ImGui::SetNextItemWidth(100);
                    ImGui::InputInt("Time periods", &parallelb_t, 1, 10, 0);
                    ImGui::SetNextItemWidth(100);
                    ImGui::InputInt("Observations per cluster-period", &parallelb_n, 1, 10, 0);
                    ImGui::SetNextItemWidth(100);
                    ImGui::InputInt("Clusters per arm", &parallelb_J, 1, 10, 0);
                    if (ImGui::Button("Set")) {
                        design.set_parallel_with_baseline(parallelb_t, parallelb_n, parallelb_J);
                        checker.updater.update_data();
                    }
           }

           if (ImGui::CollapsingHeader("Stepped-wedge")){
                    static int wedge_t = 6;
                    static int wedge_n = 10;
                    static int wedge_J = 1;
                    static bool implement = false;
                    ImGui::SetNextItemWidth(100);
                    ImGui::InputInt("Time periods", &wedge_t, 1, 10, 0);
                    ImGui::SetNextItemWidth(100);
                    ImGui::InputInt("Observations per cluster-period", &wedge_n, 1, 10, 0);
                    ImGui::SetNextItemWidth(100);
                    ImGui::InputInt("Clusters per arm", &wedge_J, 1, 10, 0);
                    ImGui::Checkbox("Implementation periods", &implement);
                    if (ImGui::Button("Set")) {
                        design.set_stepped_wedge(wedge_t, wedge_n, wedge_J, implement);
                        checker.updater.update_data();
                    }
           }

           if (ImGui::CollapsingHeader("Crossover")){
                    static int cross_n;
                    static int cross_J;
                    ImGui::Text("Trial details");
                    ImGui::SetNextItemWidth(100);
                    ImGui::InputInt("Observations per cluster-period", &cross_n, 1, 10, 0);
                    ImGui::SetNextItemWidth(100);
                    ImGui::InputInt("Clusters per arm", &cross_J, 1, 10, 0);
                    if (ImGui::Button("Set")) {
                        design.set_crossover(cross_n, cross_J);
                        checker.updater.update_data();
                    }
           }

           if (ImGui::CollapsingHeader("Staircase")){
                    static int staircase_t;
                    static int staircase_n;
                    static int staircase_J;
                    ImGui::SetNextItemWidth(100);
                    ImGui::InputInt("Time periods", &staircase_t, 1, 10, 0);
                    ImGui::SetNextItemWidth(100);
                    ImGui::InputInt("Observations per cluster-period", &staircase_n, 1, 10, 0);
                    ImGui::SetNextItemWidth(100);
                    ImGui::InputInt("Clusters per arm", &staircase_J, 1, 10, 0);
                    if (ImGui::Button("Set")) {
                        design.set_staircase(staircase_t, staircase_n, staircase_J);
                        checker.updater.update_data();
                    }
           }

           if (ImGui::CollapsingHeader("Factorial")){
            if (option.two_treatments) {
                        static int factorial_t;
                        static int factorial_n;
                        static int factorial_J;
                        ImGui::SetNextItemWidth(100);
                        ImGui::InputInt("Time periods", &factorial_t, 1, 10, 0);
                        ImGui::SetNextItemWidth(100);
                        ImGui::InputInt("Observations per cluster-period", &factorial_n, 1, 10, 0);
                        ImGui::SetNextItemWidth(100);
                        ImGui::InputInt("Clusters per arm", &factorial_J, 1, 10, 0);
                        if (ImGui::Button("Set")) {
                            design.set_factorial(factorial_t, factorial_n, factorial_J);
                            checker.updater.update_data();
                        }
                    }
                    else {
                        ImGui::Text("To enable this option, set two treatments.");
                    }
           }

            static int structure_sampling = 0;
            ImGui::Dummy(ImVec2(0.0f, 20.0f));
            ImGui::Text("Sampling Structure");ImGui::SameLine(); HelpMarker(
                "In a cross-sectional design at each measurement occasion a different sample of participants is measured. In a cohort design, participants are repeatedly measured."
            );
            if(design.time == 1){
                structure_sampling = 0;
                ImGui::RadioButton("Cross-section", &structure_sampling, 0); ImGui::SameLine();
            } else {
                ImGui::RadioButton("Cross-section", &structure_sampling, 0); ImGui::SameLine();
                ImGui::RadioButton("Closed cohort", &structure_sampling, 1); ImGui::SameLine();
                ImGui::RadioButton("Open cohort", &structure_sampling, 2); 
            }
            
            model.sampling = static_cast<ClusterApp::Sampling>(structure_sampling + 1);

            if(structure_sampling == 2){
                ImGui::SetNextItemWidth(200);
            ImGui::DragFloat("Individual replacement rate", &model.cov_pars[4], 0.01f, 0.0f, 1.0f, "%.3f", ImGuiSliderFlags_None); ImGui::SameLine(); HelpMarker(
                "This is the proportion of individuals that dropout each time period and are replaced in the cohort. A value of 1 is a closed cohort and a value of 0 is a cross-sectional study.");
            }



            ImGui::Dummy(ImVec2(0.0f, 20.0f));
            ImGui::Text("Cluster sizes");
            static int total_t = 1;
            static int total_n = 10;
            static int mean_cp_size = 50;
            static float cv_cp_size = 0.5;
            static float within_cv = 0.0;
            static int variable_cl_sizes = 0;
            static int variable_cp_sizes = 0; 

            ImGui::Text("Variable cluster sizes"); ImGui::SameLine();
            ImGui::RadioButton("No", &variable_cl_sizes, 0); ImGui::SameLine();
            ImGui::RadioButton("Yes", &variable_cl_sizes, 1);

            if(variable_cl_sizes == 0){
                ImGui::SetNextItemWidth(100);
                ImGui::InputInt("Set size for all cluster-periods", &total_n, 1, 10, 0); ImGui::SameLine();
                ImGui::PushID(200299);
                if (ImGui::Button("Set n")) {
                    if (total_n > 0) {
                        for (int j = 0; j < design.sequences; j++) {
                            for (int t = 0; t < design.time; t++) {
                                *(design.n(j, t)) = total_n;
                            }
                        }
                    }
                }
                ImGui::PopID();
            } else {
                ImGui::TextWrapped("Generate random cluster and/or cluster-period sizes with pre-defined coefficient of variation. Generated sizes of zero or less are set to one. Note that clusters  in the same sequence in the designer will be given the same size. Split the sequences to get different sizes using Edit->Design->Sequences->Split sequences into clusters");
                if(design.time > 1){
                    ImGui::Text("Constant size within clusters"); ImGui::SameLine();
                    ImGui::RadioButton("No", &variable_cp_sizes, 0); ImGui::SameLine();
                    ImGui::RadioButton("Yes", &variable_cp_sizes, 1);
                }
                
                ImGui::SetNextItemWidth(100);
                ImGui::InputInt("Mean cluster size", &mean_cp_size, 1, 50, 0);
                ImGui::SetNextItemWidth(100);
                ImGui::DragFloat("Cluster size c.v.", &cv_cp_size, 0.1f, 0.0f, 10.0f, "%.1f", ImGuiSliderFlags_None);
                if(variable_cp_sizes==1 && design.time > 1){
                    ImGui::SetNextItemWidth(100);
                    ImGui::DragFloat("Within-cluster period c.v.", &within_cv, 0.1f, 0.0f, 10.0f, "%.1f", ImGuiSliderFlags_None);
                }
                if(ImGui::Button("Generate sizes")){
                    boost::variate_generator<boost::mt19937, boost::normal_distribution<> >
                        generator(boost::mt19937(time(0)),
                                boost::normal_distribution<>());
                    VectorXd zz(design.sequences);      
                    glmmr::randomGaussian(generator, zz); 
                    zz.array() *= cv_cp_size*mean_cp_size;    
                    zz.array() += mean_cp_size;
                    if(variable_cp_sizes==0 && design.time > 1){
                        VectorXd zz_within(design.time);
                        for (int j = 0; j < design.sequences; j++) {
                            glmmr::randomGaussian(generator, zz_within); 
                            zz_within.array() *= within_cv*(zz(j));    
                            zz_within.array() += zz(j);
                            for (int t = 0; t < design.time; t++) {
                                *(design.n(j, t)) = zz_within(t) <= 0.0 ? 1 : (int)(zz_within(t)+1);
                            }
                        }
                    } else {
                        for (int j = 0; j < design.sequences; j++) {
                            for (int t = 0; t < design.time; t++) {
                                *(design.n(j, t)) = zz(j) <= 0.0 ? 1 : (int)(zz(j)+1);
                            }
                        }
                    }
                }
            }

            ImGui::Dummy(ImVec2(0.0f, 20.0f));
            ImGui::Text("Statistical Model");
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
            

            ImGui::Dummy(ImVec2(0.0f, 20.0f));
            ImGui::Text("Correlation Structure");
            const char* cluster_covariance_items[] = { "Exchangeable", "Nested exchangeable", "Autoregressive", "Exponential", "Squared exponential" };
            const char* individual_covariance_items[] = { "Exchangeable", "Autoregressive", "Exponential", "Squared exponential" };
            static int cl_cov_item_current = 0;
            static int ind_cov_item_current = 0;
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

            ImGui::Dummy(ImVec2(0.0f, 20.0f));
            ImGui::Text("Correlations");
            ImGui::SetNextItemWidth(200);
            ImGui::DragFloat(option.heterogeneous_te ? "Control ICC" : "ICC", &model.ixx_pars[0], 0.001f, 0.0f, 1.0f, "%.3f", ImGuiSliderFlags_None); ImGui::SameLine();  HelpMarker(
                option.heterogeneous_te ? "Control group intra-class correlation coefficient." : "Intra-class correlation coefficient.");
            if(option.heterogeneous_te){
                ImGui::SetNextItemWidth(200);
                ImGui::DragFloat("Treatment ICC", &model.ixx_pars[3], 0.001f, 0.0f, 1.0f, "%.3f", ImGuiSliderFlags_None); ImGui::SameLine();  HelpMarker(
                "Treatment group intra-class correlation coefficient.");
            }
            if (model.covariance == ClusterApp::Covariance::nested_exchangeable) {
                ImGui::SetNextItemWidth(200);
                ImGui::DragFloat("CAC", &model.ixx_pars[1], 0.001f, 0.0f, 1.0f, "%.3f", ImGuiSliderFlags_None); ImGui::SameLine();  HelpMarker(
                    "Cluster autocorrelation coefficient.");
            }
            if (model.covariance == ClusterApp::Covariance::autoregressive) {
                ImGui::SetNextItemWidth(200);
                //ImGui::PushFont(unifont);ImGui::PopFont();
                ImGui::DragFloat("Autoregressive", &model.cov_pars[1], 0.01f, 0.0f, 1.0f, "%.3f", ImGuiSliderFlags_None); ImGui::SameLine(); HelpMarker(
                    "Cluster-period autoregressive parameter.");
            }
            if (model.covariance == ClusterApp::Covariance::exponential || model.covariance == ClusterApp::Covariance::squared_exponential) {
                ImGui::SetNextItemWidth(200);
                ImGui::DragFloat("Denominator", &model.cov_pars[1], 0.01f, 0.0f, +FLT_MAX, "%.3f", ImGuiSliderFlags_None); ImGui::SameLine(); HelpMarker(
                    "Cluster-period denominator parameter, e.g. exp(|t-t'|/parameter) for exponential covariance.");
            }
            if (structure_sampling == 1 || structure_sampling == 2) {
                ImGui::SetNextItemWidth(200);
                ImGui::DragFloat("IAC", &model.ixx_pars[2], 0.001f, 0.0f, 1.0f, "%.3f", ImGuiSliderFlags_None); ImGui::SameLine();  HelpMarker(
                    "Individual autocorrelation coefficient. For open cohorts, set this parameter as if it were a closed cohort.");
                if (model.ind_covariance == ClusterApp::IndividualCovariance::autoregressive && structure_sampling != 2) {
                    ImGui::SetNextItemWidth(200);
                    ImGui::DragFloat("Autoregressive (individual)", &model.cov_pars[4], 0.01f, 0.0f, 1.0f, "%.3f", ImGuiSliderFlags_None); ImGui::SameLine();  HelpMarker(
                        "Individual autoregressive parameter.");
                }
                if ((model.ind_covariance == ClusterApp::IndividualCovariance::exponential || model.ind_covariance == ClusterApp::IndividualCovariance::squared_exponential) && structure_sampling != 2) {
                    ImGui::SetNextItemWidth(200);
                    ImGui::DragFloat("Denominator (individual)", &model.cov_pars[4], 0.01f, 0.0f, +FLT_MAX, "%.3f", ImGuiSliderFlags_None); ImGui::SameLine();  HelpMarker(
                        "Individual denominator parameter, e.g. exp(|t-t'|/parameter) for exponential covariance.");
                }
            }



            ImGui::Dummy(ImVec2(0.0f, 20.0f));
            ImGui::Text("Treatment Effect");
            static float control_mean = 0.5;
            static float treatment_mean = 0.5;

            if( model.family == ClusterApp::Family::gaussian || option.two_treatments){
                ImGui::SetNextItemWidth(200);
                ImGui::DragFloat("Treatment effect parameter", &model.te_pars[0], 0.005f, -FLT_MAX, +FLT_MAX, "%.3f", ImGuiSliderFlags_None); ImGui::SameLine();  HelpMarker(
                    "Use CTRL+Click to directly input the value.");
                if (option.two_treatments) {
                    ImGui::SetNextItemWidth(200);
                    ImGui::DragFloat("Treatment 2 effect parameter", &model.te_pars[1], 0.005f, -FLT_MAX, +FLT_MAX, "%.3f", ImGuiSliderFlags_None);
                    ImGui::SetNextItemWidth(200);
                    ImGui::DragFloat("Interaction effect parameter", &model.te_pars[2], 0.005f, -FLT_MAX, +FLT_MAX, "%.3f", ImGuiSliderFlags_None);
                }
            } 

            if(model.family != ClusterApp::Family::gaussian){
            if (!option.two_treatments) {
                ImGui::SameLine(); HelpMarker("Automatically set fixed effect parameters by specifying the mean outcomes in treatment and control groups. Temporal variation is assumed to be zero. For non-linear models, a correction is applied for the random effect. Parameter values can be found above.");


                    switch (model.family) {
                    case ClusterApp::Family::gaussian: // case ClusterApp::Family::quantile:
                        ImGui::SetNextItemWidth(200);
                        ImGui::DragFloat("Control group mean", &control_mean, 0.01f, -FLT_MAX, +FLT_MAX, "%.2f", ImGuiSliderFlags_None);
                        ImGui::SetNextItemWidth(200);
                        ImGui::DragFloat("Treatment group mean", &treatment_mean, 0.01f, -FLT_MAX, +FLT_MAX, "%.2f", ImGuiSliderFlags_None);
                        break;
                    case ClusterApp::Family::binomial: case ClusterApp::Family::beta:
                        ImGui::SetNextItemWidth(200);
                        ImGui::DragFloat("Control group mean", &control_mean, 0.01f, 0.0f, 1.0f, "%.2f", ImGuiSliderFlags_None);
                        ImGui::SetNextItemWidth(200);
                        ImGui::DragFloat("Treatment group mean", &treatment_mean, 0.01f, 0.0f, 1.0f, "%.2f", ImGuiSliderFlags_None);
                        break;
                    case ClusterApp::Family::poisson: case ClusterApp::Family::gamma:
                        ImGui::SetNextItemWidth(200);
                        ImGui::DragFloat("Control group mean", &control_mean, 0.01f, 0.0f, +FLT_MAX, "%.2f", ImGuiSliderFlags_None);
                        ImGui::SetNextItemWidth(200);
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
            }
        }


        if (!option.auto_update && checker.updater.requires_update){

            ImGui::PushStyleColor(ImGuiCol_Button, colours.base1());
            ImGui::PushStyleColor(ImGuiCol_ButtonHovered, colours.base1(1));
            ImGui::PushStyleColor(ImGuiCol_ButtonActive, colours.base1(2));
            ImGui::PushStyleColor(ImGuiCol_Text, colours.base02(2));

            if (ImGui::Button("Refresh", ImVec2(80, 30))){
                checker.update();
            }
            ImGui::PopStyleColor(4);
        }

        ImGui::End();
    };

    
}

