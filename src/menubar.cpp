#include "clusterapp.h"

namespace ClusterApp {

    void RenderSelector(std::array<bool,3>& isin, ClusterApp::options& option){

        ImGui::Begin("Selector", NULL,  ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoMove);
        ImGui::Text("ClusterApp");
        ImGui::PushID(301000);
        if (ImGui::Button("Add a new model", ImVec2(150, 30))){
            for(int i = 0; i < 3; i++){
                if(!isin[i]){
                    isin[i] = true;
                    break;
                }
            }
        } 
        ImGui::PopID();
        ImGui::SetItemTooltip("The app currently supports three different design windows open at the same time, click here to open a new one.");

        //ImGui::PushID(301002);
        if (ImGui::Button("Global options", ImVec2(150, 30))){
            option.show_global_options = true;
            ImGui::OpenPopup("Global options menu");
        } 
        ImGui::SetItemTooltip("Shows some options for the app");
        //ImGui::PopID();

        if (ImGui::BeginPopupModal("Global options menu", &option.show_global_options, ImGuiWindowFlags_AlwaysAutoResize)){

            ImGui::SeparatorText("Designer view options");
            ImGui::Checkbox("Show cluster-period count (n)", &option.show_n_period); 
            ImGui::Checkbox("Show intervention status", &option.show_status_period); ImGui::SameLine(); HelpMarker("This will show a 1/0 for intervention control or the dose for a dose response design.");
            ImGui::Checkbox("Show number of clusters per sequence", &option.show_J_seq);

            ImGui::SeparatorText("Treatment effect type");
            if(!option.dose_effect && !option.time_averaged_te && !option.heterogeneous_te && !option.one_arm_clustering)ImGui::Checkbox("Two treatments", &option.two_treatments);
            if(!option.two_treatments && !option.time_averaged_te && !option.heterogeneous_te && !option.one_arm_clustering)ImGui::Checkbox("Dose effect", &option.dose_effect);
            if(!option.dose_effect && !option.two_treatments && !option.time_averaged_te&& !option.one_arm_clustering)ImGui::Checkbox("Heterogeneous", &option.heterogeneous_te);
            if(!option.dose_effect && !option.two_treatments && !option.time_averaged_te && !option.heterogeneous_te)ImGui::Checkbox("Clustering in one arm only", &option.one_arm_clustering);

            if (ImGui::Button("Close")){
                option.show_global_options = false;
                ImGui::CloseCurrentPopup();
            }
                    
            ImGui::EndPopup();
        }

        ImGui::PushID(301001);
        if (ImGui::Button("Show/hide app log", ImVec2(150, 30))){
            option.log = !option.log;
        } 
        ImGui::PopID();
        ImGui::SetItemTooltip("Shows a log of what the app is doing.");

        ImGui::Text("Version 0.7.123");
        ImGui::End();
    }


    void RenderMenuBar(ClusterApp::appModel& model, ClusterApp::options& option, bool* mopen) {
        ImGui::Begin(int_to_char("Design ",model.id+1), mopen);
            static colourPicker colours;
            const char* outcome_items[] = { "Continuous", "Binary", "Count"};
            const char* estimator_items[] = { "Mixed model", "Mixed model (t-test)", "Satterthwaite", "Kenward-Roger", "GEE independence, robust","GEE independence, robust t-test","Design effect, GLM","Design effect, non-parametric", "Box correction"};
            const char* family_items[] = { "Gaussian", "Binomial", "Poisson", "Beta", "Gamma"};
            const char* gaussian_link_items[] = { "Identity", "Log" };
            const char* binomial_link_items[] = { "Identity", "Log", "Logit", "Probit" };
            const char* poisson_link_items[] = { "Identity", "Log" };
            const char* beta_link_items[] = { "Identity", "Log", "Logit", "Probit" };
            const char* gamma_link_items[] = { "Identity", "Log", "Inverse" };
            //static int family_item_current = 0;
            //static int estimator_item_current = 0;
            //static int link_item_current = 0;
            //static bool designopen = false;
            
            ImGui::SeparatorText("Menu");
            {
                ImGui::BeginChild("ToolbarChild",ImVec2(ImGui::GetContentRegionAvail().x, 50),ImGuiChildFlags_Borders,ImGuiWindowFlags_None);
                                
                if (ImGui::Button("Choose trial design")){
                    model.designopen = true;
                    ImGui::OpenPopup("Choose design");
                }
                ImGui::SetItemTooltip("Click here to select from a range of pre-specified cluster trial designs");
                 ImGui::SameLine();

                if (ImGui::BeginPopupModal("Choose design", &model.designopen, ImGuiWindowFlags_AlwaysAutoResize))
                {
                    ImGui::Text("Trial design");
                    ImGui::TextWrapped("Select from among the default trial designs below. The trial design can be modified directly in the designer below.");
                    bool parallel_modal = true;
                    bool parallel_baseline_modal = true;
                    bool stepped_wedge_modal = true;
                    bool crossover_modal = true;
                    bool staircase_modal = true;
                    bool factorial_modal = true;

                    if (ImGui::Button("Parallel trial", ImVec2(ImGui::GetContentRegionAvail().x, 30)))
                        ImGui::OpenPopup("Parallel menu");

                    if (ImGui::Button("Parallel with baseline", ImVec2(ImGui::GetContentRegionAvail().x, 30)))
                        ImGui::OpenPopup("Parallel baseline menu");

                    if (ImGui::Button("Stepped-wedge", ImVec2(ImGui::GetContentRegionAvail().x, 30)))
                        ImGui::OpenPopup("Stepped-wedge menu");

                    if (ImGui::Button("Crossover", ImVec2(ImGui::GetContentRegionAvail().x, 30)))
                        ImGui::OpenPopup("Crossover menu");

                    if (ImGui::Button("Staircase", ImVec2(ImGui::GetContentRegionAvail().x, 30)))
                        ImGui::OpenPopup("Staircase menu");

                    if (ImGui::Button("Factorial", ImVec2(ImGui::GetContentRegionAvail().x, 30)))
                        ImGui::OpenPopup("Factorial menu");

                    if(ImGui::BeginPopupModal("Parallel menu", &parallel_modal, ImGuiWindowFlags_AlwaysAutoResize)){
                        static int parallel_t = 1;
                        static int parallel_n = 10;
                        static int parallel_J = 1;
                        ImGui::SetNextItemWidth(100);
                        ImGui::PushID(33000+model.id);
                        ImGui::InputInt("Time periods", &parallel_t, 1, 10, 0);
                        ImGui::PopID();
                        ImGui::SetNextItemWidth(100);
                        ImGui::PushID(33010+model.id);
                        ImGui::InputInt("Observations per cluster-period", &parallel_n, 1, 10, 0);
                        ImGui::PopID();
                        ImGui::SetNextItemWidth(100);
                        ImGui::PushID(33020+model.id);
                        ImGui::InputInt("Clusters per arm", &parallel_J, 1, 10, 0);
                        ImGui::PopID();
                        ImGui::PushID(33920+model.id);
                        if (ImGui::Button("Update")) {
                            parallel_t = parallel_t < 1 ? 1 : parallel_t;
                            parallel_n = parallel_n < 1 ? 1 : parallel_n;
                            parallel_J = parallel_J < 1 ? 1 : parallel_J;
                            model.designs.set_parallel(parallel_t, parallel_n, parallel_J);
                            model.model.mean_size = parallel_n;
                            model.updater.update_data();
                            ImGui::CloseCurrentPopup();
                            model.designopen = false;
                        } ImGui::SameLine();
                        ImGui::PopID();
                        ImGui::PushID(33820+model.id);
                        if (ImGui::Button("Close")) {
                            ImGui::CloseCurrentPopup();
                        }
                        ImGui::PopID();
                        ImGui::EndPopup();
                    }

                    if(ImGui::BeginPopupModal("Parallel baseline menu", &parallel_baseline_modal, ImGuiWindowFlags_AlwaysAutoResize)){
                        static int parallelb_t = 1;
                        static int parallelb_n = 10;
                        static int parallelb_J = 1;
                        ImGui::SetNextItemWidth(100);
                        ImGui::PushID(33100+model.id);
                        ImGui::InputInt("Time periods", &parallelb_t, 1, 10, 0);
                        ImGui::PopID();
                        ImGui::SetNextItemWidth(100);
                        ImGui::PushID(33110+model.id);
                        ImGui::InputInt("Observations per cluster-period", &parallelb_n, 1, 10, 0);
                        ImGui::PopID();
                        ImGui::SetNextItemWidth(100);
                        ImGui::PushID(33120+model.id);
                        ImGui::InputInt("Clusters per arm", &parallelb_J, 1, 10, 0);
                        ImGui::PopID();
                        ImGui::PushID(21100+model.id);
                        if (ImGui::Button("Update")) {
                            model.designs.set_parallel_with_baseline(parallelb_t, parallelb_n, parallelb_J);
                            model.model.mean_size = parallelb_n;
                            model.updater.update_data();
                            ImGui::CloseCurrentPopup();
                            model.designopen = false;
                        } ImGui::SameLine();
                        ImGui::PopID();
                        ImGui::PushID(21190+model.id);
                        if (ImGui::Button("Close")) {
                            ImGui::CloseCurrentPopup();
                        }
                        ImGui::PopID();
                        ImGui::EndPopup();
                    }

                    if(ImGui::BeginPopupModal("Stepped-wedge menu", &stepped_wedge_modal, ImGuiWindowFlags_AlwaysAutoResize)){
                        static int wedge_t = 6;
                        static int wedge_n = 10;
                        static int wedge_J = 1;
                        static bool implement = false;
                        ImGui::SetNextItemWidth(100);
                        ImGui::PushID(33200+model.id);
                        ImGui::InputInt("Time periods", &wedge_t, 1, 10, 0);
                        ImGui::PopID();
                        ImGui::SetNextItemWidth(100);
                        ImGui::PushID(33210+model.id);
                        ImGui::InputInt("Observations per cluster-period", &wedge_n, 1, 10, 0);
                        ImGui::PopID();
                        ImGui::SetNextItemWidth(100);
                        ImGui::PushID(33220+model.id);
                        ImGui::InputInt("Clusters per sequence", &wedge_J, 1, 10, 0);
                        ImGui::PopID();
                        ImGui::PushID(33230+model.id);
                        ImGui::Checkbox("Implementation periods", &implement);
                        ImGui::PopID();
                        ImGui::PushID(21200+model.id);
                        if (ImGui::Button("Update")) {
                            model.designs.set_stepped_wedge(wedge_t, wedge_n, wedge_J, implement);
                            model.model.mean_size = wedge_n;
                            model.updater.update_data();
                            ImGui::CloseCurrentPopup();
                            model.designopen = false;
                        } ImGui::SameLine();
                        ImGui::PopID();
                        ImGui::PushID(21290+model.id);
                        if (ImGui::Button("Close")) {
                            ImGui::CloseCurrentPopup();
                        }
                        ImGui::PopID();
                        ImGui::EndPopup();
                    }

                    if(ImGui::BeginPopupModal("Crossover menu", &crossover_modal, ImGuiWindowFlags_AlwaysAutoResize)){
                        static int cross_n;
                        static int cross_J;
                        ImGui::Text("Trial details");
                        ImGui::SetNextItemWidth(100);
                        ImGui::PushID(33300+model.id);
                        ImGui::InputInt("Observations per cluster-period", &cross_n, 1, 10, 0);
                        ImGui::PopID();
                        ImGui::SetNextItemWidth(100);
                        ImGui::PushID(33310+model.id);
                        ImGui::InputInt("Clusters per arm", &cross_J, 1, 10, 0);
                        ImGui::PopID();
                        ImGui::PushID(21300+model.id);
                        if (ImGui::Button("Update")) {
                            model.designs.set_crossover(cross_n, cross_J);
                            model.model.mean_size = cross_n;
                            model.updater.update_data();
                            ImGui::CloseCurrentPopup();
                            model.designopen = false;
                        } ImGui::SameLine();
                        ImGui::PopID();
                        ImGui::PushID(21390+model.id);
                        if (ImGui::Button("Close")) {
                            ImGui::CloseCurrentPopup();
                        }
                        ImGui::PopID();
                        ImGui::EndPopup();
                    }

                    if(ImGui::BeginPopupModal("Staircase menu", &staircase_modal, ImGuiWindowFlags_AlwaysAutoResize)){
                        static int staircase_t;
                        static int staircase_n;
                        static int staircase_J;
                        ImGui::SetNextItemWidth(100);
                        ImGui::PushID(33400+model.id);
                        ImGui::InputInt("Time periods", &staircase_t, 1, 10, 0);
                        ImGui::PopID();
                        ImGui::SetNextItemWidth(100);
                        ImGui::PushID(33410+model.id);
                        ImGui::InputInt("Observations per cluster-period", &staircase_n, 1, 10, 0);
                        ImGui::PopID();
                        ImGui::SetNextItemWidth(100);
                        ImGui::PushID(33420+model.id);
                        ImGui::InputInt("Clusters per sequence", &staircase_J, 1, 10, 0);
                        ImGui::PopID();
                        ImGui::PushID(21400+model.id);
                        if (ImGui::Button("Update")) {
                            model.designs.set_staircase(staircase_t, staircase_n, staircase_J);
                            model.model.mean_size = staircase_n;
                            model.updater.update_data();
                            ImGui::CloseCurrentPopup();
                            model.designopen = false;
                        } ImGui::SameLine();
                        ImGui::PopID();
                        ImGui::PushID(21490+model.id);
                        if (ImGui::Button("Close")) {
                            ImGui::CloseCurrentPopup();
                        }
                        ImGui::PopID();
                        ImGui::EndPopup();
                    }

                    if(ImGui::BeginPopupModal("Factorial menu", &factorial_modal, ImGuiWindowFlags_AlwaysAutoResize)){
                        if (option.two_treatments) {
                            static int factorial_t;
                            static int factorial_n;
                            static int factorial_J;
                            ImGui::SetNextItemWidth(100);
                            ImGui::PushID(33500+model.id);
                            ImGui::InputInt("Time periods", &factorial_t, 1, 10, 0);
                            ImGui::PopID();
                            ImGui::SetNextItemWidth(100);
                            ImGui::PushID(33510+model.id);
                            ImGui::InputInt("Observations per cluster-period", &factorial_n, 1, 10, 0);
                            ImGui::PopID();
                            ImGui::SetNextItemWidth(100);
                            ImGui::PushID(33520+model.id);
                            ImGui::InputInt("Clusters per arm", &factorial_J, 1, 10, 0);
                            ImGui::PopID();
                            ImGui::PushID(21500+model.id);
                            if (ImGui::Button("Update")) {
                                model.designs.set_factorial(factorial_t, factorial_n, factorial_J);
                                model.model.mean_size = factorial_n;
                                model.updater.update_data();
                                ImGui::CloseCurrentPopup();
                                model.designopen = false;
                            } ImGui::SameLine();
                            ImGui::PopID();
                        }
                        else {
                            ImGui::Text("To enable this option, set two treatments.");
                        }
                        ImGui::PushID(21590+model.id);
                        if (ImGui::Button("Close")) {
                            ImGui::CloseCurrentPopup();
                        }
                        ImGui::PopID();
                        ImGui::EndPopup();
                    }

                    if (ImGui::Button("Close")){
                        model.designopen = false;
                        ImGui::CloseCurrentPopup();
                    }
                            
                        ImGui::EndPopup();
                }


                if (ImGui::Button("Choose statistical model")){
                    model.modelopen = true;
                    ImGui::OpenPopup("Choose stat model");
                }
                ImGui::SetItemTooltip("Click here to set the statistical model, including outcome, ICC, effect size, and more");

                if (ImGui::BeginPopupModal("Choose stat model", &model.modelopen, ImGuiWindowFlags_AlwaysAutoResize))
                {
                    ImGui::Text("Statistical model");
                    static int structure_sampling = 0;
                    // ImGui::Dummy(ImVec2(0.0f, 20.0f));
                    ImGui::Text("Sampling Structure");ImGui::SameLine(); HelpMarker(
                        "In a cross-sectional design at each measurement occasion a different sample of participants is measured. In a cohort design, participants are repeatedly measured."
                    );
                    if(model.designs.time == 1){
                        structure_sampling = 0;
                        ImGui::RadioButton("Cross-section", &structure_sampling, 0); ImGui::SameLine();
                    } else {
                        ImGui::RadioButton("Cross-section", &structure_sampling, 0); ImGui::SameLine();
                        ImGui::RadioButton("Closed cohort", &structure_sampling, 1); ImGui::SameLine();
                        ImGui::RadioButton("Open cohort", &structure_sampling, 2); 
                    }
                    
                    model.model.sampling = static_cast<ClusterApp::Sampling>(structure_sampling + 1);

                    /*if(structure_sampling == 2){
                        ImGui::SetNextItemWidth(200);
                    ImGui::DragFloat("Individual replacement rate", &model.model.cov_pars[4], 0.01f, 0.0f, 1.0f, "%.3f", ImGuiSliderFlags_None); ImGui::SameLine(); HelpMarker(
                        "This is the proportion of individuals that dropout each time period and are replaced in the cohort. A value of 1 is a closed cohort and a value of 0 is a cross-sectional study.");
                    }*/



                    ImGui::Dummy(ImVec2(0.0f, 20.0f));
                    ImGui::Text("Cluster sizes");
                    static int total_t = 1;
                    static int total_n = 10;
                    static int mean_cp_size = 50;
                    static float cv_cp_size = 0.5;
                    static float within_cv = 0.0;
                    static int variable_cl_sizes = 0;
                    static int variable_cp_sizes = 0; 

                    //ImGui::Text("Variable cluster sizes"); ImGui::SameLine();
                    //ImGui::RadioButton("No", &variable_cl_sizes, 0); ImGui::SameLine();
                    //ImGui::RadioButton("Yes", &variable_cl_sizes, 1);

                    ImGui::Text("Use exact cluster-period sizes"); 
                    ImGui::SetItemTooltip("If using exact cluster-period sizes, change the cluster-period sizes by clicking the cells and rows in the designer. If not, then an approximation will be used to determine the change in power.");
                    ImGui::SameLine();
                    ImGui::RadioButton("No", &model.option.use_exact_cell_sizes, 0); ImGui::SameLine();
                    ImGui::RadioButton("Yes", &model.option.use_exact_cell_sizes, 1);

                    ImGui::SetNextItemWidth(100);
                    ImGui::InputInt("Mean cluster size", &model.model.mean_size, 1, 1000, 0);
                    if(!model.option.use_exact_cell_sizes){
                        ImGui::SetNextItemWidth(100);
                        ImGui::DragFloat("Cluster size c.v.", &model.model.cv_size, 0.1f, 0.0f, 3.0f, "%.2f", ImGuiSliderFlags_None);
                        if(model.designs.time > 1){
                            ImGui::SetNextItemWidth(100);
                            ImGui::DragFloat("Within-cluster period c.v.", &model.model.cv_size_within, 0.1f, 0.0f, 3.0f, "%.2f", ImGuiSliderFlags_None);
                        }
                    }
                    
                    ImGui::Dummy(ImVec2(0.0f, 20.0f));
                    ImGui::Text("Statistical Model");
                    ImGui::SetNextItemWidth(200);
                    ImGui::Combo("Outcome", &model.outcome_item_current, outcome_items, IM_ARRAYSIZE(outcome_items));

                    if(model.outcome_item_current == 0){
                        model.model.family = ClusterApp::Family::gaussian;
                        model.model.link = ClusterApp::Link::identity;
                    } else if(model.outcome_item_current == 1) {
                        model.model.family = ClusterApp::Family::binomial;
                        model.model.link = ClusterApp::Link::logit;
                    } else {
                        model.model.family = ClusterApp::Family::poisson;
                        model.model.link = ClusterApp::Link::log;
                    }

                    /*
                    ImGui::Combo("Family", &model.family_item_current, family_items, IM_ARRAYSIZE(family_items));

                    switch (model.family_item_current) {
                    case 0:
                    {
                        model.model.family = ClusterApp::Family::gaussian;
                        ImGui::SetNextItemWidth(200);
                        ImGui::Combo("Link", &model.link_item_current, gaussian_link_items, IM_ARRAYSIZE(gaussian_link_items));
                        switch (model.link_item_current) {
                        case 0:
                            model.model.link = ClusterApp::Link::identity;
                            break;
                        case 1:
                            model.model.link = ClusterApp::Link::log;
                            break;
                        }
                        break;
                    }
                    case 1:
                    {
                        model.model.family = ClusterApp::Family::binomial;
                        ImGui::SetNextItemWidth(200);
                        ImGui::Combo("Link", &model.link_item_current, binomial_link_items, IM_ARRAYSIZE(binomial_link_items));
                        switch (model.link_item_current) {
                        case 0:
                            model.model.link = ClusterApp::Link::identity;
                            break;
                        case 1:
                            model.model.link = ClusterApp::Link::log;
                            break;
                        case 2:
                            model.model.link = ClusterApp::Link::logit;
                            break;
                        case 3:
                            model.model.link = ClusterApp::Link::probit;
                            break;
                        }
                        break;
                    }
                    case 2:
                    {
                        model.model.family = ClusterApp::Family::poisson;
                        ImGui::SetNextItemWidth(200);
                        ImGui::Combo("Link", &model.link_item_current, poisson_link_items, IM_ARRAYSIZE(poisson_link_items));
                        switch (model.link_item_current) {
                        case 0:
                            model.model.link = ClusterApp::Link::identity;
                            break;
                        case 1:
                            model.model.link = ClusterApp::Link::log;
                            break;
                        }
                        break;
                    }
                    case 3:
                    {
                        model.model.family = ClusterApp::Family::beta;
                        ImGui::SetNextItemWidth(200);
                        ImGui::Combo("Link", &model.link_item_current, beta_link_items, IM_ARRAYSIZE(beta_link_items));
                        switch (model.link_item_current) {
                        case 0:
                            model.model.link = ClusterApp::Link::identity;
                            break;
                        case 1:
                            model.model.link = ClusterApp::Link::log;
                            break;
                        case 2:
                            model.model.link = ClusterApp::Link::logit;
                            break;
                        case 3:
                            model.model.link = ClusterApp::Link::probit;
                            break;
                        }
                        break;
                    }
                    case 4:
                    {
                        model.model.family = ClusterApp::Family::gamma;
                        ImGui::SetNextItemWidth(200);
                        ImGui::Combo("Link", &model.link_item_current, gamma_link_items, IM_ARRAYSIZE(gamma_link_items));
                        switch (model.link_item_current) {
                        case 0:
                            model.model.link = ClusterApp::Link::identity;
                            break;
                        case 1:
                            model.model.link = ClusterApp::Link::log;
                            break;
                        case 2:
                            model.model.link = ClusterApp::Link::inverse;
                            break;
                        }
                        break;
                    }
                    } */
                    

                    ImGui::Dummy(ImVec2(0.0f, 20.0f));
                    ImGui::Text("Correlation Structure");
                    const char* cluster_covariance_items[] = { "Exchangeable", "Nested exchangeable", "Autoregressive", "Exponential", "Squared exponential" };
                    const char* individual_covariance_items[] = { "Exchangeable", "Autoregressive", "Exponential", "Squared exponential" };
                    static int cl_cov_item_current = 0;
                    static int ind_cov_item_current = 0;
                    int options_size = 1;
                    if (model.designs.time > 1)options_size = IM_ARRAYSIZE(cluster_covariance_items);

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
                        ImGui::DragFloat("Individual replacement rate", &model.model.cov_pars[4], 0.01f, 0.0f, 1.0f, "%.3f", ImGuiSliderFlags_None); ImGui::SameLine(); HelpMarker(
                            "This is the proportion of individuals that dropout each time period and are replaced in the cohort. A value of 1 is a closed cohort and a value of 0 is a cross-sectional study.");
                    }

                model.model.covariance = static_cast<ClusterApp::Covariance>(cl_cov_item_current + 1);
                model.model.ind_covariance = static_cast<ClusterApp::IndividualCovariance>(ind_cov_item_current + 1);

                    ImGui::Dummy(ImVec2(0.0f, 20.0f));
                    ImGui::Text("Correlations");
                    ImGui::SetNextItemWidth(200);
                    ImGui::DragFloat(option.heterogeneous_te ? "Control ICC" : "ICC", &model.model.ixx_pars[0], 0.001f, 0.0f, 1.0f, "%.3f", ImGuiSliderFlags_None); ImGui::SameLine();  HelpMarker(
                        option.heterogeneous_te ? "Control group intra-class correlation coefficient." : "Intra-class correlation coefficient.");
                    if(option.heterogeneous_te){
                        ImGui::SetNextItemWidth(200);
                        ImGui::DragFloat("Treatment ICC", &model.model.ixx_pars[3], 0.001f, 0.0f, 1.0f, "%.3f", ImGuiSliderFlags_None); ImGui::SameLine();  HelpMarker(
                        "Treatment group intra-class correlation coefficient.");
                    }
                    if (model.model.covariance == ClusterApp::Covariance::nested_exchangeable) {
                        ImGui::SetNextItemWidth(200);
                        ImGui::DragFloat("CAC", &model.model.ixx_pars[1], 0.001f, 0.0f, 1.0f, "%.3f", ImGuiSliderFlags_None); ImGui::SameLine();  HelpMarker(
                            "Cluster autocorrelation coefficient.");
                    }
                    if (model.model.covariance == ClusterApp::Covariance::autoregressive) {
                        ImGui::SetNextItemWidth(200);
                        //ImGui::PushFont(unifont);ImGui::PopFont();
                        ImGui::DragFloat("Autoregressive", &model.model.cov_pars[1], 0.01f, 0.0f, 1.0f, "%.3f", ImGuiSliderFlags_None); ImGui::SameLine(); HelpMarker(
                            "Cluster-period autoregressive parameter.");
                    }
                    if (model.model.covariance == ClusterApp::Covariance::exponential || model.model.covariance == ClusterApp::Covariance::squared_exponential) {
                        ImGui::SetNextItemWidth(200);
                        ImGui::DragFloat("Denominator", &model.model.cov_pars[1], 0.01f, 0.0f, +FLT_MAX, "%.3f", ImGuiSliderFlags_None); ImGui::SameLine(); HelpMarker(
                            "Cluster-period denominator parameter, e.g. exp(|t-t'|/parameter) for exponential covariance.");
                    }
                    if (structure_sampling == 1 || structure_sampling == 2) {
                        ImGui::SetNextItemWidth(200);
                        ImGui::DragFloat("IAC", &model.model.ixx_pars[2], 0.001f, 0.0f, 1.0f, "%.3f", ImGuiSliderFlags_None); ImGui::SameLine();  HelpMarker(
                            "Individual autocorrelation coefficient. For open cohorts, set this parameter as if it were a closed cohort.");
                        if (model.model.ind_covariance == ClusterApp::IndividualCovariance::autoregressive && structure_sampling != 2) {
                            ImGui::SetNextItemWidth(200);
                            ImGui::DragFloat("Autoregressive (individual)", &model.model.cov_pars[4], 0.01f, 0.0f, 1.0f, "%.3f", ImGuiSliderFlags_None); ImGui::SameLine();  HelpMarker(
                                "Individual autoregressive parameter.");
                        }
                        if ((model.model.ind_covariance == ClusterApp::IndividualCovariance::exponential || model.model.ind_covariance == ClusterApp::IndividualCovariance::squared_exponential) && structure_sampling != 2) {
                            ImGui::SetNextItemWidth(200);
                            ImGui::DragFloat("Denominator (individual)", &model.model.cov_pars[4], 0.01f, 0.0f, +FLT_MAX, "%.3f", ImGuiSliderFlags_None); ImGui::SameLine();  HelpMarker(
                                "Individual denominator parameter, e.g. exp(|t-t'|/parameter) for exponential covariance.");
                        }
                    }



                    ImGui::Dummy(ImVec2(0.0f, 20.0f));
                    ImGui::Text("Treatment Effect");
                    static float control_mean = 0.5;
                    static float treatment_mean = 0.4;

                    if( model.model.family == ClusterApp::Family::gaussian || option.two_treatments){
                        ImGui::SetNextItemWidth(200);
                        ImGui::DragFloat("Treatment effect parameter", &model.model.te_pars[0], 0.005f, -FLT_MAX, +FLT_MAX, "%.3f", ImGuiSliderFlags_None); ImGui::SameLine();  HelpMarker(
                            "Use CTRL+Click to directly input the value.");
                        if (option.two_treatments) {
                            ImGui::SetNextItemWidth(200);
                            ImGui::DragFloat("Treatment 2 effect parameter", &model.model.te_pars[1], 0.005f, -FLT_MAX, +FLT_MAX, "%.3f", ImGuiSliderFlags_None);
                            ImGui::SetNextItemWidth(200);
                            ImGui::DragFloat("Interaction effect parameter", &model.model.te_pars[2], 0.005f, -FLT_MAX, +FLT_MAX, "%.3f", ImGuiSliderFlags_None);
                        }
                    } 

                    if(model.model.family != ClusterApp::Family::gaussian){
                    if (!option.two_treatments) {
                        ImGui::SameLine(); HelpMarker("Automatically set fixed effect parameters by specifying the mean outcomes in treatment and control groups. Temporal variation is assumed to be zero. For non-linear models, a correction is applied for the random effect. Parameter values can be found above.");


                            switch (model.model.family) {
                            case ClusterApp::Family::gaussian: // case ClusterApp::Family::quantile:
                                ImGui::SetNextItemWidth(200);
                                ImGui::DragFloat("Control group mean", &control_mean, 0.01f, -FLT_MAX, +FLT_MAX, "%.3f", ImGuiSliderFlags_None);
                                ImGui::SetNextItemWidth(200);
                                ImGui::DragFloat("Treatment group mean", &treatment_mean, 0.01f, -FLT_MAX, +FLT_MAX, "%.3f", ImGuiSliderFlags_None);
                                break;
                            case ClusterApp::Family::binomial: case ClusterApp::Family::beta:
                                ImGui::SetNextItemWidth(200);
                                ImGui::DragFloat("Control group mean", &control_mean, 0.01f, 0.0f, 1.0f, "%.3f", ImGuiSliderFlags_None);
                                ImGui::SetNextItemWidth(200);
                                ImGui::DragFloat("Treatment group mean", &treatment_mean, 0.01f, 0.0f, 1.0f, "%.3f", ImGuiSliderFlags_None);
                                break;
                            case ClusterApp::Family::poisson: case ClusterApp::Family::gamma:
                                ImGui::SetNextItemWidth(200);
                                ImGui::DragFloat("Control group mean", &control_mean, 0.01f, 0.0f, +FLT_MAX, "%.3f", ImGuiSliderFlags_None);
                                ImGui::SetNextItemWidth(200);
                                ImGui::DragFloat("Treatment group mean", &treatment_mean, 0.01f, 0.0f, +FLT_MAX, "%.3f", ImGuiSliderFlags_None);
                                break;
                            }


                            //if (ImGui::Button("Update")) {
                            //    
                            //}
                        }
                    }
                    if(ImGui::Button("Update")){
                        if(model.model.family != ClusterApp::Family::gaussian){
                            float re_adj = model.model.cov_pars[0];
                            if (model.model.covariance == ClusterApp::Covariance::nested_exchangeable)re_adj += model.model.cov_pars[1];
                            if (model.model.sampling == ClusterApp::Sampling::cohort) re_adj += model.model.cov_pars[3] / model.designs.mean_n();
                            re_adj *= 0.5;
                            switch (model.model.link) {
                            case ClusterApp::Link::identity:
                            {
                                if (model.model.include_intercept == 1) {
                                    model.model.beta_pars[0] = control_mean;
                                }
                                else {
                                    for (int l = 0; l < model.model.beta_pars.size(); l++) {
                                        model.model.beta_pars[l] = control_mean;
                                    }
                                }
                                model.model.te_pars[0] = treatment_mean - control_mean;
                                break;
                            }
                            case ClusterApp::Link::log:
                            {
                                if (model.model.include_intercept == 1) {
                                    model.model.beta_pars[0] = log(control_mean) - re_adj;
                                }
                                else {
                                    for (int l = 0; l < model.model.beta_pars.size(); l++) {
                                        model.model.beta_pars[l] = log(control_mean) - re_adj;
                                    }
                                }
                                model.model.te_pars[0] = log(treatment_mean) - log(control_mean);
                                break;
                            }
                            case ClusterApp::Link::logit:
                            {
                                if (model.model.include_intercept == 1) {
                                    model.model.beta_pars[0] = log(control_mean / (1 - control_mean)) - re_adj;
                                }
                                else {
                                    for (int l = 0; l < model.model.beta_pars.size(); l++) {
                                        model.model.beta_pars[l] = log(control_mean / (1 - control_mean)) - re_adj;
                                    }
                                }
                                model.model.te_pars[0] = log(treatment_mean / (1 - treatment_mean)) - log(control_mean / (1 - control_mean));
                                break;
                            }
                            case ClusterApp::Link::probit:
                            {
                                boost::math::normal norm = boost::math::normal(0.0, 1.0);
                                if (model.model.include_intercept == 1) {
                                    model.model.beta_pars[0] = boost::math::quantile(norm, control_mean) - re_adj;
                                }
                                else {
                                    for (int l = 0; l < model.model.beta_pars.size(); l++) {
                                        model.model.beta_pars[l] = boost::math::quantile(norm, control_mean) - re_adj;
                                    }
                                }
                                model.model.te_pars[0] = boost::math::quantile(norm, treatment_mean) - model.model.beta_pars[0];
                                break;
                            }
                            case ClusterApp::Link::inverse:
                            {
                                if (model.model.include_intercept == 1) {
                                    model.model.beta_pars[0] = 1/control_mean;
                                }
                                else {
                                    for (int l = 0; l < model.model.beta_pars.size(); l++) {
                                        model.model.beta_pars[l] = 1/control_mean;
                                    }
                                }
                                model.model.te_pars[0] = 1/(treatment_mean - control_mean);
                                break;
                            }

                            }
                        }

                        model.checker.update();
                        model.modelopen = false;
                        ImGui::CloseCurrentPopup();
                    }
                    ImGui::SameLine();
                    if (ImGui::Button("Close")){
                        model.modelopen = false;
                        ImGui::CloseCurrentPopup();
                    }
                            
                    ImGui::EndPopup();
                }


                //ImGui::Text("Outcome: "); ImGui::SameLine();
                //ImGui::Text(outcome_items[model.outcome_item_current]); 

                ImGui::EndChild();
            }

            
        // trial design graphical
        ImGui::SeparatorText("Designer");
        {
            enum Mode
            {
                Mode_Copy,
                Mode_Move,
                Mode_Swap
            };
            static int mode = 2;

            if (ImGui::TreeNode("How to use the designer")){
                ImGui::TextWrapped("Design the trial below. Rows are sequences, columns are time periods. Click + to add new sequences or time periods. Select cells to edit their details or you can drag and drop them to new positions. Select row or column headers to change the numbers of clusters."); ImGui::SameLine(); HelpMarker(
                "You can change what the cell buttons show with the buttons below. Red and blue indicate intervention and control status, respectively. Where there are two treatments, yellow is used for treatment 2, and yellow-red for both treatments");
                
                if (ImGui::TreeNode("Drag and drop mode")) {            
                    ImGui::Text("Drag and drop mode: "); ImGui::SameLine();
                    if (ImGui::RadioButton("Copy", mode == Mode_Copy)) { mode = Mode_Copy; } ImGui::SameLine();
                    if (ImGui::RadioButton("Move", mode == Mode_Move)) { mode = Mode_Move; } ImGui::SameLine();
                    if (ImGui::RadioButton("Swap", mode == Mode_Swap)) { mode = Mode_Swap; }
                    ImGui::TreePop();
                }
                
                ImGui::TreePop();
            }

            static int max_large_dim = 50;
            
            ImGui::SetNextItemWidth(100);
                ImGui::DragInt("Cell size", &max_large_dim, 5, 20, 150, "%d px", ImGuiSliderFlags_AlwaysClamp); ImGui::SameLine();
                ImGui::SetItemTooltip("Change the size of the cells in the diagram below by dragging this");


            if (ImGui::Button("Optimal design weights")){
                model.optimalopen = true;
                ImGui::OpenPopup("Optimum design");
            } ImGui::SameLine();
            ImGui::SetItemTooltip("Shows the cluster-period weights that would be most efficient with the current sample size");

            if (ImGui::BeginPopupModal("Optimum design", &model.optimalopen, ImGuiWindowFlags_AlwaysAutoResize))
            {
                ImGui::TextWrapped("The design below shows the optimum weights per cluster in the design selected in the main window. Where sequences contain more than one cluster, these have been split out "); ImGui::SameLine(); HelpMarker(
                    "The colour of each cell indicates the weight. Click on the cell to see the weight and the rounded total number of observation. Rounding uses Hamilton's method. The design can be applied to the designer using the button below.");
        
                if (option.two_treatments) {
                    ImGui::TextWrapped("For two treatments, the design minimises the variance of a weighted combination of the two treatment effect parameters and the interaction. You can change the weights in the options menu of this window.");
                }
        
                static ClusterApp::colourPicker colors;
        
                if (model.updater.optimum_data.size() != model.updater.designs.total_clusters()) {
                    model.updater.update_data();
                    model.updater.update_optimum();
                }
        
                int T = model.updater.optimum_data[0].size();
                int N = model.updater.optimum_data.size();
                int dim = 35;
                ImGuiStyle& style = ImGui::GetStyle();
                int horizontal_align = dim;
        
                ImGui::Dummy(ImVec2(dim * 0.4, dim * 0.4)); ImGui::SameLine(horizontal_align);
                for (int t = 0; t < T; t++) {
                    ImGui::Text("%i", t + 1);
                    if (t < T - 1)ImGui::SameLine(horizontal_align + (t + 1) * (dim + style.ItemSpacing[0]));
                }
        
                for (int i = 0; i < N; i++) {
                    int seq = model.updater.designs.seq_by_cluster(i);
                    ImGui::Text("%i", i + 1); ImGui::SameLine(horizontal_align);
                    for (int t = 0; t < T; t++) {
                        float w = model.updater.optimum_data[i][t];
                        ImVec4 newcol;
                        if (model.option.light_mode) {
                            newcol = colors.base3();
                        }
                        else {
                            newcol = colors.base03();
                        }
                        bool both_two = model.option.two_treatments && *model.designs.intervention(seq, t) && *model.designs.intervention_2(seq, t);
                        if (w > 1e-3) {
                            if (*model.designs.intervention(seq, t) && !(model.option.two_treatments && *model.designs.intervention_2(seq, t))) {
                                newcol = colors.red();
                                newcol.w *= 2 * w;
                            }
                            else if (model.option.two_treatments && *model.designs.intervention_2(seq, t)) {
                                newcol = colors.yellow();
                                newcol.w *= 2 * w;
                            }
                            else {
                                newcol = colors.blue();
                                newcol.w *= 2 * w;
                            }
                        }
        
                        ImGui::PushStyleColor(ImGuiCol_Button, newcol);
                        newcol.y += 0.1;
                        ImGui::PushStyleColor(ImGuiCol_ButtonHovered, newcol);
                        newcol.y += 0.1;
                        ImGui::PushStyleColor(ImGuiCol_ButtonActive, newcol);
                        int newid = t + i * T;
                        ImGui::PushID(10 * model.id * model.designs.sequences * model.designs.time + newid);
                        if (!both_two) {
                            ImGui::Button("", ImVec2(dim, dim));
                            if (ImGui::BeginPopupContextItem(NULL, ImGuiPopupFlags_MouseButtonLeft))
                            {
                                ImGui::Text("Weight: "); ImGui::SameLine();
                                ImGui::Text("%.3f", w);
                                ImGui::Text("n: "); ImGui::SameLine();
                                ImGui::Text("%d", model.updater.optimum_n[i][t]);
                                ImGui::EndPopup();
                            }
                        }
                        else {
                            ImVec4 yel = colors.yellow();
                            ImVec4 red = colors.yellow();
                            yel.w *= 2 * w;
                            red.w *= 2 * w;
                            ClusterApp::MultiColorButton("", yel, yel, red, red, 0, 0,
                                ImVec2(dim, dim),
                                yel, yel, red, red,
                                yel, yel, red, red);
                            if (ImGui::BeginPopupContextItem(NULL, ImGuiPopupFlags_MouseButtonLeft))
                            {
                                ImGui::Text("Weight: "); ImGui::SameLine();
                                ImGui::Text("%.3f", w);
                                ImGui::Text("n: "); ImGui::SameLine();
                                ImGui::Text("%d", model.updater.optimum_n[i][t]);
                                ImGui::EndPopup();
                            }
                        }
        
                        ImGui::PopStyleColor(3);
                        ImGui::PopID();
                        if (t < model.updater.optimum_data[0].size() - 1)ImGui::SameLine();
                    }
                }
                ImGui::Dummy(ImVec2(dim, dim));
                ImGui::Dummy(ImVec2(dim * 0.5, dim * 0.5)); ImGui::SameLine();
                ImGui::PushID(78700 + model.id);
                if (ImGui::SmallButton("Apply Design")) {
                    if(model.option.use_exact_cell_sizes){
                        model.designs.apply_design(model.updater.optimum_n);
                    } else {
                        ImGui::OpenPopup("Set exact sizes");
                    }
                    
                }

                if (ImGui::BeginPopupModal("Set exact sizes", NULL, ImGuiWindowFlags_AlwaysAutoResize))
                {
                    ImGui::Text("The designer will use exact cell sizes and not coefficients of variation.\nTo change this setting, go to the statistical model.");
                    ImGui::Separator();
        
                    if (ImGui::Button("OK", ImVec2(120, 0))) { 
                        model.designs.apply_design(model.updater.optimum_n);
                        model.option.use_exact_cell_sizes = 1;
                        model.optimalopen = false;
                        ImGui::CloseCurrentPopup(); 
                    }
                    ImGui::SetItemDefaultFocus();
                    ImGui::SameLine();
                    if (ImGui::Button("Cancel", ImVec2(120, 0))) { ImGui::CloseCurrentPopup(); }
                    ImGui::EndPopup();
                }

                ImGui::PopID();

                if (ImGui::Button("Close")){
                    model.optimalopen = false;
                    ImGui::CloseCurrentPopup();
                }
                        
                ImGui::EndPopup();
            }

            if (ImGui::Button("Sample size summary")){
                model.samplesizeopen = true;
                ImGui::OpenPopup("Sample size total");
            } 
            ImGui::SetItemTooltip("Shows a table with the total numbers of observations in the design");

            if (ImGui::BeginPopupModal("Sample size total", &model.samplesizeopen, ImGuiWindowFlags_AlwaysAutoResize)){
                static ImGuiTableFlags flags = ImGuiTableFlags_Borders | ImGuiTableFlags_Resizable | ImGuiTableFlags_Reorderable | ImGuiTableFlags_Hideable;// | ImGuiTableFlags_NoHostExtendX;

                ImGui::Text("SAMPLE SIZE");
                if (ImGui::BeginTable("summary", 4, flags))
                {
                    ImGui::TableSetupColumn("Sequence");
                    ImGui::TableSetupColumn("N. clusters");
                    ImGui::TableSetupColumn("N. observations");
                    ImGui::TableSetupColumn("Total observations");
                    ImGui::TableHeadersRow();
                    ImGui::TableNextColumn();
                    int total_total_n = 0;
                    for (int i = 0; i < model.designs.sequences; i++) {
                        ImGui::Text("%i",i + 1);
                        ImGui::TableNextColumn();
                        ImGui::Text("%i", *model.designs.n_clusters(i));
                        ImGui::TableNextColumn();
                        int total_n = 0;
                        if (ImGui::BeginTable("summaryn", model.designs.time, ImGuiTableFlags_Borders | ImGuiTableFlags_Resizable | ImGuiTableFlags_Reorderable | ImGuiTableFlags_Hideable)) {
                            for (int t = 0; t < model.designs.time; t++) {
                                ImGui::TableSetupColumn(int_to_char(t + 1));
                            }
                            ImGui::TableHeadersRow();
                            ImGui::TableNextColumn();
                            for (int j = 0; j < *model.designs.n_clusters(i); j++) {
                                for (int t = 0; t < model.designs.time; t++) {
                                    int nper = 0;
                                    if (*model.designs.active(i, t))nper = *model.designs.n(i, t);
                                    ImGui::Text("%i", nper);
                                    total_n += nper;
                                    ImGui::TableNextColumn();
                                }
                            }
                            ImGui::EndTable();
                        }
                        ImGui::TableNextColumn();
                        ImGui::Text("%i", total_n);
                        total_total_n += total_n;
                        ImGui::TableNextColumn();
                    }
                    ImGui::Text("TOTAL");
                    ImGui::TableNextColumn();
                    ImGui::Text("%i", model.designs.total_clusters());
                    ImGui::TableNextColumn(); ImGui::TableNextColumn();
                    ImGui::Text("%i", total_total_n);
                    ImGui::EndTable();
                }

                if (ImGui::Button("Close")){
                    model.samplesizeopen = false;
                    ImGui::CloseCurrentPopup();
                }
                        
                ImGui::EndPopup();
            }

            ImGuiStyle& style = ImGui::GetStyle();
            static colourPicker colours;
            static int horiztonal_align = 40;
            const short   s16_zero = 0, s16_one = 1;
            static bool remove_seq = false;
            static int which_remove = 0;
            static int set_all_n = 0;
            // static int max_large_dim = 50;
            static int large_dim;
            static int small_dim = 20;
            float window_visible_x2 = ImGui::GetContentRegionAvail().x; //ImGui::GetWindowPos().x + 
            float large_test_size = (window_visible_x2 - (3 * small_dim) - ((model.designs.time + 5) * style.ItemSpacing.x)) * (1.0f / model.designs.time);
            if (large_test_size > 20) {
                if (large_test_size < max_large_dim) {
                    large_dim = large_test_size;
                }
                else {
                    large_dim = max_large_dim;
                }
            }
            else {
                large_dim = small_dim;
            }
            

            ImGui::BeginChild("DesignerChild",ImVec2(ImGui::GetContentRegionAvail().x, model.designs.sequences*large_dim + 5*small_dim),ImGuiChildFlags_Borders,ImGuiWindowFlags_None);
            
            // Variables and other setup for the designer
        
        // View options for the designer
        

        ImGui::Spacing();

        // The designer 
        ImGui::Dummy(ImVec2(small_dim, small_dim)); ImGui::SameLine();
        ImGui::Dummy(ImVec2(small_dim, small_dim)); ImGui::SameLine();
        if (option.show_J_seq) {
            ImGui::Dummy(ImVec2(small_dim * 1.5, small_dim)); ImGui::SameLine();
        }
        ImGui::PushStyleColor(ImGuiCol_Button, colours.base1());
        ImGui::PushStyleColor(ImGuiCol_ButtonHovered, colours.base1(1));
        ImGui::PushStyleColor(ImGuiCol_ButtonActive, colours.base1(2));
        ImGui::PushStyleColor(ImGuiCol_Text, colours.base02(2));

        // First iterate over the time periods to add the + buttons to extend the columns
        for (int t = 0; t <= model.designs.time; t++) {
            ImGui::PushID(model.designs.time * model.designs.sequences + t);
            if (ImGui::Button("+", ImVec2(small_dim, small_dim))) {
                if (t == model.designs.time) {
                    model.designs.add_period();
                }
                else {
                    model.designs.add_period(t);
                }
            };
            ImGui::SetItemTooltip("Add new time period");
            ImGui::PopID();
            ImGui::SameLine();
            ImGui::Dummy(ImVec2(large_dim - small_dim - style.ItemSpacing[0], small_dim));
            if (t < (model.designs.time))ImGui::SameLine();
        }

        ImGui::Dummy(ImVec2(small_dim, small_dim)); ImGui::SameLine();
        ImGui::Dummy(ImVec2(small_dim, small_dim)); ImGui::SameLine();
        if (option.show_J_seq) {
            ImGui::Dummy(ImVec2(small_dim * 1.5, small_dim)); ImGui::SameLine();
        }

        // now iterate again to create the column heads
        for (int t = 0; t < model.designs.time; t++) {
            ImGui::PushID(model.designs.time * model.designs.sequences + model.designs.time + 2 + model.designs.sequences + t);
            ImGui::Button(int_to_char(t + 1), ImVec2(large_dim, small_dim));
            ImGui::SetItemTooltip("Edit time period");
            if (ImGui::BeginPopupContextItem(NULL, ImGuiPopupFlags_MouseButtonLeft))
            {
                ImGui::Text("OPTIONS");
                if (ImGui::SmallButton("Activate all")) {
                    for (int s = 0; s < model.designs.sequences; s++)*model.designs.active(s, t) = true;
                }
                if (ImGui::SmallButton("De-activate all")) {
                    for (int s = 0; s < model.designs.sequences; s++)*model.designs.active(s, t) = false;
                }
                if (ImGui::SmallButton("Set all intervention")) {
                    for (int s = 0; s < model.designs.sequences; s++) {
                        *model.designs.intervention(s, t) = true;
                        *model.designs.intervention_2(s, t) = false;
                    }
                }
                if (option.dose_effect) {
                    static float dose_seq = 1.0;
                    ImGui::SetNextItemWidth(100);
                    ImGui::DragFloat("Dose", &dose_seq, 0.1f, 0.0, +FLT_MAX, "%.1f", ImGuiSliderFlags_None); ImGui::SameLine();
                    if (ImGui::SmallButton("Set dose")) {
                        for (int s = 0; s < model.designs.sequences; s++) {
                            *model.designs.dose(s, t) = dose_seq;
                        }
                    }
                }
                if (option.two_treatments) {
                    if (ImGui::SmallButton("Set all intervention 2")) {
                        for (int s = 0; s < model.designs.sequences; s++) {
                            *model.designs.intervention(s, t) = false;
                            *model.designs.intervention_2(s, t) = true;
                        }
                    }
                    if (ImGui::SmallButton("Set all intervention 1 & 2")) {
                        for (int s = 0; s < model.designs.sequences; s++) {
                            *model.designs.intervention(s, t) = true;
                            *model.designs.intervention_2(s, t) = true;
                        }
                    }
                }
                if (ImGui::SmallButton("Set all control")) {
                    for (int s = 0; s < model.designs.sequences; s++)*model.designs.intervention(s, t) = false;
                }
                if (model.designs.time > 1) {
                    if (ImGui::SmallButton("Remove column")) {
                        model.designs.remove_period(t);
                    }
                }
                ImGui::EndPopup();
            }
            ImGui::PopID();
            if (t < (model.designs.time - 1))ImGui::SameLine();
        }
        ImGui::PopStyleColor(4);

        // iterate over the sequences to produce the row heads and main buttons

        for (int n = 0; n < model.designs.sequences; n++)
        {
            if (option.show_J_seq) {
                std::string label_j = std::to_string(*model.designs.n_clusters(n));
                char* char_array_j = new char[label_j.length() + 1];
                strcpy(char_array_j, label_j.c_str());
                ImGui::PushID(20000 + n);
                ImGui::Button(char_array_j, ImVec2(small_dim * 1.5, small_dim));
                ImGui::SetItemTooltip("Number of clusters per sequence");
                if (ImGui::BeginPopupContextItem(NULL, ImGuiPopupFlags_MouseButtonLeft)) {
                    static int n_clusters = 10;
                    ImGui::Text("Number of clusters");
                    ImGui::SetNextItemWidth(100);
                    ImGui::InputScalar("N", ImGuiDataType_S16, &n_clusters, &s16_one, NULL, "%d");
                    if (ImGui::Button("Set")) {
                        *(model.designs.n_clusters(n)) = n_clusters;
                    }
                    ImGui::EndPopup();
                }
                ImGui::PopID();
                ImGui::SameLine();
            }

            ImGui::PushStyleColor(ImGuiCol_Button, colours.base1());
            ImGui::PushStyleColor(ImGuiCol_ButtonHovered, colours.base1(1));
            ImGui::PushStyleColor(ImGuiCol_ButtonActive, colours.base1(2));
            ImGui::PushStyleColor(ImGuiCol_Text, colours.base02(2));

            ImGui::PushID(model.designs.time * model.designs.sequences + model.designs.time + 1 + n);
            if (ImGui::Button("+", ImVec2(small_dim, small_dim))) {
                model.designs.add_sequence(n);
            }
            ImGui::SetItemTooltip("Add a new sequence");
            ImGui::SameLine();
            ImGui::PopID();
            // Row heads
            ImGui::PushID(model.designs.time * (model.designs.sequences + 1) + model.designs.time + model.designs.sequences + 2 + n);
            ImGui::Button(int_to_char(n + 1), ImVec2(small_dim, large_dim)); 
            ImGui::SetItemTooltip("Click to edit sequence");
            ImGui::SameLine();
            if (ImGui::BeginPopupContextItem(NULL, ImGuiPopupFlags_MouseButtonLeft))
            {
                static int seq_n = 10;
                static float dose = 1.0;
                ImGui::Text("Number of clusters");
                ImGui::SetNextItemWidth(100);
                ImGui::InputScalar("N", ImGuiDataType_S16, &seq_n, &s16_one, NULL, "%d"); ImGui::SameLine();
                if (ImGui::Button("Set")) {
                    *(model.designs.n_clusters(n)) = seq_n;
                }
                ImGui::Text("OPTIONS");
                if (ImGui::SmallButton("Activate all")) {
                    for (int s = 0; s < model.designs.time; s++)*model.designs.active(n, s) = true;
                }
                if (ImGui::SmallButton("De-activate all")) {
                    for (int s = 0; s < model.designs.time; s++)*model.designs.active(n, s) = false;
                }
                if (ImGui::SmallButton("Set all intervention")) {
                    for (int s = 0; s < model.designs.time; s++) {
                        *model.designs.intervention(n, s) = true;
                        *model.designs.intervention_2(n, s) = false;
                    }
                }
                if (option.dose_effect) {
                    ImGui::SetNextItemWidth(100);
                    ImGui::DragFloat("Dose", &dose, 0.1f, 0.0, +FLT_MAX, "%.1f", ImGuiSliderFlags_None); ImGui::SameLine();
                    if (ImGui::SmallButton("Set dose")) {
                        for (int s = 0; s < model.designs.time; s++) {
                            *model.designs.dose(n, s) = dose;
                        }
                    }
                }
                if (option.two_treatments) {
                    if (ImGui::SmallButton("Set all intervention 2")) {
                        for (int s = 0; s < model.designs.time; s++) {
                            *model.designs.intervention(n, s) = false;
                            *model.designs.intervention_2(n, s) = true;
                        }
                    }
                    if (ImGui::SmallButton("Set all intervention 1 & 2")) {
                        for (int s = 0; s < model.designs.time; s++) {
                            *model.designs.intervention(n, s) = true;
                            *model.designs.intervention_2(n, s) = true;
                        }
                    }
                }
                if (ImGui::SmallButton("Set all control")) {
                    for (int s = 0; s < model.designs.time; s++)*model.designs.intervention(n, s) = false;
                }
                if (model.designs.sequences > 1) {
                    if (ImGui::SmallButton("Remove sequence")) {
                        remove_seq = true;
                        which_remove = n;
                    }
                }
                if(model.option.use_exact_cell_sizes){
                    ImGui::SetNextItemWidth(100);
                    ImGui::InputScalar("n", ImGuiDataType_S16, &set_all_n, NULL, NULL, "%d"); ImGui::SameLine();
                    if (ImGui::SmallButton("Set all n")) {
                        for (int s = 0; s < model.designs.time; s++)*model.designs.n(n, s) = set_all_n;
                    }
                }                

                ImGui::EndPopup();
            }
            ImGui::PopID();
            ImGui::PopStyleColor(4);
            
            // cluster-period buttons
            for (int t = 0; t < model.designs.time; t++) {
                int id = t + n * model.designs.time;
                ImGui::PushID(id);
                std::string label = "";
                if (*model.designs.active(n, t)) {
                    if (option.show_n_period || option.use_exact_cell_sizes) {
                        label += "n=";
                        label += std::to_string(*model.designs.n(n, t));
                        if (option.show_status_period) label += "\n"; //replace with | for single line
                    }
                    if (option.show_status_period) {                        
                        label += option.dose_effect ? std::format("{:.2f}", (*model.designs.intervention(n, t)) * (*model.designs.dose(n, t))) : std::to_string(*model.designs.intervention(n, t));
                        if (option.two_treatments)label += "/" + std::to_string(*model.designs.intervention_2(n, t));
                    }

                    if (*model.designs.intervention(n, t) && !*model.designs.intervention_2(n, t)) {
                        ImGui::PushStyleColor(ImGuiCol_Button, colours.red());
                        ImGui::PushStyleColor(ImGuiCol_ButtonHovered, colours.red(1));
                        ImGui::PushStyleColor(ImGuiCol_ButtonActive, colours.red(2));
                    }
                    else if (!*model.designs.intervention(n, t) && *model.designs.intervention_2(n, t)) {
                        ImGui::PushStyleColor(ImGuiCol_Button, colours.yellow());
                        ImGui::PushStyleColor(ImGuiCol_ButtonHovered, colours.yellow(1));
                        ImGui::PushStyleColor(ImGuiCol_ButtonActive, colours.yellow(2));
                    }
                    else {
                        ImGui::PushStyleColor(ImGuiCol_Button, colours.blue());
                        ImGui::PushStyleColor(ImGuiCol_ButtonHovered, colours.blue(1));
                        ImGui::PushStyleColor(ImGuiCol_ButtonActive, colours.blue(2));
                    }
                }
                else {
                    // inactive cluster period
                    if (option.light_mode) {
                        ImGui::PushStyleColor(ImGuiCol_Button, colours.base2());
                        ImGui::PushStyleColor(ImGuiCol_ButtonHovered, colours.base2(1));
                        ImGui::PushStyleColor(ImGuiCol_ButtonActive, colours.base2(2));
                    }
                    else {
                        ImGui::PushStyleColor(ImGuiCol_Button, colours.base02());
                        ImGui::PushStyleColor(ImGuiCol_ButtonHovered, colours.base02(1));
                        ImGui::PushStyleColor(ImGuiCol_ButtonActive, colours.base02(2));
                    }

                }
                char* char_array = new char[label.length() + 1];
                strcpy(char_array, label.c_str());
                if (!option.two_treatments || (option.two_treatments && !(*model.designs.intervention(n, t) && *model.designs.intervention_2(n, t)))) {
                    ImGui::Button(char_array, ImVec2(large_dim, large_dim));
                    ImGui::SetItemTooltip("Cluster-period, blue = control, red = intervention");
                    // popup context menu for cluster-period cells
                    if (ImGui::BeginPopupContextItem(NULL, ImGuiPopupFlags_MouseButtonLeft))
                    {
                        static int cell_n = 10;
                        static float dose_n = 1.0;
                        ImGui::Checkbox("Active", model.designs.active(n, t));
                        ImGui::Checkbox("Intervention", model.designs.intervention(n, t));
                        if (option.two_treatments)ImGui::Checkbox("Intervention 2", model.designs.intervention_2(n, t));
                        if (option.dose_effect && *model.designs.intervention(n, t)) {
                            ImGui::SetNextItemWidth(100);
                            ImGui::DragFloat("Dose", &dose_n, 0.1f, 0.0, +FLT_MAX, "%.1f", ImGuiSliderFlags_None);
                        }
                        if(model.option.use_exact_cell_sizes){
                            ImGui::SetNextItemWidth(100);
                            ImGui::InputScalar("N", ImGuiDataType_S16, &cell_n, &s16_one, NULL, "%d"); ImGui::SameLine();
                            if (ImGui::Button("Set")) {
                                *(model.designs.n(n, t)) = cell_n;
                                if (option.dose_effect && *model.designs.intervention(n, t)) *model.designs.dose(n, t) = dose_n;
                            }
                        }                        
                        ImGui::EndPopup();
                    }
                }
                else {
                    ClusterApp::MultiColorButton(char_array, colours.yellow(), colours.yellow(), colours.red(), colours.red(), 0, 0, ImVec2(large_dim, large_dim),
                        colours.yellow(1), colours.yellow(1), colours.red(1), colours.red(1),
                        colours.yellow(2), colours.yellow(2), colours.red(2), colours.red(2));
                    if (ImGui::BeginPopupContextItem(NULL, ImGuiPopupFlags_MouseButtonLeft))
                    {
                        ImGui::Checkbox("Active", model.designs.active(n, t));
                        ImGui::Checkbox("Intervention", model.designs.intervention(n, t));
                        ImGui::Checkbox("Intervention 2", model.designs.intervention_2(n, t));
                        if(model.option.use_exact_cell_sizes){
                            ImGui::SetNextItemWidth(100);
                            ImGui::InputScalar("N", ImGuiDataType_S16, model.designs.n(n, t), &s16_one, NULL, "%d");
                        }                        
                        ImGui::EndPopup();
                    }
                }
                ImGui::PopStyleColor(3);

                // Our buttons are both drag sources and drag targets here!
                if (ImGui::BeginDragDropSource(ImGuiDragDropFlags_None))
                {
                    ImGui::SetDragDropPayload("DND_DEMO_CELL", &id, sizeof(int));
                    if (mode == Mode_Copy) { ImGui::Text("Copy"); }
                    if (mode == Mode_Move) { ImGui::Text("Move"); }
                    if (mode == Mode_Swap) { ImGui::Text("Swap"); }
                    ImGui::EndDragDropSource();
                }
                if (ImGui::BeginDragDropTarget())
                {
                    if (const ImGuiPayload* payload = ImGui::AcceptDragDropPayload("DND_DEMO_CELL"))
                    {
                        IM_ASSERT(payload->DataSize == sizeof(int));
                        int payload_id = *(const int*)payload->Data;
                        if (mode == Mode_Copy)
                        {
                            model.designs.copy_cells(payload_id, id);
                        }
                        if (mode == Mode_Move)
                        {
                            model.designs.move_cells(payload_id, id);
                        }
                        if (mode == Mode_Swap)
                        {
                            model.designs.swap_cells(payload_id, id);
                        }
                    }
                    ImGui::EndDragDropTarget();
                }

                ImGui::PopID();
                if (t < (model.designs.time - 1))ImGui::SameLine();
            }
        }

        ImGui::PushStyleColor(ImGuiCol_Button, colours.base1());
        ImGui::PushStyleColor(ImGuiCol_ButtonHovered, colours.base1(1));
        ImGui::PushStyleColor(ImGuiCol_ButtonActive, colours.base1(2));
        ImGui::PushStyleColor(ImGuiCol_Text, colours.base02(2));
        if (option.show_J_seq) {
            ImGui::Dummy(ImVec2(small_dim * 1.5, small_dim)); ImGui::SameLine();
        }

        ImGui::PushID(model.designs.time * model.designs.sequences + model.designs.time + 1 + model.designs.sequences);
        if (ImGui::Button("+", ImVec2(small_dim, small_dim))) {
            model.designs.add_sequence();
        }
        ImGui::SetItemTooltip("Add a new sequence");
        ImGui::PopID();
        ImGui::PopStyleColor(4);

        if (remove_seq) {
            model.designs.remove_sequence(which_remove);
            remove_seq = false;
        }


            ImGui::EndChild();
        }

        ImGui::SeparatorText("Estimator");
        {
            ImGui::BeginChild("EstimatorChild",ImVec2(ImGui::GetContentRegionAvail().x, 100),ImGuiChildFlags_Borders,ImGuiWindowFlags_None);
            ImGui::SetNextItemWidth(250);
            ImGui::Combo("Estimator", &model.estimator_item_current, estimator_items, IM_ARRAYSIZE(estimator_items));
            ImGui::SetItemTooltip("Choose the statistical estimator of the treatment effect.");
            model.model.powertype = static_cast<ClusterApp::PowerType>(model.estimator_item_current);

            ImGui::SetNextItemWidth(250);
            ImGui::PushID(77760+model.id);
            ImGui::DragFloat("Target Power (%)", &model.model.target_power, 1.0f, 0.0f, 100.0f, "%.0f", ImGuiSliderFlags_None);  
            ImGui::PopID();
            ImGui::SetItemTooltip("Set target power for minimum detectable effect sizes and sample size searches");

            ImGui::SetNextItemWidth(250);
            ImGui::PushID(77790+model.id);
            ImGui::DragFloat("Alpha", &model.model.alpha, 0.01f, 0.0f, 1.0f, "%.3f", ImGuiSliderFlags_None);
            ImGui::PopID();
            ImGui::SetItemTooltip("Set type I error rate");

            ImGui::EndChild();
        }
        
        ImGui::SeparatorText("Results");
        {
            ImGuiTabBarFlags tab_bar_flags = ImGuiTabBarFlags_None;

            if (!option.auto_update && model.updater.requires_update && !model.updater.is_updating){

                //ImGui::PushStyleColor(ImGuiCol_Button, colours.base1());
                //ImGui::PushStyleColor(ImGuiCol_ButtonHovered, colours.base1(1));
                //ImGui::PushStyleColor(ImGuiCol_ButtonActive, colours.base1(2));
                //ImGui::PushStyleColor(ImGuiCol_Text, colours.base02(2));
    
                if (ImGui::Button("Click to update results", ImVec2(150, 30))){
                    model.checker.update();
                }
                ImGui::SetItemTooltip("As the design and model settings have changed, click this button to refresh the results");
                //ImGui::PopStyleColor(4);
            } else if(model.updater.is_updating) {
                //const ImU32 bg = ImGui::GetColorU32(ImGuiCol_Button);
                //ClusterApp::Spinner("##spinner", 15, 6, bg);
                ImGui::Text("Updating...");
            } else {
                if (ImGui::BeginTabBar("ResultsBar", tab_bar_flags))
                {
                    if (ImGui::BeginTabItem("Power (%)"))
                        {
                            if (ImGui::Button("Show power plot")){
                                model.checker.plotopen = true;
                                ImGui::OpenPopup("Power plot");
                            }
                            ImGui::SetItemTooltip("Click here to show plots of the power versus sample size and other parameters");
                             
            
                            if (ImGui::BeginPopupModal("Power plot", &model.checker.plotopen, ImGuiWindowFlags_AlwaysAutoResize)){

                                ImGui::Text("Plot settings");

                                const char* xaxis_items_exchangeable[] = { "Clusters", "N per cluster-period", "ICC", "Treatment effect", "Baseline" };
                                const char* xaxis_items_other[] = { "Clusters", "N per cluster-period", "ICC", "Treatment effect", "Baseline", "CAC" };
                                static int xaxis_item_current = 2;
                                static int series_item_current = 4;
                                const char* yaxis_items[] = { "Power", "Minimum detectable effect"};
                                static int yaxis_item_current = 0;
                                static int series_clusters[] = { 1,2,3 };
                                static int series_n[] = { 10,20,30 };
                                static float series_icc[] = { 0.01,0.02,0.05 };
                                static float series_te[] = { 0.0, 0.5, 1.0 };
                                static float series_baseline[] = { 0.0, 0.5, 1.0 };
                                static float series_cac[] = { 0.2,0.5,0.8 };
                                const short s16_one = 1;
                                static std::string series_label = "";
                                colourPicker colours;
                                static int print_prec = 3;
                                std::string x_label = "";
                                std::string ylabel = "";
                        
                                ImGui::SetNextItemWidth(200);
                                if (model.glmm.statmodel.covariance == ClusterApp::Covariance::exchangeable) {
                                    ImGui::Combo("X-Axis", &xaxis_item_current, xaxis_items_exchangeable, IM_ARRAYSIZE(xaxis_items_exchangeable));
                                }
                                else {
                                    ImGui::Combo("X-Axis", &xaxis_item_current, xaxis_items_other, IM_ARRAYSIZE(xaxis_items_other));
                                }
                                ImGui::SameLine(); HelpMarker("Clusters is the total number of clusters in the trial. The clusters are assumed to be allocated in the same ratio as the design (subject to rounding, which uses Hamilton's method).");
                                ImGui::SameLine();
                                ImGui::SetNextItemWidth(200);
                                ImGui::Combo("Y-Axis", &yaxis_item_current, yaxis_items, IM_ARRAYSIZE(yaxis_items));
                                ImGui::Text("X-axis limits: "); ImGui::SameLine();
                        
                                switch (xaxis_item_current) {
                                case 0:
                                    model.plot.xaxis = ClusterApp::XAxis::clusters;
                                    ImGui::SetNextItemWidth(100);
                                    ImGui::InputScalar("Lower", ImGuiDataType_S16, &model.plot.lower_int[0], &s16_one, NULL, "%d"); ImGui::SameLine();
                                    ImGui::SetNextItemWidth(100);
                                    ImGui::InputScalar("Upper", ImGuiDataType_S16, &model.plot.upper_int[0], &s16_one, NULL, "%d");
                                    print_prec = 0;
                                    x_label = "Clusters";
                                    break;
                                case 1:
                                    model.plot.xaxis = ClusterApp::XAxis::individual_n;
                                    ImGui::SetNextItemWidth(100);
                                    ImGui::InputScalar("Lower", ImGuiDataType_S16, &model.plot.lower_int[1], &s16_one, NULL, "%d"); ImGui::SameLine();
                                    ImGui::SetNextItemWidth(100);
                                    ImGui::InputScalar("Upper", ImGuiDataType_S16, &model.plot.upper_int[1], &s16_one, NULL, "%d");
                                    print_prec = 0;
                                    x_label = "n";
                                    break;
                                case 2:
                                    model.plot.xaxis = ClusterApp::XAxis::icc;
                                    ImGui::SetNextItemWidth(100);
                                    ImGui::DragFloat("Lower", &model.plot.lower_float[0], 0.001f, 0.0f, 1.0f, "%.3f", ImGuiSliderFlags_None); ImGui::SameLine();
                                    ImGui::SetNextItemWidth(100);
                                    ImGui::DragFloat("Upper", &model.plot.upper_float[0], 0.001f, 0.0f, 1.0f, "%.3f", ImGuiSliderFlags_None);
                                    print_prec = 3;
                                    x_label = "ICC";
                                    break;
                                case 3:
                                    model.plot.xaxis = ClusterApp::XAxis::treatment_effect;
                                    ImGui::SetNextItemWidth(100);
                                    ImGui::DragFloat("Lower", &model.plot.lower_float[1], 0.01f, -FLT_MAX, -FLT_MAX, "%.2f", ImGuiSliderFlags_None); ImGui::SameLine();
                                    ImGui::SetNextItemWidth(100);
                                    ImGui::DragFloat("Upper", &model.plot.upper_float[1], 0.01f, -FLT_MAX, -FLT_MAX, "%.2f", ImGuiSliderFlags_None);
                                    print_prec = 3;
                                    x_label = "Treatment effect";
                                    break;
                                case 4:
                                    model.plot.xaxis = ClusterApp::XAxis::baseline;
                                    ImGui::SetNextItemWidth(100);
                                    ImGui::DragFloat("Lower", &model.plot.lower_float[2], 0.01f, -FLT_MAX, -FLT_MAX, "%.2f", ImGuiSliderFlags_None); ImGui::SameLine();
                                    ImGui::SetNextItemWidth(100);
                                    ImGui::DragFloat("Upper", &model.plot.upper_float[2], 0.01f, -FLT_MAX, -FLT_MAX, "%.2f", ImGuiSliderFlags_None);
                                    print_prec = 3;
                                    x_label = "Baseline";
                                    break;
                                case 5:
                                    model.plot.xaxis = ClusterApp::XAxis::cac;
                                    ImGui::SetNextItemWidth(100);
                                    ImGui::DragFloat("Lower", &model.plot.lower_float[3], 0.001f, 0.0f, 1.0f, "%.3f", ImGuiSliderFlags_None); ImGui::SameLine();
                                    ImGui::SetNextItemWidth(100);
                                    ImGui::DragFloat("Upper", &model.plot.upper_float[3], 0.001f, 0.0f, 1.0f, "%.3f", ImGuiSliderFlags_None);
                                    print_prec = 3;
                                    x_label = "CAC";
                                    break;
                                }
                        
                                switch (yaxis_item_current) {
                                case 0:
                                    model.plot.yaxis = ClusterApp::YAxis::power;
                                    ylabel = "Power (%%)";
                                    break;
                                case 1:
                                    model.plot.yaxis = ClusterApp::YAxis::min_eff;
                                    ylabel = "Minimum detectable effect";
                                    break;
                                }
                        
                                ImGui::Checkbox("Multiple series?", &model.plot.multiple_series);
                                if (model.plot.multiple_series) {
                                    ImGui::Text("Add up to three series");
                                    ImGui::SetNextItemWidth(200);
                                    if (model.glmm.statmodel.covariance == ClusterApp::Covariance::exchangeable) {
                                        ImGui::Combo("Series variable", &series_item_current, xaxis_items_exchangeable, IM_ARRAYSIZE(xaxis_items_exchangeable));
                                    }
                                    else {
                                        ImGui::Combo("Series variable", &series_item_current, xaxis_items_other, IM_ARRAYSIZE(xaxis_items_other));
                                    }
                                    switch (series_item_current) {
                                    case 0:
                                        model.plot.series = ClusterApp::XAxis::clusters;
                                        series_label = "Clusters";
                                        model.plot.x_series[0] = (float)series_clusters[0];
                                        model.plot.x_series[1] = (float)series_clusters[1];
                                        model.plot.x_series[2] = (float)series_clusters[2];
                                        break;
                                    case 1:
                                        model.plot.series = ClusterApp::XAxis::individual_n;
                                        series_label = "n";
                                        model.plot.x_series[0] = (float)series_n[0];
                                        model.plot.x_series[1] = (float)series_n[1];
                                        model.plot.x_series[2] = (float)series_n[2];
                                        break;
                                    case 2:
                                        model.plot.series = ClusterApp::XAxis::icc;
                                        series_label = "ICC";
                                        model.plot.x_series[0] = series_icc[0];
                                        model.plot.x_series[1] = series_icc[1];
                                        model.plot.x_series[2] = series_icc[2];
                                        break;
                                    case 3:
                                        model.plot.series = ClusterApp::XAxis::treatment_effect;
                                        series_label = "Effect size";
                                        model.plot.x_series[0] = series_te[0];
                                        model.plot.x_series[1] = series_te[1];
                                        model.plot.x_series[2] = series_te[2];
                                        break;
                                    case 4:
                                        model.plot.series = ClusterApp::XAxis::baseline;
                                        series_label = "Baseline";
                                        model.plot.x_series[0] = series_baseline[0];
                                        model.plot.x_series[1] = series_baseline[1];
                                        model.plot.x_series[2] = series_baseline[2];
                                        break;
                                    case 5:
                                        model.plot.series = ClusterApp::XAxis::cac;
                                        series_label = "CAC";
                                        model.plot.x_series[0] = series_cac[0];
                                        model.plot.x_series[1] = series_cac[1];
                                        model.plot.x_series[2] = series_cac[2];
                                        break;
                                    }
                        
                                    //std::string label1 = std::to_string(model.plot.x_series[0]);
                        
                                    std::stringstream stream;
                                    stream << std::fixed << std::setprecision(3) << model.plot.x_series[0];
                                    std::string label1 = stream.str();
                                    char* char_array = new char[label1.length() + 1];
                                    strcpy(char_array, label1.c_str());
                                    char* label_array = new char[series_label.length() + 1];
                                    strcpy(label_array, series_label.c_str());
                        
                                    ImGui::PushID(10001);
                                    ImGui::PushStyleColor(ImGuiCol_Button, colours.red());
                                    ImGui::PushStyleColor(ImGuiCol_ButtonHovered, colours.red(1));
                                    ImGui::PushStyleColor(ImGuiCol_ButtonActive, colours.red(2));
                                    if (ImGui::Button(char_array, ImVec2(60, 40)))
                                        ImGui::OpenPopup("Series val 1");
                        
                                    if (ImGui::BeginPopupContextItem("Series val 1"))
                                    {
                                        ImGui::SetNextItemWidth(100);
                                        switch (series_item_current) {
                                        case 0:
                                            ImGui::InputScalar(label_array, ImGuiDataType_S16, &series_clusters[0], &s16_one, NULL, "%d");
                                            break;
                                        case 1:
                                            ImGui::InputScalar(label_array, ImGuiDataType_S16, &series_n[0], &s16_one, NULL, "%d");
                                            break;
                                        case 2:
                                            ImGui::DragFloat(label_array, &series_icc[0], 0.001f, 0.0f, 1.0f, "%.3f", ImGuiSliderFlags_None);
                                            break;
                                        case 3:
                                            ImGui::DragFloat(label_array, &series_te[0], 0.001f, 0.0f, 1.0f, "%.3f", ImGuiSliderFlags_None);
                                            break;
                                        case 4:
                                            ImGui::DragFloat(label_array, &series_baseline[0], 0.001f, 0.0f, 1.0f, "%.3f", ImGuiSliderFlags_None);
                                            break;
                                        case 5:
                                            ImGui::DragFloat(label_array, &series_cac[0], 0.001f, 0.0f, 1.0f, "%.3f", ImGuiSliderFlags_None);
                                            break;
                                        }
                        
                                        ImGui::EndPopup();
                                    }
                                    ImGui::PopID();
                                    ImGui::PopStyleColor(3);
                        
                                    if (model.plot.n_series > 1) {
                                        std::stringstream stream2;
                                        stream2 << std::fixed << std::setprecision(3) << model.plot.x_series[1];
                                        std::string label2 = stream2.str();
                                        char* char_array2 = new char[label2.length() + 1];
                                        strcpy(char_array2, label2.c_str());
                                        ImGui::SameLine();
                                        ImGui::PushID(10002);
                                        ImGui::PushStyleColor(ImGuiCol_Button, colours.blue());
                                        ImGui::PushStyleColor(ImGuiCol_ButtonHovered, colours.blue(1));
                                        ImGui::PushStyleColor(ImGuiCol_ButtonActive, colours.blue(2));
                                        if (ImGui::Button(char_array2, ImVec2(60, 40)))
                                            ImGui::OpenPopup("Series val 2");
                        
                                        if (ImGui::BeginPopupContextItem("Series val 2"))
                                        {
                                            ImGui::SetNextItemWidth(100);
                                            switch (series_item_current) {
                                            case 0:
                                                ImGui::InputScalar(label_array, ImGuiDataType_S16, &series_clusters[1], &s16_one, NULL, "%d");
                                                break;
                                            case 1:
                                                ImGui::InputScalar(label_array, ImGuiDataType_S16, &series_n[1], &s16_one, NULL, "%d");
                                                break;
                                            case 2:
                                                ImGui::DragFloat(label_array, &series_icc[1], 0.001f, 0.0f, 1.0f, "%.3f", ImGuiSliderFlags_None);
                                                break;
                                            case 3:
                                                ImGui::DragFloat(label_array, &series_te[1], 0.001f, 0.0f, 1.0f, "%.3f", ImGuiSliderFlags_None);
                                                break;
                                            case 4:
                                                ImGui::DragFloat(label_array, &series_baseline[1], 0.001f, 0.0f, 1.0f, "%.3f", ImGuiSliderFlags_None);
                                                break;
                                            case 5:
                                                ImGui::DragFloat(label_array, &series_cac[1], 0.001f, 0.0f, 1.0f, "%.3f", ImGuiSliderFlags_None);
                                                break;
                                            }
                        
                                            ImGui::EndPopup();
                                        }
                                        ImGui::PopID();
                                        ImGui::PopStyleColor(3);
                        
                                        if (model.plot.n_series == 3) {
                                            std::stringstream stream3;
                                            stream3 << std::fixed << std::setprecision(3) << model.plot.x_series[2];
                                            std::string label3 = stream3.str();
                                            char* char_array3 = new char[label3.length() + 1];
                                            strcpy(char_array3, label3.c_str());
                        
                                            ImGui::SameLine();
                                            ImGui::PushID(10003);
                                            ImGui::PushStyleColor(ImGuiCol_Button, colours.green());
                                            ImGui::PushStyleColor(ImGuiCol_ButtonHovered, colours.green(1));
                                            ImGui::PushStyleColor(ImGuiCol_ButtonActive, colours.green(2));
                                            if (ImGui::Button(char_array3, ImVec2(60, 40)))
                                                ImGui::OpenPopup("Series val 3");
                        
                                            if (ImGui::BeginPopupContextItem("Series val 3"))
                                            {
                                                ImGui::SetNextItemWidth(100);
                                                switch (series_item_current) {
                                                case 0:
                                                    ImGui::InputScalar(label_array, ImGuiDataType_S16, &series_clusters[2], &s16_one, NULL, "%d");
                                                    break;
                                                case 1:
                                                    ImGui::InputScalar(label_array, ImGuiDataType_S16, &series_n[2], &s16_one, NULL, "%d");
                                                    break;
                                                case 2:
                                                    ImGui::DragFloat(label_array, &series_icc[2], 0.001f, 0.0f, 1.0f, "%.3f", ImGuiSliderFlags_None);
                                                    break;
                                                case 3:
                                                    ImGui::DragFloat(label_array, &series_te[2], 0.001f, 0.0f, 1.0f, "%.3f", ImGuiSliderFlags_None);
                                                    break;
                                                case 4:
                                                    ImGui::DragFloat(label_array, &series_baseline[2], 0.001f, 0.0f, 1.0f, "%.3f", ImGuiSliderFlags_None);
                                                    break;
                                                case 5:
                                                    ImGui::DragFloat(label_array, &series_cac[2], 0.001f, 0.0f, 1.0f, "%.3f", ImGuiSliderFlags_None);
                                                    break;
                                                }
                        
                                                ImGui::EndPopup();
                                            }
                                            ImGui::PopID();
                                            ImGui::PopStyleColor(3);
                        
                                            ImGui::SameLine();
                                            ImGui::PushID(10223);
                                            if (ImGui::Button("-", ImVec2(20, 20))) {
                                                model.plot.n_series = 2;
                                            }
                                            ImGui::PopID();
                                        }
                                        else {
                        
                                            ImGui::SameLine();
                                            ImGui::PushID(10022);
                                            if (ImGui::Button("-", ImVec2(20, 20))) {
                                                model.plot.n_series = 1;
                                            }
                                            ImGui::PopID();
                        
                                            ImGui::SameLine();
                                            ImGui::PushID(10032);
                                            if (ImGui::Button("+", ImVec2(20, 20))) {
                                                model.plot.n_series = 3;
                                            }
                                            ImGui::PopID();
                                        }
                                    }
                                    else {
                                        ImGui::PushID(10012);
                                        ImGui::SameLine();
                                        if (ImGui::Button("+", ImVec2(20, 20))) {
                                            model.plot.n_series = 2;
                                        }
                                        ImGui::PopID();
                                    }
                        
                                }
                        
                                if (!model.plot.initialised)model.plot.update_data();
                        
                                char* x_char_array = new char[x_label.length() + 1];
                                strcpy(x_char_array, x_label.c_str());
                                char* y_char_array = new char[ylabel.length() + 1];
                                strcpy(y_char_array, ylabel.c_str());
                        
                                static ImPlotAxisFlags xflags = ImPlotAxisFlags_AutoFit; //| ImPlotAxisFlags_RangeFit
                                static ImPlotAxisFlags yflags = ImPlotAxisFlags_AutoFit;
                        
                                if (!option.auto_update && (model.updater.plot_requires_update)){
                                                           
                                    if (ImGui::Button("Refresh", ImVec2(80, 30))){
                                        model.checker.update();
                                    }
                                } else {
                                    if (!model.plot.updating) {
                                        if (ImPlot::BeginPlot("Cluster trial plot")) {
                                            ImPlot::SetupAxis(ImAxis_X1, x_char_array, xflags);
                                            ImPlot::SetupAxis(ImAxis_Y1, y_char_array, yflags);
                                            ImPlot::PushStyleColor(ImPlotCol_Line, colours.red());
                                            ImPlot::PlotLine("Series 1", model.plot.x_data, model.plot.y_data_1, model.plot.n_data_points);
                                            ImPlot::PopStyleColor();
                                            if (model.plot.multiple_series) {
                                                if (model.plot.n_series > 1) {
                                                    ImPlot::PushStyleColor(ImPlotCol_Line, colours.blue());
                                                    ImPlot::PlotLine("Series 2", model.plot.x_data, model.plot.y_data_2, model.plot.n_data_points);
                                                    ImPlot::PopStyleColor();
                                                    if (model.plot.n_series == 3) {
                                                        ImPlot::PushStyleColor(ImPlotCol_Line, colours.green());
                                                        ImPlot::PlotLine("Series 3", model.plot.x_data, model.plot.y_data_3, model.plot.n_data_points);
                                                        ImPlot::PopStyleColor();
                                                    }
                                                }
                                            }
                        
                                            ImPlot::EndPlot();
                                        }
                                }
                                else {
                                    ImGui::Text("Calculating...");
                                }
                                }

                                if (ImGui::Button("Close")){
                                    model.checker.plotopen = false;
                                }
                                ImGui::EndPopup();
                            }
                            
                            static ImGuiTableFlags flags = ImGuiTableFlags_Resizable | ImGuiTableFlags_Reorderable | ImGuiTableFlags_Hideable | ImGuiTableFlags_BordersOuter | ImGuiTableFlags_BordersV;
                            if(ImGui::BeginTable("resulttable2",2,flags)){
                                ImGui::TableSetupColumn("Value");
                                ImGui::TableSetupColumn("Design 1");
                                ImGui::TableHeadersRow();
                                ImGui::TableNextRow();
                                ImGui::TableSetColumnIndex(0);
                                ImGui::Text("Power (%%)");
                                ImGui::TableSetColumnIndex(1);
                                ImGui::Text("%.1f", model.updater.summary.power);

                                ImGui::TableNextRow();
                                ImGui::TableSetColumnIndex(0);
                                ImGui::Text("SE");
                                ImGui::TableSetColumnIndex(1);
                                ImGui::Text("%.2f", model.updater.summary.se);

                                ImGui::TableNextRow();
                                ImGui::TableSetColumnIndex(0);
                                ImGui::Text("DoF");
                                ImGui::TableSetColumnIndex(1);
                                ImGui::Text("%.1f", model.updater.summary.dof);

                                ImGui::TableNextRow();
                                ImGui::TableSetColumnIndex(0);
                                ImGui::Text("CI width");
                                ImGui::TableSetColumnIndex(1);
                                ImGui::Text("%.2f", model.updater.summary.ci_width);

                                ImGui::TableNextRow();
                                ImGui::TableSetColumnIndex(0);
                                ImGui::Text("Minimum detectable effect size");
                                ImGui::TableSetColumnIndex(1);
                                ImGui::Text("%.2f", model.updater.summary.min_eff);
                                ImGui::EndTable();
                            }
                            ImGui::EndTabItem();
                        }
                        if (ImGui::BeginTabItem("Sample size"))
                        {
                            
                            ImGui::Text("Find either the mean cluster-period size, or number of clusters, to provide the desired level of power.");
                            
                            ImGui::PushID(77770+model.id);
                            if(ImGui::Button("Find cluster-period size")){
                                model.target_ind_size = model.updater.sample_size_search(false, model.model.powertype);
                                model.find_period_size_trigger = true;
                                model.find_clusters_trigger = false;
                            }
                            ImGui::SetItemTooltip("Find the mean cluster size that will provide the target power");
                            ImGui::PopID();
                            ImGui::SameLine();
                            ImGui::PushID(77780+model.id);
                            if(ImGui::Button("Find number of clusters")){
                                model.target_cl_size = model.updater.sample_size_search(true, model.model.powertype);
                                model.find_period_size_trigger = false;
                                model.find_clusters_trigger = true;
                            }
                            ImGui::SetItemTooltip("Find the mean cluster size that will provide the target power");
                            ImGui::PopID();
                            ImGui::SameLine();

                            if (ImGui::Button("Show power surface")){
                                if (option.log)model.logger.AddLog("[%05d] [%s] Trigger power surface \n", ImGui::GetFrameCount(), model.logger.cat[0]);
                                model.openkrig = true;
                                model.kriggeropen = true;
                            }
                            ImGui::SetItemTooltip("Calculates and displays the 3D power surface for number of clusters versus cluster");
                            
                            
                            if(model.find_period_size_trigger){
                                ImGui::SeparatorText("Mean cluster-period size");
                                {
                                    ImGui::Text("Target sample size: %d", model.target_ind_size); ImGui::SameLine();
                                    ImGui::SetNextItemWidth(200);
                                    ImGui::PushID(77780+model.id);
                                    if(ImGui::Button("Set mean cluster size")){
                                        model.model.mean_size = model.target_ind_size;
                                    }
                                    ImGui::SetItemTooltip("Set mean cluster size to this value");
                                    ImGui::PopID();
                                }
                            }

                            if(model.find_clusters_trigger){
                                ImGui::SeparatorText("Number of clusters");
                                {
                                    ImGui::Text("Target number of clusters: %d", model.target_cl_size); ImGui::SameLine();
                                    ImGui::SetNextItemWidth(200);
                                    ImGui::PushID(77790+model.id);
                                    if(ImGui::Button("Set number of clusters")){
                                        int J = model.designs.total_clusters();
                                        std::vector<float> n_weights(model.designs.sequences);
                                        for (int i = 0; i < model.designs.sequences; i++) {
                                            n_weights[i] = (float)(*model.designs.n_clusters(i)) / (float)J;
                                        }
                                        std::vector<int> new_n_cl = model.glmm.round_weights(n_weights, model.target_cl_size);
                                        for (int k = 0; k < model.designs.sequences; k++) {
                                            *model.designs.n_clusters(k) = new_n_cl[k];
                                        }
                                    }
                                    ImGui::SetItemTooltip("Set number of clusters to this value");
                                    ImGui::PopID();
                                }
                            }
                            
                            
                            ImGui::EndTabItem();
                        }

                        /*if (ImGui::BeginTabItem("Effect size"))
                        {
                            ImGui::Text("This is the Cucumber tab!\nblah blah blah blah blah");
                            ImGui::EndTabItem();
                        }*/
                        ImGui::EndTabBar();
                }

                
            }


            if(model.kriggeropen){
                ImGui::OpenPopup("Power surface estimate");
            }
            if (ImGui::BeginPopupModal("Power surface estimate", &model.kriggeropen, ImGuiWindowFlags_AlwaysAutoResize))
            {
                ImGui::Text("Estimates the power across the sample size space");
                ImGui::TextWrapped("This feature is experimental. The power surface is estimated using a kriging approach and estimates will improve by sampling more points. The number of clusters is divided \
        according to the proportions given by the main trial design. Individuals are the number of individuals per cluster-period.");
                ImGui::TextWrapped("The power plot may not meet exactly the specified plot limits due to rounding of the increments in the plot and currently only provides power for GLS estimators.");
                
                const short  s16_zero = 0, s16_one = 1;
                static int new_sample_size = 25;
                colourPicker colours;

                if(!model.krig.start)model.krig.update(false);
                ImGui::SetNextItemWidth(150);
                ImGui::InputScalar("Sample size", ImGuiDataType_S16, &new_sample_size, &s16_one, NULL, "%d");
                ImGui::SetNextItemWidth(150);
                ImGui::InputScalar("Clusters min.", ImGuiDataType_S16, &model.krig.lower_int[0], &s16_one, NULL, "%d"); ImGui::SameLine();
                ImGui::SetNextItemWidth(150);
                ImGui::InputScalar("Clusters max.", ImGuiDataType_S16, &model.krig.upper_int[0], &s16_one, NULL, "%d");
                ImGui::SetNextItemWidth(150);
                ImGui::InputScalar("Individual min.", ImGuiDataType_S16, &model.krig.lower_int[1], &s16_one, NULL, "%d"); ImGui::SameLine();
                ImGui::SetNextItemWidth(150);
                ImGui::InputScalar("Individual max.", ImGuiDataType_S16, &model.krig.upper_int[1], &s16_one, NULL, "%d");
                if (ImGui::Button("Generate new sample")) {
                    model.krig.generate_grid();
                    model.krig.new_sample(new_sample_size);
                }

                ImGui::Text("Sample size = %lu", model.krig.n_ind.size());
                if(option.debug_info)ImGui::Text("Mu = %.3f", model.krig.mu);
                ImGui::Text("Sampled points");

                static ImPlotAxisFlags xflags = ImPlotAxisFlags_AutoFit; 
                static ImPlotAxisFlags yflags = ImPlotAxisFlags_AutoFit;

                if (ImPlot::BeginPlot("Scatter Plot", ImVec2(250, 250), ImPlotFlags_NoLegend | ImPlotFlags_NoMouseText)) {
                    ImPlot::SetupAxis(ImAxis_X1, "Clusters", xflags);
                    ImPlot::SetupAxis(ImAxis_Y1, "Individuals", yflags);
                    ImPlot::PushStyleColor(ImPlotCol_Line, colours.red());
                    ImPlot::PlotScatter("Sampled points", model.krig.n_cl.data(), model.krig.n_ind.data(), model.krig.n_ind.size());
                    ImPlot::EndPlot();
                }
                if (option.debug_info) {
                    ImGui::Text("Last sampled point: (%i", model.krig.n_ind[model.krig.n_ind.size() - 1]); ImGui::SameLine();
                    ImGui::Text(" , %i )", model.krig.n_cl[model.krig.n_cl.size() - 1]);
                }
                ImGui::SetNextItemWidth(200);
                ImGui::DragFloat("Bandwidth", &model.krig.bandwidth, 0.01f, 0.0f, 2.0f, "%.2f", ImGuiSliderFlags_None); ImGui::SameLine(); HelpMarker(
                    "The bandwidth controls the level of smoothing for the power surface. Sample size values are scaled to [0,1] x [0,1].");
                
                ImGui::SetNextItemWidth(200);
                ImGui::DragFloat("Power threshold", &model.krig.threshold_power, 0.01f, 0.0f, 1.0f, "%0.2f", ImGuiSliderFlags_None); ImGui::SameLine(); HelpMarker(
                    "The points are sampled to minimised the uncertainty about where the contour for this level of power is.");
        
                if (ImGui::Button("Update heat map")) {
                    model.krig.update(false);
                }
                ImGui::SameLine();
                if (ImGui::Button("Sample new point")) {
                    model.krig.update(true);
                }

                if (model.krig.surface_initialised) {
                    static ImPlotColormap map = ImPlotColormap_Viridis;
                    if (ImPlot::ColormapButton(ImPlot::GetColormapName(map), ImVec2(225, 0), map)) {
                        map = (map + 1) % ImPlot::GetColormapCount();
                        ImPlot::BustColorCache("Power plot");
                    }

                    ImPlot::PushColormap(map);

                    static ImPlotAxisFlags axes_flags = ImPlotAxisFlags_Lock | ImPlotAxisFlags_NoGridLines;
                    if (ImPlot::BeginPlot("Power plot", ImVec2(800,500), ImPlotFlags_NoLegend | ImPlotFlags_NoMouseText)) {
                        ImPlot::SetupAxes("Clusters", "Individuals", axes_flags, axes_flags);
                        ImPlot::SetupAxisTicks(ImAxis_X1, 0 + 1.0 / 40.0, 1 - 1.0 / 40.0, 20, model.krig.n_cl_grid_label);
                        ImPlot::SetupAxisTicks(ImAxis_Y1, 0 + 1.0 / 40.0, 1 - 1.0 / 40.0, 20, model.krig.n_ind_grid_label);
                        ImPlot::PlotHeatmap("Power", model.krig.surface, 20,20,0,0,"%.0f");
                        ImPlot::EndPlot();
                    }
                }
                if (ImGui::Button("Close")){
                    model.kriggeropen = false;
                }
                        
                ImGui::EndPopup();
            }
            
            
        }

        ImGui::End();
    };

    
}

