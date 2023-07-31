#include "clusterapp.h"

namespace ClusterApp {

    void ShowMainMenu(ClusterApp::options& windows) {
        if (ImGui::BeginMainMenuBar())
        {
            if (ImGui::BeginMenu("File"))
            {
                if (ImGui::BeginMenu("Options"))
                {
                    ImGui::Checkbox("Light model", &windows.light_mode);
                    ImGui::EndMenu();
                }
                ImGui::EndMenu();
            }
            if (ImGui::BeginMenu("Trial options"))
            {
                ImGui::Checkbox("Two treatments", &windows.two_treatments);
                ImGui::EndMenu();
            }
            if (ImGui::BeginMenu("Windows"))
            {
                ImGui::Checkbox("Sample size", &windows.sample_size);
                ImGui::Checkbox("Statistical model", &windows.model);
                ImGui::Checkbox("Results", &windows.results);
                ImGui::Checkbox("Optimal design", &windows.optimiser);
                ImGui::EndMenu();
            }

            ImGui::EndMainMenuBar();
        }
    }

    void RenderDesigner(ClusterApp::design& designs, ClusterApp::modelUpdater& updater, ClusterApp::options& option) {
        ImGui::Begin("Trial Designer", NULL, ImGuiWindowFlags_MenuBar);
        if (ImGui::BeginMenuBar()) {
            if (ImGui::BeginMenu("Designs")) {
                if (ImGui::BeginMenu("Parallel")) {
                    static int parallel_t = 1;
                    static int parallel_n = 10;
                    static int parallel_J = 1;
                    ImGui::Text("Trial details");
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
                        designs.set_parallel(parallel_t, parallel_n, parallel_J);
                        updater.update_data();
                    }
                    ImGui::EndMenu();
                }
                if (ImGui::BeginMenu("Parallel with baseline")) {
                    static int parallelb_t = 1;
                    static int parallelb_n = 10;
                    static int parallelb_J = 1;
                    ImGui::Text("Trial details");
                    ImGui::SetNextItemWidth(100);
                    ImGui::InputInt("Time periods", &parallelb_t, 1, 10, 0);
                    ImGui::SetNextItemWidth(100);
                    ImGui::InputInt("Observations per cluster-period", &parallelb_n, 1, 10, 0);
                    ImGui::SetNextItemWidth(100);
                    ImGui::InputInt("Clusters per arm", &parallelb_J, 1, 10, 0);
                    if (ImGui::Button("Set")) {
                        designs.set_parallel_with_baseline(parallelb_t, parallelb_n, parallelb_J);
                        updater.update_data();
                    }
                    ImGui::EndMenu();
                }
                if (ImGui::BeginMenu("Stepped-wedge")) {
                    static int wedge_t = 6;
                    static int wedge_n = 10;
                    static int wedge_J = 1;
                    ImGui::Text("Trial details");
                    ImGui::SetNextItemWidth(100);
                    ImGui::InputInt("Time periods", &wedge_t, 1, 10, 0);
                    ImGui::SetNextItemWidth(100);
                    ImGui::InputInt("Observations per cluster-period", &wedge_n, 1, 10, 0);
                    ImGui::SetNextItemWidth(100);
                    ImGui::InputInt("Clusters per arm", &wedge_J, 1, 10, 0);
                    if (ImGui::Button("Set")) {
                        designs.set_stepped_wedge(wedge_t, wedge_n, wedge_J);
                        updater.update_data();
                    }
                    ImGui::EndMenu();
                }
                if (ImGui::BeginMenu("Crossover")) {
                    static int cross_n;
                    static int cross_J;
                    ImGui::Text("Trial details");
                    ImGui::SetNextItemWidth(100);
                    ImGui::InputInt("Observations per cluster-period", &cross_n, 1, 10, 0);
                    ImGui::SetNextItemWidth(100);
                    ImGui::InputInt("Clusters per arm", &cross_J, 1, 10, 0);
                    if (ImGui::Button("Set")) {
                        designs.set_crossover(cross_n, cross_J);
                        updater.update_data();
                    }
                    ImGui::EndMenu();
                }
                if (ImGui::BeginMenu("Staircase")) {
                    static int staircase_t;
                    static int staircase_n;
                    static int staircase_J;
                    ImGui::Text("Trial details");
                    ImGui::SetNextItemWidth(100);
                    ImGui::InputInt("Time periods", &staircase_t, 1, 10, 0);
                    ImGui::SetNextItemWidth(100);
                    ImGui::InputInt("Observations per cluster-period", &staircase_n, 1, 10, 0);
                    ImGui::SetNextItemWidth(100);
                    ImGui::InputInt("Clusters per arm", &staircase_J, 1, 10, 0);
                    if (ImGui::Button("Set")) {
                        designs.set_staircase(staircase_t, staircase_n, staircase_J);
                        updater.update_data();
                    }
                    ImGui::EndMenu();
                }
                if (ImGui::BeginMenu("Factorial")) {
                    if (option.two_treatments) {
                        static int factorial_t;
                        static int factorial_n;
                        static int factorial_J;
                        ImGui::Text("Trial details");
                        ImGui::SetNextItemWidth(100);
                        ImGui::InputInt("Time periods", &factorial_t, 1, 10, 0);
                        ImGui::SetNextItemWidth(100);
                        ImGui::InputInt("Observations per cluster-period", &factorial_n, 1, 10, 0);
                        ImGui::SetNextItemWidth(100);
                        ImGui::InputInt("Clusters per arm", &factorial_J, 1, 10, 0);
                        if (ImGui::Button("Set")) {
                            designs.set_factorial(factorial_t, factorial_n, factorial_J);
                            updater.update_data();
                        }
                    }
                    else {
                        ImGui::Text("To enable this option, set two treatments in trial options.");
                    }
                    ImGui::EndMenu();
                }
                ImGui::EndMenu();
            }
            if (ImGui::BeginMenu("Modify")) {
                if (ImGui::MenuItem("Activate all")) {
                    for (int i = 0; i < designs.sequences; i++) {
                        for (int t = 0; t < designs.time; t++) {
                            *designs.active(i, t) = true;
                        }
                    }
                }
                if (ImGui::MenuItem("Deactivate all")) {
                    for (int i = 0; i < designs.sequences; i++) {
                        for (int t = 0; t < designs.time; t++) {
                            *designs.active(i, t) = false;
                        }
                    }
                }
                if (ImGui::BeginMenu("Time periods")) {
                    static int total_t;
                    ImGui::SetNextItemWidth(100);
                    ImGui::InputInt("Time periods", &total_t, 1, 10, 0);
                    if (ImGui::Button("Set")) {
                        if (total_t < designs.time && total_t > 0) {
                            int difft = designs.time - total_t;
                            for (int i = 0; i < difft; i++) {
                                designs.remove_period(designs.time - 1);
                            }
                        }
                        else if (total_t > designs.time) {
                            int difft = total_t - designs.time;
                            for (int i = 0; i < difft; i++) {
                                designs.add_period();
                            }
                        }
                    }
                    ImGui::EndMenu();
                }
                ImGui::EndMenu();
            }
            if (ImGui::BeginMenu("Sequences")) {
                static int total_t;
                ImGui::SetNextItemWidth(100);
                ImGui::InputInt("Sequences", &total_t, 1, 10, 0);
                if (ImGui::Button("Set")) {
                    if (total_t < designs.sequences && total_t > 0) {
                        int difft = designs.sequences - total_t;
                        for (int i = 0; i < difft; i++) {
                            designs.remove_sequence(designs.sequences - 1);
                        }
                    }
                    else if (total_t < designs.sequences) {
                        int difft = total_t - designs.sequences;
                        for (int i = 0; i < difft; i++) {
                            designs.add_sequence();
                        }
                    }
                }
                if (ImGui::MenuItem("Split sequences into clusters")) {
                    designs.split_sequences();
                }
                if (ImGui::MenuItem("Combine clusters into sequences")) {
                    designs.combine_sequences();
                }
                ImGui::EndMenu();
            }
            if (ImGui::BeginMenu("Clusters")) {
                static int total_t;
                static int total_n;
                ImGui::SetNextItemWidth(100);
                ImGui::InputInt("Clusters per sequence", &total_t, 1, 10, 0); ImGui::SameLine();
                if (ImGui::Button("Set")) {
                    if (total_t > 0) {
                        for (int j = 0; j < designs.sequences; j++)  *designs.n_clusters(j) = total_t;
                    }
                }
                ImGui::SetNextItemWidth(100);
                ImGui::InputInt("n per cluster-period", &total_n, 1, 10, 0); ImGui::SameLine();
                if (ImGui::Button("Set")) {
                    if (total_n > 0) {
                        for (int j = 0; j < designs.sequences; j++) {
                            for (int t = 0; t < designs.time; t++) {
                                *designs.n(j, t) = total_n;
                            }
                        }
                    }
                    ImGui::EndMenu();
                }
                ImGui::EndMenu();
            }
            ImGui::EndMenuBar();
        }
        ImGuiStyle& style = ImGui::GetStyle();
        colourPicker colours;
        int horiztonal_align = 70;
        const short   s16_zero = 0, s16_one = 1;
        static bool remove_seq = false;
        static int which_remove = 0;
        static int set_all_n = 0;
        static int max_large_dim = 80;
        static int large_dim;
        static int small_dim = 20;
        float window_visible_x2 = ImGui::GetWindowContentRegionMax().x; //ImGui::GetWindowPos().x + 
        float large_test_size = (window_visible_x2 - (3 * small_dim) - ((designs.time + 5) * style.ItemSpacing.x)) * (1.0f / designs.time);
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



        ImGui::TextWrapped("Design the trial below. Rows are sequences, columns are time periods. Click + to add new sequences or time periods. Select cells to edit their details. Select row or column headers to change the numbers of clusters.");
        ImGui::Spacing();

        ImGui::Dummy(ImVec2(large_dim, small_dim)); ImGui::SameLine(horiztonal_align + 30);
        ImGui::Text("TIME");

        //ImGui::Text("SEQUENCE"); ImGui::SameLine(horiztonal_align);
        ImGui::Dummy(ImVec2(small_dim, small_dim)); ImGui::SameLine();
        ImGui::Dummy(ImVec2(small_dim, small_dim)); ImGui::SameLine();
        ImGui::PushStyleColor(ImGuiCol_Button, colours.base1());
        ImGui::PushStyleColor(ImGuiCol_ButtonHovered, colours.base1(1));
        ImGui::PushStyleColor(ImGuiCol_ButtonActive, colours.base1(2));
        ImGui::PushStyleColor(ImGuiCol_Text, colours.base02(2));
        for (int t = 0; t <= designs.time; t++) {
            ImGui::PushID(designs.time * designs.sequences + t);
            if (ImGui::Button("+", ImVec2(small_dim, small_dim))) {
                if (t == designs.time) {
                    designs.add_period();
                }
                else {
                    designs.add_period(t);
                }
            };
            ImGui::PopID();
            ImGui::SameLine();
            ImGui::Dummy(ImVec2(large_dim - small_dim - style.ItemSpacing[0], small_dim));
            if (t < (designs.time))ImGui::SameLine();
        }

        ImGui::Dummy(ImVec2(small_dim, small_dim)); ImGui::SameLine();
        ImGui::Dummy(ImVec2(small_dim, small_dim)); ImGui::SameLine();
        for (int t = 0; t < designs.time; t++) {
            ImGui::PushID(designs.time * designs.sequences + designs.time + 2 + designs.sequences + t);
            ImGui::Button(int_to_char(t + 1), ImVec2(large_dim, small_dim));
            if (ImGui::BeginPopupContextItem(NULL, ImGuiPopupFlags_MouseButtonLeft))
            {
                ImGui::Text("OPTIONS");
                if (ImGui::SmallButton("Activate all")) {
                    for (int s = 0; s < designs.sequences; s++)*designs.active(s, t) = true;
                }
                if (ImGui::SmallButton("Set all intervention")) {
                    for (int s = 0; s < designs.sequences; s++) {
                        *designs.intervention(s, t) = true;
                        *designs.intervention_2(s, t) = false;
                    }
                }
                if (option.two_treatments) {
                    if (ImGui::SmallButton("Set all intervention 2")) {
                        for (int s = 0; s < designs.sequences; s++) {
                            *designs.intervention(s, t) = false;
                            *designs.intervention_2(s, t) = true;
                        }
                    }
                    if (ImGui::SmallButton("Set all intervention 1 & 2")) {
                        for (int s = 0; s < designs.sequences; s++) {
                            *designs.intervention(s, t) = true;
                            *designs.intervention_2(s, t) = true;
                        }
                    }
                }
                if (ImGui::SmallButton("Set all control")) {
                    for (int s = 0; s < designs.sequences; s++)*designs.intervention(s, t) = false;
                }
                if (designs.time > 1) {
                    if (ImGui::SmallButton("Remove column")) {
                        designs.remove_period(t);
                    }
                }
                ImGui::EndPopup();
            }
            ImGui::PopID();
            if (t < (designs.time - 1))ImGui::SameLine();
        }
        ImGui::PopStyleColor(4);

        for (int n = 0; n < designs.sequences; n++)
        {
            //ImGui::Dummy(ImVec2(large_dim, small_dim)); ImGui::SameLine(horiztonal_align);
            ImGui::PushStyleColor(ImGuiCol_Button, colours.base1());
            ImGui::PushStyleColor(ImGuiCol_ButtonHovered, colours.base1(1));
            ImGui::PushStyleColor(ImGuiCol_ButtonActive, colours.base1(2));
            ImGui::PushStyleColor(ImGuiCol_Text, colours.base02(2));
            ImGui::PushID(designs.time * designs.sequences + designs.time + 1 + n);
            if (ImGui::Button("+", ImVec2(small_dim, small_dim))) {
                designs.add_sequence(n);
            }
            ImGui::SameLine();
            ImGui::PopID();
            ImGui::PushID(designs.time * (designs.sequences + 1) + designs.time + designs.sequences + 2 + n);
            ImGui::Button(int_to_char(n + 1), ImVec2(small_dim, large_dim)); ImGui::SameLine();
            if (ImGui::BeginPopupContextItem(NULL, ImGuiPopupFlags_MouseButtonLeft))
            {
                ImGui::Text("Number of clusters");
                ImGui::SetNextItemWidth(100);
                ImGui::InputScalar("N", ImGuiDataType_S16, designs.n_clusters(n), &s16_one, NULL, "%d");
                ImGui::Text("OPTIONS");
                if (ImGui::SmallButton("Activate all")) {
                    for (int s = 0; s < designs.time; s++)*designs.active(n, s) = true;
                }
                if (ImGui::SmallButton("Set all intervention")) {
                    for (int s = 0; s < designs.time; s++) {
                        *designs.intervention(n, s) = true;
                        *designs.intervention_2(n, s) = false;
                    }
                }
                if (option.two_treatments) {
                    if (ImGui::SmallButton("Set all intervention 2")) {
                        for (int s = 0; s < designs.time; s++) {
                            *designs.intervention(n, s) = false;
                            *designs.intervention_2(n, s) = true;
                        }
                    }
                    if (ImGui::SmallButton("Set all intervention 1 & 2")) {
                        for (int s = 0; s < designs.time; s++) {
                            *designs.intervention(n, s) = true;
                            *designs.intervention_2(n, s) = true;
                        }
                    }
                }
                if (ImGui::SmallButton("Set all control")) {
                    for (int s = 0; s < designs.time; s++)*designs.intervention(n, s) = false;
                }
                if (designs.sequences > 1) {
                    if (ImGui::SmallButton("Remove sequence")) {
                        //designs.remove_sequence(n);
                        remove_seq = true;
                        which_remove = n;
                    }
                }
                ImGui::SetNextItemWidth(100);
                ImGui::InputScalar("n", ImGuiDataType_S16, &set_all_n, NULL, NULL, "%d"); ImGui::SameLine();
                if (ImGui::SmallButton("Set all n")) {
                    for (int s = 0; s < designs.time; s++)*designs.n(n, s) = set_all_n;
                }

                ImGui::EndPopup();
            }
            ImGui::PopID();
            ImGui::PopStyleColor(4);
            for (int t = 0; t < designs.time; t++) {
                ImGui::PushID(t + n * designs.time);
                std::string label = "";
                if (*designs.active(n, t)) {
                    label = std::to_string(*designs.intervention(n, t));
                    if (option.two_treatments)label += "/" + std::to_string(*designs.intervention_2(n, t));
                    if (*designs.intervention(n, t) && !*designs.intervention_2(n, t)) {
                        ImGui::PushStyleColor(ImGuiCol_Button, colours.red());
                        ImGui::PushStyleColor(ImGuiCol_ButtonHovered, colours.red(1));
                        ImGui::PushStyleColor(ImGuiCol_ButtonActive, colours.red(2));
                    }
                    else if (!*designs.intervention(n, t) && *designs.intervention_2(n, t)) {
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
                if (!option.two_treatments || (option.two_treatments && !(*designs.intervention(n, t) && *designs.intervention_2(n, t)))) {
                    ImGui::Button(char_array, ImVec2(large_dim, large_dim));
                    // popup context menu for cluster-period cells
                    if (ImGui::BeginPopupContextItem(NULL, ImGuiPopupFlags_MouseButtonLeft))
                    {
                        ImGui::Checkbox("Active", designs.active(n, t));
                        ImGui::Checkbox("Intervention", designs.intervention(n, t));
                        if (option.two_treatments)ImGui::Checkbox("Intervention 2", designs.intervention_2(n, t));
                        ImGui::SetNextItemWidth(100);
                        ImGui::InputScalar("N", ImGuiDataType_S16, designs.n(n, t), &s16_one, NULL, "%d");
                        ImGui::EndPopup();
                    }
                }
                else {
                    ClusterApp::MultiColorButton(char_array, colours.yellow(), colours.yellow(), colours.red(), colours.red(), 0, 0, ImVec2(large_dim, large_dim),
                        colours.yellow(1), colours.yellow(1), colours.red(1), colours.red(1),
                        colours.yellow(2), colours.yellow(2), colours.red(2), colours.red(2));
                    if (ImGui::BeginPopupContextItem(NULL, ImGuiPopupFlags_MouseButtonLeft))
                    {
                        ImGui::Checkbox("Active", designs.active(n, t));
                        ImGui::Checkbox("Intervention", designs.intervention(n, t));
                        ImGui::Checkbox("Intervention 2", designs.intervention_2(n, t));
                        ImGui::SetNextItemWidth(100);
                        ImGui::InputScalar("N", ImGuiDataType_S16, designs.n(n, t), &s16_one, NULL, "%d");
                        ImGui::EndPopup();
                    }
                }
                ImGui::PopStyleColor(3);

                ImGui::PopID();
                if (t < (designs.time - 1))ImGui::SameLine();
            }
        }

        ImGui::PushStyleColor(ImGuiCol_Button, colours.base1());
        ImGui::PushStyleColor(ImGuiCol_ButtonHovered, colours.base1(1));
        ImGui::PushStyleColor(ImGuiCol_ButtonActive, colours.base1(2));
        ImGui::PushStyleColor(ImGuiCol_Text, colours.base02(2));
        ImGui::PushID(designs.time * designs.sequences + designs.time + 1 + designs.sequences);
        if (ImGui::Button("+", ImVec2(small_dim, small_dim))) {
            designs.add_sequence();
        }
        ImGui::PopID();
        ImGui::PopStyleColor(4);


        if (remove_seq) {
            designs.remove_sequence(which_remove);
            remove_seq = false;
        }

        ImGui::Text("Design Checksum"); ImGui::SameLine();
        ImGui::Text("%d", designs.crc_val);
        ImGui::End();
    }

    void RenderSampleSize(ClusterApp::design& designs) {
        ImGui::Begin("Trial Sample Size");

        static ImGuiTableFlags flags = ImGuiTableFlags_Borders | ImGuiTableFlags_Resizable | ImGuiTableFlags_Reorderable | ImGuiTableFlags_Hideable;// | ImGuiTableFlags_NoHostExtendX;
        //ImGui::CheckboxFlags("ImGuiTableFlags_NoHostExtendX", &flags, ImGuiTableFlags_NoHostExtendX);ImGuiTableFlags_Borders | ImGuiTableFlags_Resizable | ImGuiTableFlags_Reorderable | ImGuiTableFlags_Hideable

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
            for (int i = 0; i < designs.sequences; i++) {
                ImGui::Text(int_to_char(i + 1));
                ImGui::TableNextColumn();
                ImGui::Text(int_to_char(*designs.n_clusters(i)));
                ImGui::TableNextColumn();
                int total_n = 0;
                if (ImGui::BeginTable("summaryn", designs.time, ImGuiTableFlags_Borders | ImGuiTableFlags_Resizable | ImGuiTableFlags_Reorderable | ImGuiTableFlags_Hideable)) {
                    for (int t = 0; t < designs.time; t++) {
                        ImGui::TableSetupColumn(int_to_char(t + 1));
                    }
                    ImGui::TableHeadersRow();
                    ImGui::TableNextColumn();
                    for (int j = 0; j < *designs.n_clusters(i); j++) {
                        for (int t = 0; t < designs.time; t++) {
                            int nper = 0;
                            if (*designs.active(i, t))nper = *designs.n(i, t);
                            ImGui::Text(int_to_char(nper));
                            total_n += nper;
                            ImGui::TableNextColumn();
                        }
                    }
                    ImGui::EndTable();
                }
                ImGui::TableNextColumn();
                ImGui::Text(int_to_char(total_n));
                total_total_n += total_n;
                ImGui::TableNextColumn();
            }
            ImGui::Text("TOTAL");
            ImGui::TableNextColumn();
            ImGui::Text(int_to_char(designs.total_clusters()));
            ImGui::TableNextColumn(); ImGui::TableNextColumn();
            ImGui::Text(int_to_char(total_total_n));
            ImGui::EndTable();
        }

        ImGui::End();
    }

    void RenderModel(ClusterApp::design& design, ClusterApp::statisticalModel& model, ClusterApp::options& option, ImFont* unifont) {
        ImGui::Begin("Statistical Model");

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
            ImGui::RadioButton("Cohort", &structure_sampling, 1);
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

                if (family_item_current == 0 && link_item_current == 0) {
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
                        ImGui::PushFont(unifont);
                        ImGui::DragFloat(u8"\u03BB", &model.cov_pars[1], 0.01f, 0.0f, 1.0f, "%.3f", ImGuiSliderFlags_None); ImGui::SameLine(); ImGui::PopFont(); HelpMarker(
                            "Cluster-period autoregressive parameter.");
                    }
                    if (cl_cov_item_current == 3 || cl_cov_item_current == 4) {
                        ImGui::SetNextItemWidth(200);
                        ImGui::PushFont(unifont);
                        ImGui::DragFloat(u8"\u03B8", &model.cov_pars[1], 0.01f, 0.0f, +FLT_MAX, "%.3f", ImGuiSliderFlags_None); ImGui::SameLine(); ImGui::PopFont(); HelpMarker(
                            u8"Cluster-period denominator parameter, e.g. exp(\u03B4/\u03B8) for exponential covariance.");
                    }
                    if (structure_sampling == 1) {
                        ImGui::SetNextItemWidth(200);
                        ImGui::DragFloat("IAC", &model.ixx_pars[2], 0.001f, 0.0f, 1.0f, "%.3f", ImGuiSliderFlags_None); ImGui::SameLine();  HelpMarker(
                            "Individual autocorrelation coefficient.");
                        if (ind_cov_item_current == 1) {
                            ImGui::SetNextItemWidth(200);
                            ImGui::PushFont(unifont);
                            ImGui::DragFloat(u8"\u03BB (individual)", &model.cov_pars[4], 0.01f, 0.0f, 1.0f, "%.3f", ImGuiSliderFlags_None); ImGui::SameLine(); ImGui::PopFont(); HelpMarker(
                                "Individual autoregressive parameter.");
                        }
                        if (ind_cov_item_current == 2 || ind_cov_item_current == 3) {
                            ImGui::SetNextItemWidth(200);
                            ImGui::PushFont(unifont);
                            ImGui::DragFloat(u8"\u03B8 (individual)", &model.cov_pars[4], 0.01f, 0.0f, +FLT_MAX, "%.3f", ImGuiSliderFlags_None); ImGui::SameLine(); ImGui::PopFont(); HelpMarker(
                                u8"Individual denominator parameter, e.g. exp(\u03B4/\u03B8) for exponential covariance.");
                        }
                    }
                }
                else {
                    ImGui::TextWrapped("Use of ICC, CAC, and IAC values is limited to Gaussian-identity models currently, as these values are not constants. Please set the variance parameters directly.");
                    ImGui::SetNextItemWidth(200);
                    ImGui::PushFont(unifont);
                    ImGui::DragFloat(u8"\u03C4\u00B2", &model.cov_pars[0], 0.01f, 0.0f, +FLT_MAX, "%.3f", ImGuiSliderFlags_None); ImGui::SameLine(); ImGui::PopFont(); HelpMarker(
                        "Cluster-level variance parameter.");
                    if (cl_cov_item_current == 1) {
                        ImGui::SetNextItemWidth(200);
                        ImGui::PushFont(unifont);
                        ImGui::DragFloat(u8"\u03C9\u00B2", &model.cov_pars[1], 0.01f, 0.0f, +FLT_MAX, "%.3f", ImGuiSliderFlags_None); ImGui::SameLine(); ImGui::PopFont(); HelpMarker(
                            "Cluster-period level variance parameter.");
                    }
                    if (cl_cov_item_current == 2) {
                        ImGui::SetNextItemWidth(200);
                        ImGui::PushFont(unifont);
                        ImGui::DragFloat(u8"\u03BB", &model.cov_pars[1], 0.01f, 0.0f, 1.0f, "%.3f", ImGuiSliderFlags_None); ImGui::SameLine(); ImGui::PopFont(); HelpMarker(
                            "Cluster-period autoregressive parameter.");
                    }
                    if (cl_cov_item_current == 3 || cl_cov_item_current == 4) {
                        ImGui::SetNextItemWidth(200);
                        ImGui::PushFont(unifont);
                        ImGui::DragFloat(u8"\u03B8", &model.cov_pars[1], 0.01f, 0.0f, +FLT_MAX, "%.3f", ImGuiSliderFlags_None); ImGui::SameLine(); HelpMarker(
                            u8"Cluster-period denominator parameter, e.g. exp(\u03B4/\u03B8) for exponential covariance.");  ImGui::PopFont();
                    }
                    if (structure_sampling == 1) {
                        ImGui::SetNextItemWidth(200);
                        ImGui::PushFont(unifont);
                        ImGui::DragFloat(u8"\u03B7", &model.cov_pars[3], 0.01f, 0.0f, +FLT_MAX, "%.3f", ImGuiSliderFlags_None); ImGui::SameLine(); ImGui::PopFont(); HelpMarker(
                            "Individual-level variance term");
                        if (ind_cov_item_current == 2) {
                            ImGui::SetNextItemWidth(200);
                            ImGui::PushFont(unifont);
                            ImGui::DragFloat(u8"\u03BB (individual)", &model.cov_pars[4], 0.01f, 0.0f, 1.0f, "%.3f", ImGuiSliderFlags_None); ImGui::SameLine(); ImGui::PopFont(); HelpMarker(
                                "Individual autoregressive parameter.");
                        }
                        if (ind_cov_item_current == 3 || ind_cov_item_current == 4) {
                            ImGui::SetNextItemWidth(200);
                            ImGui::PushFont(unifont);
                            ImGui::DragFloat(u8"\u03B8 (individual)", &model.cov_pars[4], 0.01f, 0.0f, +FLT_MAX, "%.3f", ImGuiSliderFlags_None); ImGui::SameLine();  HelpMarker(
                                u8"Individual denominator parameter, e.g. exp(\u03B4/\u03B8) for exponential covariance.");  ImGui::PopFont();
                        }
                    }
                    if (family_item_current == 0 || family_item_current == 3 || family_item_current == 4) {
                        ImGui::SetNextItemWidth(200);
                        ImGui::DragFloat(u8"\u03C3\u00B2", &model.cov_pars[2], 0.01f, 0.0f, +FLT_MAX, "%.3f", ImGuiSliderFlags_None); ImGui::SameLine();  HelpMarker(
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
            ImGui::TreePop();
        }

        if (ImGui::TreeNode("Model info")) {
            ImGui::Text("Model checksum"); ImGui::SameLine();
            ImGui::Text("%d", model.crc_val);
            ImGui::Text("Parameters checksum"); ImGui::SameLine();
            ImGui::Text("%d", model.crc_val_pars);
            ImGui::TreePop();
        }

        ImGui::End();
    };

    void RenderResults(ClusterApp::modelUpdater& updater, ClusterApp::options& option) {
        ImGui::Begin("Results");
        ImGui::Text("RESULTS");

        if (updater.update) {
            ImGui::SameLine(); ImGui::Text("Recalculating...");
        }
        static ImGuiTableFlags flags = ImGuiTableFlags_Borders | ImGuiTableFlags_Resizable | ImGuiTableFlags_Reorderable | ImGuiTableFlags_Hideable | ImGuiTableFlags_NoHostExtendX;

        if (ImGui::BeginTable("results", 4, flags))
        {
            ImGui::TableSetupColumn("Statistic");
            ImGui::TableSetupColumn("GLS");
            ImGui::TableSetupColumn("GLS (between-within)");
            ImGui::TableSetupColumn("Kenward-Roger");
            ImGui::TableHeadersRow();
            ImGui::TableNextColumn();

            ImGui::Text("Power (%%)");
            ImGui::TableNextColumn();
            ImGui::Text("%.1f", updater.summary.power);
            ImGui::TableNextColumn();
            ImGui::Text("%.1f", updater.summary.power_bw);
            ImGui::TableNextColumn();
            ImGui::Text("%.1f", updater.summary.power_kr);
            ImGui::TableNextColumn();

            ImGui::Text("Confidence interval half-width");
            ImGui::TableNextColumn();
            ImGui::Text("%.3f", updater.summary.ci_width);
            ImGui::TableNextColumn();
            ImGui::Text("%.3f", updater.summary.ci_width_bw);
            ImGui::TableNextColumn();
            ImGui::Text("%.3f", updater.summary.ci_width_kr);
            ImGui::TableNextColumn();

            ImGui::Text("Standard error");
            ImGui::TableNextColumn();
            ImGui::Text("%.3f", updater.summary.se);
            ImGui::TableNextColumn();
            ImGui::Text("%.3f", updater.summary.se);
            ImGui::TableNextColumn();
            ImGui::Text("%.3f", updater.summary.se_kr);
            ImGui::TableNextColumn();

            ImGui::Text("Degrees of freedom");
            ImGui::TableNextColumn();
            ImGui::Text("%.1f", updater.summary.dof);
            ImGui::TableNextColumn();
            ImGui::Text("%.1f", updater.summary.dof_bw);
            ImGui::TableNextColumn();
            ImGui::Text("%.1f", updater.summary.dof_kr);
            ImGui::TableNextColumn();

            ImGui::EndTable();
        }

        if (option.two_treatments) {
            ImGui::Text("Treatment 2");
            if (ImGui::BeginTable("results 2", 4, flags))
            {
                ImGui::TableSetupColumn("Statistic");
                ImGui::TableSetupColumn("GLS");
                ImGui::TableSetupColumn("GLS (between-within)");
                ImGui::TableSetupColumn("Kenward-Roger");
                ImGui::TableHeadersRow();
                ImGui::TableNextColumn();

                ImGui::Text("Power (%%)");
                ImGui::TableNextColumn();
                ImGui::Text("%.1f", updater.summary.power_2);
                ImGui::TableNextColumn();
                ImGui::Text("%.1f", updater.summary.power_bw_2);
                ImGui::TableNextColumn();
                ImGui::Text("%.1f", updater.summary.power_kr_2);
                ImGui::TableNextColumn();

                ImGui::Text("Confidence interval half-width");
                ImGui::TableNextColumn();
                ImGui::Text("%.3f", updater.summary.ci_width_2);
                ImGui::TableNextColumn();
                ImGui::Text("%.3f", updater.summary.ci_width_bw_2);
                ImGui::TableNextColumn();
                ImGui::Text("%.3f", updater.summary.ci_width_kr_2);
                ImGui::TableNextColumn();

                ImGui::Text("Standard error");
                ImGui::TableNextColumn();
                ImGui::Text("%.3f", updater.summary.se_2);
                ImGui::TableNextColumn();
                ImGui::Text("%.3f", updater.summary.se_2);
                ImGui::TableNextColumn();
                ImGui::Text("%.3f", updater.summary.se_kr_2);
                ImGui::TableNextColumn();

                ImGui::Text("Degrees of freedom");
                ImGui::TableNextColumn();
                ImGui::Text("%.1f", updater.summary.dof_2);
                ImGui::TableNextColumn();
                ImGui::Text("%.1f", updater.summary.dof_bw_2);
                ImGui::TableNextColumn();
                ImGui::Text("%.1f", updater.summary.dof_kr_2);
                ImGui::TableNextColumn();

                ImGui::EndTable();
            }

            ImGui::Text("Interaction");
            if (ImGui::BeginTable("results interaction", 4, flags))
            {
                ImGui::TableSetupColumn("Statistic");
                ImGui::TableSetupColumn("GLS");
                ImGui::TableSetupColumn("GLS (between-within)");
                ImGui::TableSetupColumn("Kenward-Roger");
                ImGui::TableHeadersRow();
                ImGui::TableNextColumn();

                ImGui::Text("Power (%%)");
                ImGui::TableNextColumn();
                ImGui::Text("%.1f", updater.summary.power_12);
                ImGui::TableNextColumn();
                ImGui::Text("%.1f", updater.summary.power_bw_12);
                ImGui::TableNextColumn();
                ImGui::Text("%.1f", updater.summary.power_kr_12);
                ImGui::TableNextColumn();

                ImGui::Text("Confidence interval half-width");
                ImGui::TableNextColumn();
                ImGui::Text("%.3f", updater.summary.ci_width_12);
                ImGui::TableNextColumn();
                ImGui::Text("%.3f", updater.summary.ci_width_bw_12);
                ImGui::TableNextColumn();
                ImGui::Text("%.3f", updater.summary.ci_width_kr_12);
                ImGui::TableNextColumn();

                ImGui::Text("Standard error");
                ImGui::TableNextColumn();
                ImGui::Text("%.3f", updater.summary.se_12);
                ImGui::TableNextColumn();
                ImGui::Text("%.3f", updater.summary.se_12);
                ImGui::TableNextColumn();
                ImGui::Text("%.3f", updater.summary.se_kr_12);
                ImGui::TableNextColumn();

                ImGui::Text("Degrees of freedom");
                ImGui::TableNextColumn();
                ImGui::Text("%.1f", updater.summary.dof_12);
                ImGui::TableNextColumn();
                ImGui::Text("%.1f", updater.summary.dof_bw_12);
                ImGui::TableNextColumn();
                ImGui::Text("%.1f", updater.summary.dof_kr_12);
                ImGui::TableNextColumn();

                ImGui::EndTable();
            }
        }

        ImGui::End();
    }

    void ClusterApp::RenderOptimiser(ClusterApp::design& design, ClusterApp::modelUpdater& updater, ClusterApp::modelSummary& summary, ClusterApp::options& option) {
        ImGui::Begin("Optimiser", NULL, ImGuiWindowFlags_MenuBar);

        if (ImGui::BeginMenuBar()) {
            if (ImGui::BeginMenu("Options")) {
                if (ImGui::BeginMenu("Sample size")) {
                    ImGui::SetNextItemWidth(100);
                    ImGui::InputInt("Target total sample size", &summary.total_n, 1, 10, 0); ImGui::SameLine();
                    if (ImGui::SmallButton("Recalculate")) {
                        summary.total_n = design.total_n();
                        updater.update_optimum();
                    }
                    ImGui::EndMenu();
                }
                if (ImGui::BeginMenu("Parameter weights")) {
                    ImGui::SetNextItemWidth(100);
                    ImGui::DragFloat("Treatment 1", &updater.model.c_vals[0], 0.01f, -FLT_MAX, +FLT_MAX, "%.2f", ImGuiSliderFlags_None);
                    ImGui::SetNextItemWidth(100);
                    ImGui::DragFloat("Treatment 2", &updater.model.c_vals[1], 0.01f, -FLT_MAX, +FLT_MAX, "%.2f", ImGuiSliderFlags_None);
                    ImGui::SetNextItemWidth(100);
                    ImGui::DragFloat("Interaction", &updater.model.c_vals[2], 0.01f, -FLT_MAX, +FLT_MAX, "%.2f", ImGuiSliderFlags_None);
                    if (ImGui::SmallButton("Recalculate")) {
                        updater.update_optimum();
                    }
                    ImGui::EndMenu();
                }
                ImGui::EndMenu();
            }
            ImGui::EndMenuBar();
        }

        ImGui::TextWrapped("The design below shows the optimum weights per cluster in the design selected in the main window. Where sequences contain more than one cluster, these have been split out "); HelpMarker(
            "The colour of each cell indicates the weight. Click on the cell to see the weight and the rounded total number of observation. Rounding uses Hamilton's method. The design can be applied to the designer using the button below.");
  
        if (option.two_treatments) {
            ImGui::TextWrapped("For two treatments, the design minimises the variance of a weighted combination of the two treatment effect parameters and the interaction. You can change the weights in the options menu of this window.");
        }

        static ClusterApp::colourPicker colors;

        if (updater.optimum_data.size() != updater.designs.total_clusters()) {
            updater.update_data();
            updater.update_optimum();
        }

        int T = updater.optimum_data[0].size();
        int N = updater.optimum_data.size();
        int dim = 35;
        ImGuiStyle& style = ImGui::GetStyle();
        int horizontal_align = dim;

        ImGui::Dummy(ImVec2(dim*0.4, dim*0.4)); ImGui::SameLine(horizontal_align);
        for (int t = 0; t < T; t++) {
            ImGui::Text(int_to_char(t + 1));
            if (t < T - 1)ImGui::SameLine(horizontal_align + (t + 1) * (dim + style.ItemSpacing[0]));
        }

        for (int i = 0; i < N; i++) {            
            int seq = updater.designs.seq_by_cluster(i);
            ImGui::Text(int_to_char(i+1)); ImGui::SameLine(horizontal_align);
            for (int t = 0; t < T; t++) {
                float w = updater.optimum_data[i][t];
                ImVec4 newcol;
                if (option.light_mode) {
                    newcol = colors.base3();
                }
                else {
                    newcol = colors.base03();
                }
                bool both_two = option.two_treatments && *updater.designs.intervention(seq, t) && *updater.designs.intervention_2(seq, t);
                if (w > 1e-3) {
                    if (*updater.designs.intervention(seq, t) && !(option.two_treatments && *updater.designs.intervention_2(seq, t))) {
                        newcol = colors.red();
                        newcol.w *= 2*w;
                    }
                    else if(option.two_treatments && *updater.designs.intervention_2(seq, t)) {
                        newcol = colors.yellow();
                        newcol.w *= 2*w;
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
                int newid = t + i*T;
                ImGui::PushID(3*updater.designs.sequences*updater.designs.time + newid);
                if (!both_two) {
                    ImGui::Button("", ImVec2(dim, dim));
                    if (ImGui::BeginPopupContextItem(NULL, ImGuiPopupFlags_MouseButtonLeft))
                    {
                        ImGui::Text("Weight: "); ImGui::SameLine();
                        ImGui::Text("%.3f", w);
                        ImGui::Text("n: "); ImGui::SameLine();
                        ImGui::Text("%d", updater.optimum_n[i][t]);
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
                        ImGui::Text("%d", updater.optimum_n[i][t]);
                        ImGui::EndPopup();
                    }
                }
                                    
                ImGui::PopStyleColor(3);
                ImGui::PopID();
                if (t < updater.optimum_data[0].size() - 1)ImGui::SameLine();
            }
        }
        ImGui::Dummy(ImVec2(dim, dim));
        ImGui::Dummy(ImVec2(dim * 0.5, dim * 0.5)); ImGui::SameLine();
        if (ImGui::SmallButton("Apply Design")) {
            updater.designs.apply_design(updater.optimum_n);
        }

        ImGui::End();
    }
}

