#include "clusterapp.h"

namespace ClusterApp {

    void ShowMainMenu(ClusterApp::options& windows, ClusterApp::design& designs, ClusterApp::modelUpdater& updater, ClusterApp::modelSummary& summary) {
        if (ImGui::BeginMainMenuBar())
        {
            if (ImGui::BeginMenu("Trial Designs")) {
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
                    if (windows.two_treatments) {
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
                        ImGui::Text("To enable this option, set two treatments below.");
                    }
                    ImGui::EndMenu();
                }
                ImGui::Checkbox("Two treatments", &windows.two_treatments);
                ImGui::EndMenu();
            }
            if (ImGui::BeginMenu("Edit")) {

                if (ImGui::BeginMenu("Statistics")) {
                    ImGui::SetNextItemWidth(100);
                    ImGui::DragFloat("Alpha", &updater.model.alpha, 0.01f, 0.0f, 1.0f, "%.3f", ImGuiSliderFlags_None);

                    ImGui::EndMenu();
                }

                if (ImGui::BeginMenu("Design")) {
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
                        static int total_t = 1;
                        static int total_n = 10;
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

                if (ImGui::BeginMenu("Optimiser")) {
                    if (ImGui::BeginMenu("Sample size")) {
                        static int optim_total_n = designs.total_n();
                        ImGui::SetNextItemWidth(100);
                        ImGui::InputInt("Target total sample size", &optim_total_n, 1, 10, 0); ImGui::SameLine();
                        if (ImGui::SmallButton("Recalculate")) {
                            summary.total_n = optim_total_n;
                            updater.manual_n_optim = true;
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
                
                ImGui::EndMenu();
            }
            
            if (ImGui::BeginMenu("View"))
            {
                ImGui::Checkbox("Light mode", &windows.light_mode);
                ImGui::Checkbox("Sample size", &windows.sample_size);
                ImGui::Checkbox("Statistical model", &windows.model);
                ImGui::Checkbox("Results", &windows.results);
                ImGui::Checkbox("Optimal design", &windows.optimiser);
                ImGui::Checkbox("Plotting", &windows.plotter);
                ImGui::Checkbox("Dockspace", &windows.dockspace);
                ImGui::Checkbox("Debug Info", &windows.debug_info);
                ImGui::EndMenu();
            }
            if (ImGui::BeginMenu("Help")) {
                if (ImGui::Button("About"))
                    ImGui::OpenPopup("About");
                if (ImGui::BeginPopupModal("About", NULL, ImGuiWindowFlags_AlwaysAutoResize))
                {                    
                    ImGui::Text("(c) Sam Watson 2023.");
                    ImGui::Text("Version: 0.2.014");
                    ImGui::Text("glmmrBase Version: 0.4.6");
                    ImGui::Text("glmmrOptim Version: 0.3.1");
                    ImGui::Text("Code and license information is available on the GitHub repo.");

                    if (ImGui::Button("Close"))
                        ImGui::CloseCurrentPopup();
                    ImGui::EndPopup();
                }

                if (ImGui::Button("Changelog"))
                    ImGui::OpenPopup("Version info");
                if (ImGui::BeginPopupModal("Version info", NULL, ImGuiWindowFlags_AlwaysAutoResize))
                {
                    ImGui::Text("Version: 0.2.014");
                    ImGui::BulletText("Added plotting.");
                    ImGui::BulletText("Added set alpha value.");
                    ImGui::BulletText("Added parameters using group means.");
                    ImGui::BulletText("Few minor UI tweaks.");
                    ImGui::Text("Version: 0.2.012");
                    ImGui::BulletText("Added drag and drop for the designer.");
                    ImGui::Text("Version: 0.2.011");
                    ImGui::BulletText("Added design effects for exchangeable and nested exchangeable models.");
                    ImGui::Text("Version: 0.2.001");
                    ImGui::BulletText("Updated to glmmrBase 0.4.6.");
                    ImGui::BulletText("Added counts in designer for sequences and periods.");
                    ImGui::Text("Version: 0.1.123");
                    ImGui::BulletText("Set default start size and positions for windows.");
                    ImGui::BulletText("Updated cohort parameter for non-Gaussian models.");
                    ImGui::BulletText("Fixed optimum set sample size button not doing anything.");
                    ImGui::BulletText("Added placeholder for plotting.");
                    ImGui::Text("Version: 0.1.122");
                    ImGui::BulletText("Fixed aggregated model specification in call to glmmr - it wasn't dividing by number of observations.");
                    ImGui::Text("Version: 0.1.121");
                    ImGui::BulletText("Fixed covariance parameter values when using nested exchangeable function.");

                    if (ImGui::Button("Close"))
                        ImGui::CloseCurrentPopup();
                    ImGui::EndPopup();
                }
                ImGui::EndMenu();
            }

            ImGui::EndMainMenuBar();
        }
    }

    void AppDockSpace(bool* p_open)
    {

        static bool opt_fullscreen = true;
        static bool opt_padding = false;
        static ImGuiDockNodeFlags dockspace_flags = ImGuiDockNodeFlags_None;
        // We are using the ImGuiWindowFlags_NoDocking flag to make the parent window not dockable into,
        // because it would be confusing to have two docking targets within each others.
        ImGuiWindowFlags window_flags = ImGuiWindowFlags_MenuBar | ImGuiWindowFlags_NoDocking;
        if (opt_fullscreen)
        {
            const ImGuiViewport* viewport = ImGui::GetMainViewport();
            ImGui::SetNextWindowPos(viewport->WorkPos);
            ImGui::SetNextWindowSize(viewport->WorkSize);
            ImGui::SetNextWindowViewport(viewport->ID);
            ImGui::PushStyleVar(ImGuiStyleVar_WindowRounding, 0.0f);
            ImGui::PushStyleVar(ImGuiStyleVar_WindowBorderSize, 0.0f);
            window_flags |= ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoCollapse | ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoMove;
            window_flags |= ImGuiWindowFlags_NoBringToFrontOnFocus | ImGuiWindowFlags_NoNavFocus;
        }
        else
        {
            dockspace_flags &= ~ImGuiDockNodeFlags_PassthruCentralNode;
        }
        // When using ImGuiDockNodeFlags_PassthruCentralNode, DockSpace() will render our background
        // and handle the pass-thru hole, so we ask Begin() to not render a background.
        if (dockspace_flags & ImGuiDockNodeFlags_PassthruCentralNode)
            window_flags |= ImGuiWindowFlags_NoBackground;
        // Important: note that we proceed even if Begin() returns false (aka window is collapsed).
        // This is because we want to keep our DockSpace() active. If a DockSpace() is inactive,
        // all active windows docked into it will lose their parent and become undocked.
        // We cannot preserve the docking relationship between an active window and an inactive docking, otherwise
        // any change of dockspace/settings would lead to windows being stuck in limbo and never being visible.
        if (!opt_padding)
            ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding, ImVec2(0.0f, 0.0f));
        ImGui::Begin("DockSpace", p_open, window_flags);
        if (!opt_padding)
            ImGui::PopStyleVar();
        if (opt_fullscreen)
            ImGui::PopStyleVar(2);
        // Submit the DockSpace
        ImGuiIO& io = ImGui::GetIO();
        if (io.ConfigFlags & ImGuiConfigFlags_DockingEnable)
        {
            ImGuiID dockspace_id = ImGui::GetID("MainDockSpace");
            ImGui::DockSpace(dockspace_id, ImVec2(0.0f, 0.0f), dockspace_flags);
        }

        ImGui::End();
    }

    void RenderDesigner(ClusterApp::design& designs, ClusterApp::modelUpdater& updater, ClusterApp::options& option) {
        ImGui::Begin("Trial Designer");
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

        ImGui::TextWrapped("Design the trial below. Rows are sequences, columns are time periods. Click + to add new sequences or time periods. Select cells to edit their details or you can drag and drop them to new positions. Select row or column headers to change the numbers of clusters."); ImGui::SameLine(); HelpMarker(
            "You can change what the cell buttons show with the buttons below. Red and blue indicate intervention and control status, respectively. Where there are two treatments, yellow is used for treatment 2, and yellow-red for both treatments");
        ImGui::Text("For cluster-periods show:"); ImGui::SameLine();
        ImGui::Checkbox("Count (n)", &option.show_n_period); ImGui::SameLine();
        ImGui::Checkbox("Intervention status", &option.show_status_period);
        ImGui::Checkbox("Show number of cluster per sequence", &option.show_J_seq);
        ImGui::Text("Drag and drop mode: "); ImGui::SameLine();
        enum Mode
        {
            Mode_Copy,
            Mode_Move,
            Mode_Swap
        };
        static int mode = 0;
        if (ImGui::RadioButton("Copy", mode == Mode_Copy)) { mode = Mode_Copy; } ImGui::SameLine();
        if (ImGui::RadioButton("Move", mode == Mode_Move)) { mode = Mode_Move; } ImGui::SameLine();
        if (ImGui::RadioButton("Swap", mode == Mode_Swap)) { mode = Mode_Swap; }
        ImGui::Spacing();

        //ImGui::Dummy(ImVec2(large_dim, small_dim)); ImGui::SameLine(horiztonal_align + 30);
        //ImGui::Text("TIME");

        //ImGui::Text("SEQUENCE"); ImGui::SameLine(horiztonal_align);
        ImGui::Dummy(ImVec2(small_dim, small_dim)); ImGui::SameLine();
        ImGui::Dummy(ImVec2(small_dim, small_dim)); ImGui::SameLine();
        if (option.show_J_seq) {
            ImGui::Dummy(ImVec2(small_dim*1.5, small_dim)); ImGui::SameLine();
        }
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
        if (option.show_J_seq) {
            ImGui::Dummy(ImVec2(small_dim*1.5, small_dim)); ImGui::SameLine();
        }
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
            if (option.show_J_seq) {
                std::string label_j = std::to_string(*designs.n_clusters(n));
                char* char_array_j = new char[label_j.length() + 1];
                strcpy(char_array_j, label_j.c_str());
                ImGui::Button(char_array_j, ImVec2(small_dim * 1.5, small_dim)); ImGui::SameLine();
            }

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
                int id = t + n * designs.time;
                ImGui::PushID(id);
                std::string label = "";
                if (*designs.active(n, t)) {
                    if (option.show_n_period) {
                        label += "n=";
                        label += std::to_string(*designs.n(n, t));
                        if (option.show_status_period) label += " | ";
                    }
                    if (option.show_status_period) {
                        label += std::to_string(*designs.intervention(n, t));
                        if (option.two_treatments)label += "/" + std::to_string(*designs.intervention_2(n, t));
                    }

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

                // Our buttons are both drag sources and drag targets here!
                if (ImGui::BeginDragDropSource(ImGuiDragDropFlags_None))
                {
                    // Set payload to carry the index of our item (could be anything)
                    ImGui::SetDragDropPayload("DND_DEMO_CELL", &id, sizeof(int));

                    // Display preview (could be anything, e.g. when dragging an image we could decide to display
                    // the filename and a small preview of the image, etc.)
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
                            designs.copy_cells(payload_id, id);
                        }
                        if (mode == Mode_Move)
                        {
                            designs.move_cells(payload_id, id);
                        }
                        if (mode == Mode_Swap)
                        {
                            designs.swap_cells(payload_id, id);
                        }
                    }
                    ImGui::EndDragDropTarget();
                }

                ImGui::PopID();
                if (t < (designs.time - 1))ImGui::SameLine();
            }
        }

        ImGui::PushStyleColor(ImGuiCol_Button, colours.base1());
        ImGui::PushStyleColor(ImGuiCol_ButtonHovered, colours.base1(1));
        ImGui::PushStyleColor(ImGuiCol_ButtonActive, colours.base1(2));
        ImGui::PushStyleColor(ImGuiCol_Text, colours.base02(2));
        if (option.show_J_seq) {
            ImGui::Dummy(ImVec2(small_dim * 1.5, small_dim)); ImGui::SameLine();
        }

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

        if (option.debug_info) {
            ImGui::Text("Design Checksum"); ImGui::SameLine();
            ImGui::Text("%d", designs.crc_val);
        }
        
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

    void RenderModel(ClusterApp::design& design, ClusterApp::statisticalModel& model, ClusterApp::options& option) {
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
                    if (structure_sampling == 1) {
                        ImGui::SetNextItemWidth(200);
                        ImGui::DragFloat("IAC", &model.ixx_pars[2], 0.001f, 0.0f, 1.0f, "%.3f", ImGuiSliderFlags_None); ImGui::SameLine();  HelpMarker(
                            "Individual autocorrelation coefficient.");
                        if (ind_cov_item_current == 1) {
                            ImGui::SetNextItemWidth(200);
                            ImGui::DragFloat("Autoregressive (individual)", &model.cov_pars[4], 0.01f, 0.0f, 1.0f, "%.3f", ImGuiSliderFlags_None); ImGui::SameLine();  HelpMarker(
                                "Individual autoregressive parameter.");
                        }
                        if (ind_cov_item_current == 2 || ind_cov_item_current == 3) {
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
                    if (structure_sampling == 1) {
                        ImGui::SetNextItemWidth(200);
                        ImGui::DragFloat("Individual-level variance", &model.cov_pars[3], 0.01f, 0.0f, +FLT_MAX, "%.3f", ImGuiSliderFlags_None); ImGui::SameLine(); HelpMarker(
                            "Individual-level variance term");
                        if (ind_cov_item_current == 2) {
                            ImGui::SetNextItemWidth(200);
                            ImGui::DragFloat("Autoregressive (individual)", &model.cov_pars[4], 0.01f, 0.0f, 1.0f, "%.3f", ImGuiSliderFlags_None); ImGui::SameLine(); HelpMarker(
                                "Individual autoregressive parameter.");
                        }
                        if (ind_cov_item_current == 3 || ind_cov_item_current == 4) {
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
                            model.te_pars[0] = treatment_mean;
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
                                model.beta_pars[0] = log(control_mean/(1-control_mean)) - re_adj;
                            }
                            else {
                                for (int l = 0; l < model.beta_pars.size(); l++) {
                                    model.beta_pars[l] = log(control_mean / (1 - control_mean)) - re_adj;
                                }
                            }
                            model.te_pars[0] = log(treatment_mean/(1-treatment_mean)) - log(control_mean / (1 - control_mean));
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
                                    model.beta_pars[l] = boost::math::quantile(norm, control_mean) -re_adj;
                                }
                            }
                            model.te_pars[0] = boost::math::quantile(norm, treatment_mean) - model.beta_pars[0];
                            break; 
                        }

                        }
                    }
                                       

                    ImGui::TreePop();
                }
            }
            

            ImGui::TreePop();
        }
        if (option.debug_info) {
            if (ImGui::TreeNode("Model info")) {
                ImGui::Text("Model checksum"); ImGui::SameLine();
                ImGui::Text("%d", model.crc_val);
                ImGui::Text("Parameters checksum"); ImGui::SameLine();
                ImGui::Text("%d", model.crc_val_pars);
                ImGui::TreePop();
            }
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
            if (updater.summary.power_kr == 909) {
                ImGui::Text("ERR");
            }
            else {
                ImGui::Text("%.1f", updater.summary.power_kr);
            }            
            ImGui::TableNextColumn();

            ImGui::Text("Confidence interval half-width");
            ImGui::TableNextColumn();
            ImGui::Text("%.3f", updater.summary.ci_width);
            ImGui::TableNextColumn();
            ImGui::Text("%.3f", updater.summary.ci_width_bw);
            ImGui::TableNextColumn();
            if (updater.summary.power_kr == 909) {
                ImGui::Text("ERR");
            }
            else {
                ImGui::Text("%.1f", updater.summary.ci_width_kr);
            }
            ImGui::TableNextColumn();

            ImGui::Text("Standard error");
            ImGui::TableNextColumn();
            ImGui::Text("%.3f", updater.summary.se);
            ImGui::TableNextColumn();
            ImGui::Text("%.3f", updater.summary.se);
            ImGui::TableNextColumn();
            if (updater.summary.power_kr == 909) {
                ImGui::Text("ERR");
            }
            else {
                ImGui::Text("%.1f", updater.summary.se_kr);
            }
            ImGui::TableNextColumn();

            ImGui::Text("Degrees of freedom");
            ImGui::TableNextColumn();
            ImGui::Text("%.1f", updater.summary.dof);
            ImGui::TableNextColumn();
            ImGui::Text("%.1f", updater.summary.dof_bw);
            ImGui::TableNextColumn();
            if (updater.summary.power_kr == 909) {
                ImGui::Text("ERR");
            }
            else {
                ImGui::Text("%.1f", updater.summary.dof_kr);
            }
            ImGui::TableNextColumn();

            ImGui::EndTable();
        }

        if (!option.two_treatments && (updater.model.covariance == ClusterApp::Covariance::exchangeable || updater.model.covariance == ClusterApp::Covariance::nested_exchangeable)) {
            ImGui::Text("Design Effect Analysis"); ImGui::SameLine(); HelpMarker(
                "For comparison we report the power and confidence interval half-width for exchangeable and nested exchangeable models using an adapted approach of Hooper et al (2016). Use of ICC, CAC, and IAC values for non - Gaussian - identity models uses the mean individual - level variance, which is approximated using the GLM weights, to convert to covariance parameter values.");

            if (updater.summary.power_de == 909) {
                ImGui::TextWrapped("Error. Check parameter values. For example, if the log-binomial model produces probabilities outside [0,1] then this error will appear.");
            }
            else {
                ImGui::Text("Total observations: "); ImGui::SameLine();
                ImGui::Text("%.0f", updater.summary.individual_n); ImGui::SameLine(); HelpMarker(
                    "The sample size for the single period, individual-level RCT. This value is corrected for the allocation ratio between treatment and control conditions.");
                ImGui::Text("Design effect: "); ImGui::SameLine();
                ImGui::Text("%.3f", updater.summary.design_effect);
                ImGui::Text("Power: "); ImGui::SameLine();
                ImGui::Text("%.1f", updater.summary.power_de);
                ImGui::Text("95%% Confidence-interval half-width: "); ImGui::SameLine();
                ImGui::Text("%.3f", updater.summary.ci_width_de);
                ImGui::Text("Standard error: "); ImGui::SameLine();
                ImGui::Text("%.3f", updater.summary.se_de);
                ImGui::Text("Assumed observation level variance: "); ImGui::SameLine();
                ImGui::Text("%.3f", updater.summary.individual_var);
            }

            
            
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
                if (updater.summary.power_kr_2 == 909) {
                    ImGui::Text("ERR");
                }
                else {
                    ImGui::Text("%.1f", updater.summary.power_kr_2);
                }
                ImGui::TableNextColumn();

                ImGui::Text("Confidence interval half-width");
                ImGui::TableNextColumn();
                ImGui::Text("%.3f", updater.summary.ci_width_2);
                ImGui::TableNextColumn();
                ImGui::Text("%.3f", updater.summary.ci_width_bw_2);
                ImGui::TableNextColumn();
                if (updater.summary.power_kr_2 == 909) {
                    ImGui::Text("ERR");
                }
                else {
                    ImGui::Text("%.1f", updater.summary.ci_width_kr_2);
                }
                ImGui::TableNextColumn();

                ImGui::Text("Standard error");
                ImGui::TableNextColumn();
                ImGui::Text("%.3f", updater.summary.se_2);
                ImGui::TableNextColumn();
                ImGui::Text("%.3f", updater.summary.se_2);
                ImGui::TableNextColumn();
                if (updater.summary.power_kr_2 == 909) {
                    ImGui::Text("ERR");
                }
                else {
                    ImGui::Text("%.1f", updater.summary.se_kr_2);
                }
                ImGui::TableNextColumn();

                ImGui::Text("Degrees of freedom");
                ImGui::TableNextColumn();
                ImGui::Text("%.1f", updater.summary.dof_2);
                ImGui::TableNextColumn();
                ImGui::Text("%.1f", updater.summary.dof_bw_2);
                ImGui::TableNextColumn();
                if (updater.summary.power_kr_2 == 909) {
                    ImGui::Text("ERR");
                }
                else {
                    ImGui::Text("%.1f", updater.summary.dof_kr_2);
                }
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
                if (updater.summary.power_kr_12 == 909) {
                    ImGui::Text("ERR");
                }
                else {
                    ImGui::Text("%.1f", updater.summary.power_kr_12);
                }
                ImGui::TableNextColumn();

                ImGui::Text("Confidence interval half-width");
                ImGui::TableNextColumn();
                ImGui::Text("%.3f", updater.summary.ci_width_12);
                ImGui::TableNextColumn();
                ImGui::Text("%.3f", updater.summary.ci_width_bw_12);
                ImGui::TableNextColumn();
                if (updater.summary.power_kr_12 == 909) {
                    ImGui::Text("ERR");
                }
                else {
                    ImGui::Text("%.1f", updater.summary.ci_width_kr_12);
                }
                ImGui::TableNextColumn();

                ImGui::Text("Standard error");
                ImGui::TableNextColumn();
                ImGui::Text("%.3f", updater.summary.se_12);
                ImGui::TableNextColumn();
                ImGui::Text("%.3f", updater.summary.se_12);
                ImGui::TableNextColumn();
                if (updater.summary.power_kr_12 == 909) {
                    ImGui::Text("ERR");
                }
                else {
                    ImGui::Text("%.1f", updater.summary.se_kr_12);
                }
                ImGui::TableNextColumn();

                ImGui::Text("Degrees of freedom");
                ImGui::TableNextColumn();
                ImGui::Text("%.1f", updater.summary.dof_12);
                ImGui::TableNextColumn();
                ImGui::Text("%.1f", updater.summary.dof_bw_12);
                ImGui::TableNextColumn();
                if (updater.summary.power_kr_12 == 909) {
                    ImGui::Text("ERR");
                }
                else {
                    ImGui::Text("%.1f", updater.summary.dof_kr_12);
                }
                ImGui::TableNextColumn();

                ImGui::EndTable();
            }
        }

        ImGui::End();
    }

    void RenderOptimiser(ClusterApp::design& design, ClusterApp::modelUpdater& updater, ClusterApp::modelSummary& summary, ClusterApp::options& option) {
        ImGui::Begin("Optimiser");//, NULL, ImGuiWindowFlags_MenuBar


        ImGui::TextWrapped("The design below shows the optimum weights per cluster in the design selected in the main window. Where sequences contain more than one cluster, these have been split out "); ImGui::SameLine(); HelpMarker(
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

    void RenderPlotter(ClusterApp::plotData& plot) {
        ImGui::Begin("Plotter");
        ImGui::Text("Plot settings");

        const char* xaxis_items_exchangeable[] = { "Clusters", "N per cluster-period", "ICC", "Treatment effect", "Baseline"};
        const char* xaxis_items_other[] = { "Clusters", "N per cluster-period", "ICC", "Treatment effect", "Baseline", "CAC" };
        static int xaxis_item_current = 2;
        static int series_item_current = 4;
        const char* yaxis_items[] = { "Power (GLS)", "CI width (GLS)", "Power (GLS-BW)", "CI width (GLS-BW)", "Power (KR)", "CI width (KR)", "Power (Design effect)", "CI width (Design effect)" };
        static int yaxis_item_current = 0;
        static int series_clusters[] = { 1,2,3 };
        static int series_n[] = { 10,20,30 };
        static float series_icc[] = {0.01,0.02,0.05};
        static float series_te[] = {0.0, 0.5, 1.0};
        static float series_baseline[] = {0.0, 0.5, 1.0};
        static float series_cac[] = {0.2,0.5,0.8};
        const short s16_one = 1;
        static std::string series_label = "";
        colourPicker colours;
        static int print_prec = 3;
        std::string x_label = "";
        std::string ylabel = "";

        ImGui::SetNextItemWidth(200);
        if (plot.glmm.statmodel.covariance == ClusterApp::Covariance::exchangeable) {            
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
            plot.xaxis = ClusterApp::XAxis::clusters;
            ImGui::SetNextItemWidth(100);
            ImGui::InputScalar("Lower", ImGuiDataType_S16, &plot.lower_int[0], &s16_one, NULL, "%d"); ImGui::SameLine();
            ImGui::SetNextItemWidth(100);
            ImGui::InputScalar("Upper", ImGuiDataType_S16, &plot.upper_int[0], &s16_one, NULL, "%d");
            print_prec = 0;
            x_label = "Clusters";
            break;
        case 1:
            plot.xaxis = ClusterApp::XAxis::individual_n;
            ImGui::SetNextItemWidth(100);
            ImGui::InputScalar("Lower", ImGuiDataType_S16, &plot.lower_int[1], &s16_one, NULL, "%d"); ImGui::SameLine();
            ImGui::SetNextItemWidth(100);
            ImGui::InputScalar("Upper", ImGuiDataType_S16, &plot.upper_int[1], &s16_one, NULL, "%d");
            print_prec = 0;
            x_label = "n";
            break;
        case 2:
            plot.xaxis = ClusterApp::XAxis::icc;
            ImGui::SetNextItemWidth(100);
            ImGui::DragFloat("Lower", &plot.lower_float[0], 0.001f, 0.0f, 1.0f, "%.3f", ImGuiSliderFlags_None); ImGui::SameLine();
            ImGui::SetNextItemWidth(100);
            ImGui::DragFloat("Upper", &plot.upper_float[0], 0.001f, 0.0f, 1.0f, "%.3f", ImGuiSliderFlags_None);
            print_prec = 3;
            x_label = "ICC";
            break;
        case 3:
            plot.xaxis = ClusterApp::XAxis::treatment_effect;
            ImGui::SetNextItemWidth(100);
            ImGui::DragFloat("Lower", &plot.lower_float[1], 0.01f, -FLT_MAX, -FLT_MAX, "%.2f", ImGuiSliderFlags_None); ImGui::SameLine();
            ImGui::SetNextItemWidth(100);
            ImGui::DragFloat("Upper", &plot.upper_float[1], 0.01f, -FLT_MAX, -FLT_MAX, "%.2f", ImGuiSliderFlags_None);
            print_prec = 3;
            x_label = "Treatment effect";
            break;
        case 4:
            plot.xaxis = ClusterApp::XAxis::baseline;
            ImGui::SetNextItemWidth(100);
            ImGui::DragFloat("Lower", &plot.lower_float[2], 0.01f, -FLT_MAX, -FLT_MAX, "%.2f", ImGuiSliderFlags_None); ImGui::SameLine();
            ImGui::SetNextItemWidth(100);
            ImGui::DragFloat("Upper", &plot.upper_float[2], 0.01f, -FLT_MAX, -FLT_MAX, "%.2f", ImGuiSliderFlags_None);
            print_prec = 3;
            x_label = "Baseline";
            break;
        case 5:
            plot.xaxis = ClusterApp::XAxis::cac;
            ImGui::SetNextItemWidth(100);
            ImGui::DragFloat("Lower", &plot.lower_float[3], 0.001f, 0.0f, 1.0f, "%.3f", ImGuiSliderFlags_None); ImGui::SameLine();
            ImGui::SetNextItemWidth(100);
            ImGui::DragFloat("Upper", &plot.upper_float[3], 0.001f, 0.0f, 1.0f, "%.3f", ImGuiSliderFlags_None);
            print_prec = 3;
            x_label = "CAC";
            break;
        }               

        switch (yaxis_item_current) {
        case 0:
            plot.yaxis = ClusterApp::YAxis::power;
            ylabel = "Power (%%)";
            break;
        case 1:
            plot.yaxis = ClusterApp::YAxis::ci_width;
            ylabel = "95%% CI half-width";
            break;
        case 2:
            plot.yaxis = ClusterApp::YAxis::power_bw;
            ylabel = "Power (%%)";
            break;
        case 3:
            plot.yaxis = ClusterApp::YAxis::ci_width_bw;
            ylabel = "95%% CI half-width";
            break;
        case 4:
            plot.yaxis = ClusterApp::YAxis::power_kr;
            ylabel = "Power (%%)";
            break;
        case 5:
            plot.yaxis = ClusterApp::YAxis::ci_width_kr;
            ylabel = "95%% CI half-width";
            break;
        case 6:
            plot.yaxis = ClusterApp::YAxis::power_de;
            ylabel = "Power (%%)";
            break;
        case 7:
            plot.yaxis = ClusterApp::YAxis::ci_width_de;
            ylabel = "95%% CI half-width";
            break;
        }
 
        ImGui::Checkbox("Multiple series?", &plot.multiple_series);
        if (plot.multiple_series) {
            ImGui::Text("Add up to three series");
            ImGui::SetNextItemWidth(200);
            if (plot.glmm.statmodel.covariance == ClusterApp::Covariance::exchangeable) {
                ImGui::Combo("Series variable", &series_item_current, xaxis_items_exchangeable, IM_ARRAYSIZE(xaxis_items_exchangeable));
            }
            else {
                ImGui::Combo("Series variable", &series_item_current, xaxis_items_other, IM_ARRAYSIZE(xaxis_items_other));
            }
            switch (series_item_current) {
            case 0:
                plot.series = ClusterApp::XAxis::clusters;
                series_label = "Clusters";
                plot.x_series[0] = (float)series_clusters[0];
                plot.x_series[1] = (float)series_clusters[1];
                plot.x_series[2] = (float)series_clusters[2];
                break;
            case 1:
                plot.series = ClusterApp::XAxis::individual_n;
                series_label = "n";
                plot.x_series[0] = (float)series_n[0];
                plot.x_series[1] = (float)series_n[1];
                plot.x_series[2] = (float)series_n[2];
                break;
            case 2:
                plot.series = ClusterApp::XAxis::icc;
                series_label = "ICC";
                plot.x_series[0] = series_icc[0];
                plot.x_series[1] = series_icc[1];
                plot.x_series[2] = series_icc[2];
                break;
            case 3:
                plot.series = ClusterApp::XAxis::treatment_effect;
                series_label = "Effect size";
                plot.x_series[0] = series_te[0];
                plot.x_series[1] = series_te[1];
                plot.x_series[2] = series_te[2];
                break;
            case 4:
                plot.series = ClusterApp::XAxis::baseline;
                series_label = "Baseline";
                plot.x_series[0] = series_baseline[0];
                plot.x_series[1] = series_baseline[1];
                plot.x_series[2] = series_baseline[2];
                break;
            case 5:
                plot.series = ClusterApp::XAxis::cac;
                series_label = "CAC";
                plot.x_series[0] = series_cac[0];
                plot.x_series[1] = series_cac[1];
                plot.x_series[2] = series_cac[2];
                break;
            }

            //std::string label1 = std::to_string(plot.x_series[0]);

            std::stringstream stream;
            stream << std::fixed << std::setprecision(3) << plot.x_series[0];
            std::string label1 = stream.str();
            char* char_array = new char[label1.length() + 1];
            strcpy(char_array, label1.c_str());
            char* label_array = new char[series_label.length() + 1];
            strcpy(label_array, series_label.c_str());
            
            ImGui::PushID(10001);
            ImGui::PushStyleColor(ImGuiCol_Button, colours.red());
            ImGui::PushStyleColor(ImGuiCol_ButtonHovered, colours.red(1));
            ImGui::PushStyleColor(ImGuiCol_ButtonActive, colours.red(2));
            if(ImGui::Button(char_array, ImVec2(60, 40)))
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

            if (plot.n_series > 1) {
                std::stringstream stream2;
                stream2 << std::fixed << std::setprecision(3) << plot.x_series[1];
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

                if (plot.n_series == 3) {
                    std::stringstream stream3;
                    stream3 << std::fixed << std::setprecision(3) << plot.x_series[2];
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
                        plot.n_series = 2;
                    }
                    ImGui::PopID();
                }
                else {                   

                    ImGui::SameLine();
                    ImGui::PushID(10022);
                    if (ImGui::Button("-", ImVec2(20, 20))) {
                        plot.n_series = 1;
                    }
                    ImGui::PopID();

                    ImGui::SameLine();
                    ImGui::PushID(10032);
                    if (ImGui::Button("+", ImVec2(20, 20))) {
                        plot.n_series = 3;
                    }
                    ImGui::PopID();
                }
            }
            else {
                ImGui::PushID(10012);
                ImGui::SameLine();
                if (ImGui::Button("+", ImVec2(20, 20))) {
                    plot.n_series = 2;
                }
                ImGui::PopID();
            }

        }

        if (!plot.initialised)plot.update_data();

        char* x_char_array = new char[x_label.length() + 1];
        strcpy(x_char_array, x_label.c_str());
        char* y_char_array = new char[ylabel.length() + 1];
        strcpy(y_char_array, ylabel.c_str());

        static ImPlotAxisFlags xflags = ImPlotAxisFlags_AutoFit; //| ImPlotAxisFlags_RangeFit
        static ImPlotAxisFlags yflags = ImPlotAxisFlags_AutoFit;

        if (!plot.updating) {
            if (ImPlot::BeginPlot("Cluster trial plot")) {
                ImPlot::SetupAxis(ImAxis_X1, x_char_array, xflags);
                ImPlot::SetupAxis(ImAxis_Y1, y_char_array, yflags);
                ImPlot::PushStyleColor(ImPlotCol_Line, colours.red());
                ImPlot::PlotLine("Series 1", plot.x_data, plot.y_data_1, plot.n_data_points);
                ImPlot::PopStyleColor();
                if (plot.multiple_series) {
                    if (plot.n_series > 1) {
                        ImPlot::PushStyleColor(ImPlotCol_Line, colours.blue());
                        ImPlot::PlotLine("Series 2", plot.x_data, plot.y_data_2, plot.n_data_points);
                        ImPlot::PopStyleColor();
                        if (plot.n_series == 3) {
                            ImPlot::PushStyleColor(ImPlotCol_Line, colours.green());
                            ImPlot::PlotLine("Series 3", plot.x_data, plot.y_data_3, plot.n_data_points);
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

        ImGui::Text("X-axis limits:");
        ImGui::Text("%.3f", plot.x_axis_limits.first);
        ImGui::Text("%.3f", plot.x_axis_limits.second);

        ImGui::End();
    }
}

