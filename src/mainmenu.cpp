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
                        ImGui::Text("To enable this option, set two treatments.");
                    }
                    ImGui::EndMenu();
                }
                
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
                        ImGui::PushID(2000);
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
                        ImGui::PopID();
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
                        ImGui::PushID(2001);
                        if (ImGui::Button("Set")) {
                            if (total_t > 0) {
                                for (int j = 0; j < designs.sequences; j++)  *(designs.n_clusters(j)) = total_t;
                            }
                        }
                        ImGui::PopID();
                        ImGui::SetNextItemWidth(100);
                        ImGui::InputInt("n per cluster-period", &total_n, 1, 10, 0); ImGui::SameLine();
                        ImGui::PushID(2002);
                        if (ImGui::Button("Set n")) {
                            if (total_n > 0) {
                                for (int j = 0; j < designs.sequences; j++) {
                                    for (int t = 0; t < designs.time; t++) {
                                        *(designs.n(j, t)) = total_n;
                                    }
                                }
                            }
                        }
                        ImGui::PopID();
                        ImGui::EndMenu();
                    }
                    if (ImGui::MenuItem("Activate all")) {
                        for (int i = 0; i < designs.sequences; i++) {
                            for (int t = 0; t < designs.time; t++) {
                                *(designs.active(i, t)) = true;
                            }
                        }
                    }
                    if (ImGui::MenuItem("Deactivate all")) {
                        for (int i = 0; i < designs.sequences; i++) {
                            for (int t = 0; t < designs.time; t++) {
                                *(designs.active(i, t)) = false;
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

                if (ImGui::BeginMenu("Treatment effect")) {
                    if(!windows.dose_effect)ImGui::Checkbox("Two treatments", &windows.two_treatments);
                    if(!windows.two_treatments)ImGui::Checkbox("Dose effect", &windows.dose_effect);
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

            if (ImGui::BeginMenu("Options"))
            {
                if (ImGui::BeginMenu("Windows")) {
                    ImGui::Checkbox("Light mode", &windows.light_mode);
                    ImGui::Checkbox("Sample size", &windows.sample_size);
                    ImGui::Checkbox("Statistical model", &windows.model);
                    ImGui::Checkbox("Results", &windows.results);
                    ImGui::Checkbox("Optimal design", &windows.optimiser);
                    ImGui::Checkbox("Plotting", &windows.plotter);
                    ImGui::Checkbox("Sample minimiser", &windows.krigger);
                    ImGui::Checkbox("Data simulate", &windows.simulate);
                    ImGui::Checkbox("Dockspace", &windows.dockspace);
                    
                    //ImGui::Checkbox("Debug Info", &windows.debug_info);
                    ImGui::Checkbox("Show Box Correction", &windows.show_box);
                    ImGui::EndMenu();
                }
                if (ImGui::BeginMenu("Designer")) {
                    ImGui::Checkbox("Show cluster-period count (n)", &windows.show_n_period); 
                    ImGui::Checkbox("Show intervention status", &windows.show_status_period); ImGui::SameLine(); HelpMarker("This will show a 1/0 for intervention control or the dose for a dose response design.");
                    ImGui::Checkbox("Show number of clusters per sequence", &windows.show_J_seq);
                    ImGui::EndMenu();
                }
                if (ImGui::BeginMenu("Settings")) {
                    ImGui::Checkbox("Auto-update", &windows.auto_update);
                    ImGui::Checkbox("Console", &windows.log);
                    ImGui::EndMenu();
                }
                ImGui::EndMenu();
            }
            if (ImGui::BeginMenu("Help")) {
                if (ImGui::Button("About"))
                    ImGui::OpenPopup("About");
                if (ImGui::BeginPopupModal("About", NULL, ImGuiWindowFlags_AlwaysAutoResize))
                {
                    ImGui::Text("(c) Sam Watson 2023");
                    ImGui::Text("Version: 0.4.2");
                    ImGui::Text("glmmrBase Version: 0.8.1");
                    ImGui::Text("SparseChol Version: 0.2.2");
                    ImGui::Text("Code and license information is available on the GitHub repo.");

                    if (ImGui::Button("Close"))
                        ImGui::CloseCurrentPopup();
                    ImGui::EndPopup();
                }

                if (ImGui::Button("Changelog"))
                    ImGui::OpenPopup("Version info");
                if (ImGui::BeginPopupModal("Version info", NULL, ImGuiWindowFlags_AlwaysAutoResize))
                {
                    ImGui::Text("Version 0.4.2");
                    ImGui::BulletText("Updated to glmmrBase version 0.8.1");
                    ImGui::BulletText("Updating results is now manual by default, this provides significant improvement to performance when manipulating designs");
                    ImGui::BulletText("Hidden some text behind a tree to make the UI cleaner");
                    ImGui::Text("Version 0.4.1");
                    ImGui::BulletText("Updated to glmmrBase version 0.5.4");
                    ImGui::BulletText("Added Satterthwaite, KR Improved, and Box corrections (although Box not exposed)");
                    ImGui::BulletText("Added logging of events and console (see View->Window->Console)");
                    ImGui::BulletText("Added descriptive text to result and model windows");
                    ImGui::BulletText("Fixed an error with open cohort covariance parameter specification");
                    ImGui::Text("Version 0.3.5");
                    ImGui::BulletText("Moved design view options to menu bar");
                    ImGui::Text("Version 0.3.4");
                    ImGui::BulletText("Added dose response model");
                    ImGui::Text("Version: 0.3.3");
                    ImGui::BulletText("Fixed integer rounding in calculating plot increments");
                    ImGui::BulletText("Improved designer UI for modifying sample sizes");
                    ImGui::BulletText("Tidied up some compiler warnings");
                    ImGui::Text("Version: 0.3.2");
                    ImGui::BulletText("Added sample minimiser (experimental).");
                    ImGui::Text("Version: 0.3.1");
                    ImGui::BulletText("Added open cohort designs.");
                    ImGui::BulletText("Added data simulator.");
                    ImGui::BulletText("Simplified version numbering!");
                    ImGui::Text("Version: 0.2.015/2.016");
                    ImGui::BulletText("Bug fixes and UI improvements.");
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

    
}

