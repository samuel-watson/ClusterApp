#include "clusterapp.h"

namespace ClusterApp {

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

        if (option.debug_info) {
            ImGui::Text("Total N: "); ImGui::SameLine();
            ImGui::Text("%d", summary.total_n);
            ImGui::Text("Covariance parameters: "); ImGui::SameLine();
            for (auto i: updater.glmm.model->model.covariance.parameters_) {
                ImGui::Text("%.3f ", i); ImGui::SameLine();
            }
            ImGui::Text("Var par: "); ImGui::SameLine();
            ImGui::Text("%.3f", updater.glmm.model->model.data.var_par);
            ImGui::Text("Weights: "); ImGui::SameLine();
            for (auto i : updater.glmm.optimal_weights) {
                ImGui::Text("%.3f  ", i); ImGui::SameLine();
            }
        }

        int T = updater.optimum_data[0].size();
        int N = updater.optimum_data.size();
        int dim = 35;
        ImGuiStyle& style = ImGui::GetStyle();
        int horizontal_align = dim;

        ImGui::Dummy(ImVec2(dim * 0.4, dim * 0.4)); ImGui::SameLine(horizontal_align);
        for (int t = 0; t < T; t++) {
            ImGui::Text("%i", t + 1);
            if (t < T - 1)ImGui::SameLine(horizontal_align + (t + 1) * (dim + style.ItemSpacing[0]));
        }

        for (int i = 0; i < N; i++) {
            int seq = updater.designs.seq_by_cluster(i);
            ImGui::Text("%i", i + 1); ImGui::SameLine(horizontal_align);
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
                        newcol.w *= 2 * w;
                    }
                    else if (option.two_treatments && *updater.designs.intervention_2(seq, t)) {
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
                ImGui::PushID(3 * updater.designs.sequences * updater.designs.time + newid);
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

