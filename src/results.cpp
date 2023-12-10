#include "clusterapp.h"

namespace ClusterApp {

    void RenderResults(ClusterApp::modelUpdater& updater, ClusterApp::options& option) {
        ImGui::Begin("Results");
        ImGui::Text("RESULTS");

        ImGui::TextWrapped("The tables below show the estimated power, confidence interval widths, and other statistics from a mixed model with various standard error estimators, \
including Kenward-Roger. For covariance functions non-linear in parameters, the improved Kenward-Roger correction is also provided. The second tab provides the design effect and \
resulting power estimates, currently only t-test and Chi-squared test bases for these calculations are available");
        ImGui::TextWrapped("Note that the Kenward-Roger correction (and its improved variant) can provide a poor approximation in very small sample sizes, including both small numbers of \
clusters, and small numbers of time periods with temporal covariance functions, especially exponential decay functions. It may also display ERR if the calculation returns a non-positive definite \
matrix, which can usually be resolved by fiddling with the parameter values and forcing it to rerun. See the log for more information.");

        if (updater.update) {
            ImGui::SameLine(); ImGui::Text("Recalculating...");
        }
        static ImGuiTableFlags flags = ImGuiTableFlags_Borders | ImGuiTableFlags_Resizable | ImGuiTableFlags_Reorderable | ImGuiTableFlags_Hideable | ImGuiTableFlags_NoHostExtendX;

        if (ImGui::TreeNode("Model-based results")) {
            if (ImGui::BeginTable("results", 6, flags))
            {
                ImGui::TableSetupColumn("Statistic");
                ImGui::TableSetupColumn("GLS");
                ImGui::TableSetupColumn("GLS (B-W)");
                ImGui::TableSetupColumn("Satterthwaite");
                ImGui::TableSetupColumn("Kenward-Roger");
                ImGui::TableSetupColumn("Kenward-Roger (improved)");
                //if(option.show_box && updater.model.family == Family::gaussian && updater.model.link == Link::identity)ImGui::TableSetupColumn("Box");
                ImGui::TableHeadersRow();
                ImGui::TableNextColumn();

                ImGui::Text("Power (%%)");
                ImGui::TableNextColumn();
                ImGui::Text("%.1f", updater.summary.power);
                ImGui::TableNextColumn();
                ImGui::Text("%.1f", updater.summary.power_bw);
                ImGui::TableNextColumn();
                if (updater.summary.power_sat == 909) {
                    ImGui::Text("ERR");
                }
                else {
                    ImGui::Text("%.1f", updater.summary.power_sat);
                }
                ImGui::TableNextColumn();
                if (updater.summary.power_kr == 909) {
                    ImGui::Text("ERR");
                }
                else {
                    ImGui::Text("%.1f", updater.summary.power_kr);
                }
                ImGui::TableNextColumn();
                if (updater.model.covariance != Covariance::exchangeable && updater.model.covariance != Covariance::nested_exchangeable) {
                    if (updater.summary.power_kr2 == 909) {
                        ImGui::Text("ERR");
                    }
                    else {
                        ImGui::Text("%.1f", updater.summary.power_kr2);
                    }
                    
                }
                else {
                    ImGui::Text("N/A");
                }
                ImGui::TableNextColumn();
                /*if (option.show_box && updater.model.family == Family::gaussian && updater.model.link == Link::identity) {
                    ImGui::Text("N/A");
                    ImGui::TableNextColumn();
                }*/

                ImGui::Text("95%% CI half-width");
                ImGui::TableNextColumn();
                ImGui::Text("%.3f", updater.summary.ci_width);
                ImGui::TableNextColumn();
                ImGui::Text("%.3f", updater.summary.ci_width_bw);
                ImGui::TableNextColumn();
                if (updater.summary.ci_width_sat == 909) {
                    ImGui::Text("ERR");
                }
                else {
                    ImGui::Text("%.3f", updater.summary.ci_width_sat);
                }
                ImGui::TableNextColumn();
                if (updater.summary.power_kr == 909) {
                    ImGui::Text("ERR");
                }
                else {
                    ImGui::Text("%.3f", updater.summary.ci_width_kr);
                }
                ImGui::TableNextColumn();
                if (updater.model.covariance != Covariance::exchangeable && updater.model.covariance != Covariance::nested_exchangeable) {
                    if (updater.summary.power_kr2 == 909) {
                        ImGui::Text("ERR");
                    }
                    else {
                        ImGui::Text("%.3f", updater.summary.ci_width_kr2);
                    }
                }
                else {
                    ImGui::Text("N/A");
                }
                ImGui::TableNextColumn();
                /*if (option.show_box && updater.model.family == Family::gaussian && updater.model.link == Link::identity) {
                    ImGui::Text("N/A");
                    ImGui::TableNextColumn();
                }*/

                ImGui::Text("Standard error");
                ImGui::TableNextColumn();
                ImGui::Text("%.3f", updater.summary.se);
                ImGui::TableNextColumn();
                ImGui::Text("%.3f", updater.summary.se);
                ImGui::TableNextColumn();
                ImGui::Text("%.3f", updater.summary.se);
                ImGui::TableNextColumn();
                if (updater.summary.power_kr == 909) {
                    ImGui::Text("ERR");
                }
                else {
                    ImGui::Text("%.3f", updater.summary.se_kr);
                }
                ImGui::TableNextColumn();
                if (updater.model.covariance != Covariance::exchangeable && updater.model.covariance != Covariance::nested_exchangeable) {
                    if (updater.summary.power_kr2 == 909) {
                        ImGui::Text("ERR");
                    }
                    else {
                        ImGui::Text("%.3f", updater.summary.se_kr2);
                    }
                }
                else {
                    ImGui::Text("N/A");
                }
                ImGui::TableNextColumn();
                /*if (option.show_box && updater.model.family == Family::gaussian && updater.model.link == Link::identity) {
                    ImGui::Text("N/A");
                    ImGui::TableNextColumn();
                }*/

                ImGui::Text("DoF");
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
                if (updater.summary.power_kr == 909) {
                    ImGui::Text("ERR");
                }
                else {
                    ImGui::Text("%.1f", updater.summary.dof_kr);
                }
                ImGui::TableNextColumn();
                if (updater.model.covariance != Covariance::exchangeable && updater.model.covariance != Covariance::nested_exchangeable) {
                    if (updater.summary.power_kr2 == 909) {
                        ImGui::Text("ERR");
                    }
                    else {
                        ImGui::Text("%.1f", updater.summary.dof_kr);
                    }
                }
                else {
                    ImGui::Text("N/A");
                }
                ImGui::TableNextColumn();
                /*if (option.show_box && updater.model.family == Family::gaussian && updater.model.link == Link::identity) {
                    ImGui::Text("%.3f", updater.summary.dof_box);
                    ImGui::TableNextColumn();
                }*/

                ImGui::EndTable();
            }



            if (option.two_treatments) {
                ImGui::Text("Treatment 2");
                if (ImGui::BeginTable("results 2", 6, flags))
                {
                    ImGui::TableSetupColumn("Statistic");
                    ImGui::TableSetupColumn("GLS");
                    ImGui::TableSetupColumn("GLS (B-W)");
                    ImGui::TableSetupColumn("Satterthwaite");
                    ImGui::TableSetupColumn("Kenward-Roger");
                    ImGui::TableSetupColumn("Kenward-Roger (improved)");
                    //if (option.show_box && updater.model.family == Family::gaussian && updater.model.link == Link::identity)ImGui::TableSetupColumn("Box");
                    ImGui::TableHeadersRow();
                    ImGui::TableNextColumn();

                    ImGui::Text("Power (%%)");
                    ImGui::TableNextColumn();
                    ImGui::Text("%.1f", updater.summary.power_2);
                    ImGui::TableNextColumn();
                    ImGui::Text("%.1f", updater.summary.power_bw_2);
                    ImGui::TableNextColumn();
                    if (updater.summary.power_sat_2 == 909) {
                        ImGui::Text("ERR");
                    }
                    else {
                        ImGui::Text("%.1f", updater.summary.power_sat_2);
                    }
                    ImGui::TableNextColumn();
                    if (updater.summary.power_kr_2 == 909) {
                        ImGui::Text("ERR");
                    }
                    else {
                        ImGui::Text("%.1f", updater.summary.power_kr_2);
                    }
                    ImGui::TableNextColumn();
                    if (updater.model.covariance != Covariance::exchangeable && updater.model.covariance != Covariance::nested_exchangeable) {
                        if (updater.summary.power_kr2_2 == 909) {
                            ImGui::Text("ERR");
                        }
                        else {
                            ImGui::Text("%.1f", updater.summary.power_kr2_2);
                        }
                    }
                    else {
                        ImGui::Text("N/A");
                    }
                    ImGui::TableNextColumn();
                    /*if (option.show_box && updater.model.family == Family::gaussian && updater.model.link == Link::identity) {
                        ImGui::Text("N/A");
                        ImGui::TableNextColumn();
                    }*/

                    ImGui::Text("95%% half-width");
                    ImGui::TableNextColumn();
                    ImGui::Text("%.3f", updater.summary.ci_width_2);
                    ImGui::TableNextColumn();
                    ImGui::Text("%.3f", updater.summary.ci_width_bw_2);
                    ImGui::TableNextColumn();
                    if (updater.summary.ci_width_sat_2 == 909) {
                        ImGui::Text("ERR");
                    }
                    else {
                        ImGui::Text("%.3f", updater.summary.ci_width_sat_2);
                    }
                    ImGui::TableNextColumn();
                    if (updater.summary.power_kr_2 == 909) {
                        ImGui::Text("ERR");
                    }
                    else {
                        ImGui::Text("%.3f", updater.summary.ci_width_kr_2);
                    }
                    ImGui::TableNextColumn();
                    if (updater.model.covariance != Covariance::exchangeable && updater.model.covariance != Covariance::nested_exchangeable) {
                        if (updater.summary.power_kr2_2 == 909) {
                            ImGui::Text("ERR");
                        }
                        else {
                            ImGui::Text("%.3f", updater.summary.ci_width_kr2_2);
                        }
                    }
                    else {
                        ImGui::Text("N/A");
                    }
                    ImGui::TableNextColumn();
                    /*if (option.show_box && updater.model.family == Family::gaussian && updater.model.link == Link::identity) {
                        ImGui::Text("N/A");
                        ImGui::TableNextColumn();
                    }*/

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
                        ImGui::Text("%.3f", updater.summary.se_kr_2);
                    }
                    ImGui::TableNextColumn();
                    if (updater.model.covariance != Covariance::exchangeable && updater.model.covariance != Covariance::nested_exchangeable) {
                        if (updater.summary.power_kr2_2 == 909) {
                            ImGui::Text("ERR");
                        }
                        else {
                            ImGui::Text("%.3f", updater.summary.se_kr2_2);
                        }
                    }
                    else {
                        ImGui::Text("N/A");
                    }
                    ImGui::TableNextColumn();
                    /*if (option.show_box && updater.model.family == Family::gaussian && updater.model.link == Link::identity) {
                        ImGui::Text("N/A");
                        ImGui::TableNextColumn();
                    }*/

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
                    if (updater.summary.power_kr_2 == 909) {
                        ImGui::Text("ERR");
                    }
                    else {
                        ImGui::Text("%.1f", updater.summary.dof_kr_2);
                    }
                    ImGui::TableNextColumn();
                    if (updater.model.covariance != Covariance::exchangeable && updater.model.covariance != Covariance::nested_exchangeable) {
                        if (updater.summary.power_kr2_2 == 909) {
                            ImGui::Text("ERR");
                        }
                        else {
                            ImGui::Text("%.3f", updater.summary.dof_kr_2);
                        }
                    }
                    else {
                        ImGui::Text("N/A");
                    }
                    ImGui::TableNextColumn();
                    /*if (option.show_box && updater.model.family == Family::gaussian && updater.model.link == Link::identity) {
                        ImGui::Text("%.3f", updater.summary.dof_box_2);
                        ImGui::TableNextColumn();
                    }*/

                    ImGui::EndTable();
                }

                ImGui::Text("Interaction");
                if (ImGui::BeginTable("results interaction", 6, flags))
                {
                    ImGui::TableSetupColumn("Statistic");
                    ImGui::TableSetupColumn("GLS");
                    ImGui::TableSetupColumn("GLS (B-W)");
                    ImGui::TableSetupColumn("Satterthwaite");
                    ImGui::TableSetupColumn("Kenward-Roger");
                    ImGui::TableSetupColumn("Kenward-Roger (improved)");
                    //if (option.show_box && updater.model.family == Family::gaussian && updater.model.link == Link::identity)ImGui::TableSetupColumn("Box");
                    ImGui::TableHeadersRow();
                    ImGui::TableNextColumn();

                    ImGui::Text("Power (%%)");
                    ImGui::TableNextColumn();
                    ImGui::Text("%.1f", updater.summary.power_12);
                    ImGui::TableNextColumn();
                    ImGui::Text("%.1f", updater.summary.power_bw_12);
                    ImGui::TableNextColumn();
                    if (updater.summary.power_sat_12 == 909) {
                        ImGui::Text("ERR");
                    }
                    else {
                        ImGui::Text("%.1f", updater.summary.power_sat_12);
                    }
                    ImGui::TableNextColumn();
                    if (updater.summary.power_kr_12 == 909) {
                        ImGui::Text("ERR");
                    }
                    else {
                        ImGui::Text("%.1f", updater.summary.power_kr_12);
                    }
                    ImGui::TableNextColumn();
                    if (updater.model.covariance != Covariance::exchangeable && updater.model.covariance != Covariance::nested_exchangeable) {
                        if (updater.summary.power_kr2_12 == 909) {
                            ImGui::Text("ERR");
                        }
                        else {
                            ImGui::Text("%.1f", updater.summary.power_kr2_12);
                        }
                    }
                    else {
                        ImGui::Text("N/A");
                    }
                    ImGui::TableNextColumn();
                    /*if (option.show_box && updater.model.family == Family::gaussian && updater.model.link == Link::identity) {
                        ImGui::Text("N/A");
                        ImGui::TableNextColumn();
                    }*/

                    ImGui::Text("95%% half-width");
                    ImGui::TableNextColumn();
                    ImGui::Text("%.3f", updater.summary.ci_width_12);
                    ImGui::TableNextColumn();
                    ImGui::Text("%.3f", updater.summary.ci_width_bw_12);
                    ImGui::TableNextColumn();
                    if (updater.summary.ci_width_sat_12 == 909) {
                        ImGui::Text("ERR");
                    }
                    else {
                        ImGui::Text("%.3f", updater.summary.ci_width_sat_12);
                    }
                    ImGui::TableNextColumn();
                    if (updater.summary.power_kr_12 == 909) {
                        ImGui::Text("ERR");
                    }
                    else {
                        ImGui::Text("%.3f", updater.summary.ci_width_kr_12);
                    }
                    ImGui::TableNextColumn();
                    if (updater.model.covariance != Covariance::exchangeable && updater.model.covariance != Covariance::nested_exchangeable) {
                        if (updater.summary.power_kr2_12 == 909) {
                            ImGui::Text("ERR");
                        }
                        else {
                            ImGui::Text("%.3f", updater.summary.ci_width_kr2_12);
                        }
                    }
                    else {
                        ImGui::Text("N/A");
                    }
                    ImGui::TableNextColumn();
                    /*if (option.show_box && updater.model.family == Family::gaussian && updater.model.link == Link::identity) {
                        ImGui::Text("N/A");
                        ImGui::TableNextColumn();
                    }*/

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
                        ImGui::Text("%.3f", updater.summary.se_kr_12);
                    }
                    ImGui::TableNextColumn();
                    if (updater.model.covariance != Covariance::exchangeable && updater.model.covariance != Covariance::nested_exchangeable) {
                        if (updater.summary.power_kr2_12 == 909) {
                            ImGui::Text("ERR");
                        }
                        else {
                            ImGui::Text("%.3f", updater.summary.se_kr2_12);
                        }
                    }
                    else {
                        ImGui::Text("N/A");
                    }
                    ImGui::TableNextColumn();
                   /* if (option.show_box && updater.model.family == Family::gaussian && updater.model.link == Link::identity) {
                        ImGui::Text("N/A");
                        ImGui::TableNextColumn();
                    }*/

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
                    if (updater.summary.power_kr_12 == 909) {
                        ImGui::Text("ERR");
                    }
                    else {
                        ImGui::Text("%.1f", updater.summary.dof_kr_12);
                    }
                    ImGui::TableNextColumn();
                    if (updater.model.covariance != Covariance::exchangeable && updater.model.covariance != Covariance::nested_exchangeable) {
                        if (updater.summary.power_kr2_12 == 909) {
                            ImGui::Text("ERR");
                        }
                        else {
                            ImGui::Text("%.3f", updater.summary.dof_kr_12);
                        }
                        ImGui::TableNextColumn();
                    }
                    else {
                        ImGui::Text("N/A");
                    }
                    ImGui::TableNextColumn();
                    /*if (option.show_box && updater.model.family == Family::gaussian && updater.model.link == Link::identity) {
                        ImGui::Text("%.3f", updater.summary.dof_box_12);
                        ImGui::TableNextColumn();
                    }*/

                    ImGui::EndTable();
                }
            }
            ImGui::TreePop();
        }

        enum Dmethod {
            Dmethod_glm,
            Dmethod_chisq
        };
        static int de_mode = 0;
        static int de_mode_curr_state = 0;

        if (!option.two_treatments && (updater.model.covariance == ClusterApp::Covariance::exchangeable || updater.model.covariance == ClusterApp::Covariance::nested_exchangeable)) {
            if (ImGui::TreeNode("Design Effects")) {
                ImGui::Text("Design Effect Analysis"); ImGui::SameLine(); HelpMarker(
                    "For comparison we report the power and confidence interval half-width for exchangeable and nested exchangeable models using an adapted approach of Hooper et al (2016). Use of ICC, CAC, and IAC values for non - Gaussian - identity models uses the mean individual - level variance, which is approximated using the GLM weights, to convert to covariance parameter values.");


                ImGui::Text("Design effect: "); ImGui::SameLine();
                ImGui::Text("%.3f", updater.summary.design_effect);
                ImGui::Text("Total observations: "); ImGui::SameLine();
                ImGui::Text("%.0f", updater.summary.individual_n); ImGui::SameLine(); HelpMarker(
                    "The sample size for the single period, individual-level RCT. This value is corrected for the allocation ratio between treatment and control conditions.");
                if (updater.model.family == ClusterApp::Family::binomial) {
                    if (ImGui::RadioButton("GLM variance", de_mode == Dmethod_glm)) { de_mode = Dmethod_glm; } ImGui::SameLine();
                    if (ImGui::RadioButton("Chi-squared test", de_mode == Dmethod_chisq)) { de_mode = Dmethod_chisq; }
                }
                else {
                    de_mode = Dmethod_glm;
                }
                updater.de_mode = de_mode;

                if (de_mode_curr_state != de_mode) {
                    updater.glmm.power_de(updater.summary, de_mode);
                    de_mode_curr_state = de_mode;
                }
                if (updater.summary.power_de == 909) {
                    ImGui::TextWrapped("Error. Check parameter values. For example, if the model produces probabilities outside [0,1] or if there are cells with zero probability for 2x2 test then this error will appear.");
                    updater.log.AddLog("[%05d] [%s] Design effect power error. Int group mean: %.3f, Control group mean: %.3f, ratio: %.3f \n", ImGui::GetFrameCount(), updater.log.cat[2], updater.summary.se_de, updater.summary.ci_width_de, updater.summary.individual_se);
                    /*if (option.debug_info) {
                        ImGui::Text("Control group mean: "); ImGui::SameLine();
                        ImGui::Text("%.3f", updater.summary.se_de);
                        ImGui::Text("Intervention group mean: "); ImGui::SameLine();
                        ImGui::Text("%.3f", updater.summary.ci_width_de);
                        ImGui::Text("Ratio: "); ImGui::SameLine();
                        ImGui::Text("%.3f", updater.summary.individual_se);
                    }*/
                }
                else {
                    if (de_mode == Dmethod_glm) {
                        ImGui::Text("Power: "); ImGui::SameLine();
                        ImGui::Text("%.1f", updater.summary.power_de);
                        ImGui::Text("95%% Confidence-interval half-width: "); ImGui::SameLine();
                        ImGui::Text("%.3f", updater.summary.ci_width_de);
                        ImGui::Text("Standard error: "); ImGui::SameLine();
                        ImGui::Text("%.3f", updater.summary.se_de);
                        ImGui::Text("Assumed observation level variance: "); ImGui::SameLine();
                        ImGui::Text("%.3f", updater.summary.individual_var);
                    }
                    else if (de_mode == Dmethod_chisq) {
                        ImGui::Text("Power: "); ImGui::SameLine();
                        ImGui::Text("%.1f", updater.summary.power_de);
                    }
                }

                ImGui::TreePop();
            }

        }




        ImGui::End();
    }

}

