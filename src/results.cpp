#include "clusterapp.h"

namespace ClusterApp {

    void RenderResults(ClusterApp::modelChecker& checker, ClusterApp::options& option) {
        ImGui::Begin("Results");
        ImGui::Text("RESULTS");
        static colourPicker colours;
        if (ImGui::TreeNode("Info")){
            ImGui::TextWrapped("The tables below show the estimated power, confidence interval widths, and other statistics from a mixed model with various standard error estimators, \
including Kenward-Roger. For covariance functions non-linear in parameters, the improved Kenward-Roger correction is also provided. The second tab provides the design effect and \
resulting power estimates, currently only t-test and Chi-squared test bases for these calculations are available");
        ImGui::TextWrapped("Note that the Kenward-Roger correction (and its improved variant) can provide a poor approximation in very small sample sizes, including both small numbers of \
clusters, and small numbers of time periods with temporal covariance functions, especially exponential decay functions. It may also display ERR if the calculation returns a non-positive definite \
matrix, which can usually be resolved by fiddling with the parameter values and forcing it to rerun. See the log for more information.");
        ImGui::TreePop();
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
        } else {
            if (checker.updater.update) {
                ImGui::SameLine(); ImGui::Text("Recalculating...");
            }
        static ImGuiTableFlags flags = ImGuiTableFlags_Borders | ImGuiTableFlags_Resizable | ImGuiTableFlags_Reorderable | ImGuiTableFlags_Hideable | ImGuiTableFlags_NoHostExtendX;

        // make this into collapsible headers
        if(ImGui::CollapsingHeader("Power (%)")){

            ImGui::BulletText("Mixed model: "); ImGui::SameLine();
            ImGui::Text("%.1f", checker.updater.summary.power);
            ImGui::BulletText("Mixed model t-test (between-within): "); ImGui::SameLine();
            ImGui::Text("%.1f", checker.updater.summary.power_bw);
            ImGui::BulletText("Satterthwaite: "); ImGui::SameLine();
            if (checker.updater.summary.power_sat == 909) {
                ImGui::Text("ERR");
            }
            else {
                ImGui::Text("%.1f", checker.updater.summary.power_sat);
            }
            ImGui::BulletText("Kenward-Roger: "); ImGui::SameLine();
            if (checker.updater.summary.power_kr == 909) {
                ImGui::Text("ERR");
            }
            else {
                ImGui::Text("%.1f", checker.updater.summary.power_kr);
            }
            ImGui::BulletText("GEE Independence working correlation, robust: "); ImGui::SameLine();
            if (checker.updater.summary.power_gee_indep == 909) {
                ImGui::Text("ERR");
            }
            else {
                ImGui::Text("%.1f", checker.updater.summary.power_gee_indep);
            }
            ImGui::BulletText("GEE Independence working correlation, robust, t-test: "); ImGui::SameLine();
            if (checker.updater.summary.power_gee_indep == 909) {
                ImGui::Text("ERR");
            }
            else {
                ImGui::Text("%.1f", checker.updater.summary.power_gee_indep_bw);
            }
            if (checker.updater.model.covariance != Covariance::exchangeable && checker.updater.model.covariance != Covariance::nested_exchangeable){
                ImGui::BulletText("Kenward-Roger (improved): "); ImGui::SameLine();
                if (checker.updater.summary.power_kr2 == 909) {
                    ImGui::Text("ERR");
                }
                else {
                    ImGui::Text("%.1f", checker.updater.summary.power_kr2);
                }
            }

            if(option.two_treatments){
                ImGui::Text("Treatment 2");
                ImGui::BulletText("Mixed model: "); ImGui::SameLine();
                ImGui::Text("%.1f", checker.updater.summary.power_2);
                ImGui::BulletText("Mixed model t-test (between-within): "); ImGui::SameLine();
                ImGui::Text("%.1f", checker.updater.summary.power_bw_2);
                ImGui::BulletText("Satterthwaite: "); ImGui::SameLine();
                if (checker.updater.summary.power_sat_2 == 909) {
                    ImGui::Text("ERR");
                }
                else {
                    ImGui::Text("%.1f", checker.updater.summary.power_sat_2);
                }
                ImGui::BulletText("Kenward-Roger: "); ImGui::SameLine();
                if (checker.updater.summary.power_kr_2 == 909) {
                    ImGui::Text("ERR");
                }
                else {
                    ImGui::Text("%.1f", checker.updater.summary.power_kr_2);
                }
                ImGui::BulletText("GEE Independence working correlation, robust: "); ImGui::SameLine();
                if (checker.updater.summary.power_gee_indep_2 == 909) {
                    ImGui::Text("ERR");
                }
                else {
                    ImGui::Text("%.1f", checker.updater.summary.power_gee_indep_2);
                }
                ImGui::BulletText("GEE Independence working correlation, robust, t-test: "); ImGui::SameLine();
                if (checker.updater.summary.power_gee_indep_2 == 909) {
                    ImGui::Text("ERR");
                }
                else {
                    ImGui::Text("%.1f", checker.updater.summary.power_gee_indep_bw_2);
                }
                if (checker.updater.model.covariance != Covariance::exchangeable && checker.updater.model.covariance != Covariance::nested_exchangeable){
                    ImGui::BulletText("Kenward-Roger (improved): "); ImGui::SameLine();
                    if (checker.updater.summary.power_kr2_2 == 909) {
                        ImGui::Text("ERR");
                    }
                    else {
                        ImGui::Text("%.1f", checker.updater.summary.power_kr2_2);
                    }
                }

                ImGui::Text("Treatment Interaction");
                ImGui::BulletText("Mixed model: "); ImGui::SameLine();
                ImGui::Text("%.1f", checker.updater.summary.power_12);
                ImGui::BulletText("Mixed model t-test (between-within): "); ImGui::SameLine();
                ImGui::Text("%.1f", checker.updater.summary.power_bw_12);
                ImGui::BulletText("Satterthwaite: "); ImGui::SameLine();
                if (checker.updater.summary.power_sat_12 == 909) {
                    ImGui::Text("ERR");
                }
                else {
                    ImGui::Text("%.1f", checker.updater.summary.power_sat_12);
                }
                ImGui::BulletText("Kenward-Roger: "); ImGui::SameLine();
                if (checker.updater.summary.power_kr_12 == 909) {
                    ImGui::Text("ERR");
                }
                else {
                    ImGui::Text("%.1f", checker.updater.summary.power_kr_12);
                }
                ImGui::BulletText("GEE Independence working correlation, robust: "); ImGui::SameLine();
                if (checker.updater.summary.power_gee_indep_12 == 909) {
                    ImGui::Text("ERR");
                }
                else {
                    ImGui::Text("%.1f", checker.updater.summary.power_gee_indep_12);
                }
                ImGui::BulletText("GEE Independence working correlation, robust, t-test: "); ImGui::SameLine();
                if (checker.updater.summary.power_gee_indep_12 == 909) {
                    ImGui::Text("ERR");
                }
                else {
                    ImGui::Text("%.1f", checker.updater.summary.power_gee_indep_bw_12);
                }
                if (checker.updater.model.covariance != Covariance::exchangeable && checker.updater.model.covariance != Covariance::nested_exchangeable){
                    ImGui::BulletText("Kenward-Roger (improved): "); ImGui::SameLine();
                    if (checker.updater.summary.power_kr2_12 == 909) {
                        ImGui::Text("ERR");
                    }
                    else {
                        ImGui::Text("%.1f", checker.updater.summary.power_kr2_12);
                    }
                }
            }
        }


        if(ImGui::CollapsingHeader("Minimum detectable effect size")){

            ImGui::BulletText("Mixed model: "); ImGui::SameLine();
            ImGui::Text("%.2f", checker.updater.summary.ci_width);
            ImGui::BulletText("Mixed model t-test (between-within): "); ImGui::SameLine();
            ImGui::Text("%.2f", checker.updater.summary.ci_width_bw);
            ImGui::BulletText("Satterthwaite: "); ImGui::SameLine();
            if (checker.updater.summary.power_sat == 909) {
                ImGui::Text("ERR");
            }
            else {
                ImGui::Text("%.2f", checker.updater.summary.ci_width_sat);
            }
            ImGui::BulletText("Kenward-Roger: "); ImGui::SameLine();
            if (checker.updater.summary.power_kr == 909) {
                ImGui::Text("ERR");
            }
            else {
                ImGui::Text("%.2f", checker.updater.summary.ci_width_kr);
            }
            ImGui::BulletText("GEE Independence working correlation, robust: "); ImGui::SameLine();
            if (checker.updater.summary.power_gee_indep == 909) {
                ImGui::Text("ERR");
            }
            else {
                ImGui::Text("%.2f", checker.updater.summary.ci_width_gee_indep);
            }
            ImGui::BulletText("GEE Independence working correlation, robust, t-test: "); ImGui::SameLine();
            if (checker.updater.summary.power_gee_indep == 909) {
                ImGui::Text("ERR");
            }
            else {
                ImGui::Text("%.2f", checker.updater.summary.ci_width_gee_indep_bw);
            }
            if (checker.updater.model.covariance != Covariance::exchangeable && checker.updater.model.covariance != Covariance::nested_exchangeable){
                ImGui::BulletText("Kenward-Roger (improved): "); ImGui::SameLine();
                if (checker.updater.summary.power_kr2 == 909) {
                    ImGui::Text("ERR");
                }
                else {
                    ImGui::Text("%.2f", checker.updater.summary.ci_width_kr2);
                }
            }

            if(option.two_treatments){
                ImGui::Text("Treatment 2");
                ImGui::BulletText("Mixed model: "); ImGui::SameLine();
                ImGui::Text("%.2f", checker.updater.summary.ci_width_2);
                ImGui::BulletText("Mixed model t-test (between-within): "); ImGui::SameLine();
                ImGui::Text("%.2f", checker.updater.summary.ci_width_bw_2);
                ImGui::BulletText("Satterthwaite: "); ImGui::SameLine();
                if (checker.updater.summary.power_sat_2 == 909) {
                    ImGui::Text("ERR");
                }
                else {
                    ImGui::Text("%.2f", checker.updater.summary.ci_width_sat_2);
                }
                ImGui::BulletText("Kenward-Roger: "); ImGui::SameLine();
                if (checker.updater.summary.power_kr_2 == 909) {
                    ImGui::Text("ERR");
                }
                else {
                    ImGui::Text("%.2f", checker.updater.summary.ci_width_kr_2);
                }
                ImGui::BulletText("GEE Independence working correlation, robust: "); ImGui::SameLine();
                if (checker.updater.summary.power_gee_indep_2 == 909) {
                    ImGui::Text("ERR");
                }
                else {
                    ImGui::Text("%.2f", checker.updater.summary.ci_width_gee_indep_2);
                }
                ImGui::BulletText("GEE Independence working correlation, robust, t-test: "); ImGui::SameLine();
                if (checker.updater.summary.power_gee_indep_2 == 909) {
                    ImGui::Text("ERR");
                }
                else {
                    ImGui::Text("%.2f", checker.updater.summary.ci_width_gee_indep_bw_2);
                }
                if (checker.updater.model.covariance != Covariance::exchangeable && checker.updater.model.covariance != Covariance::nested_exchangeable){
                    ImGui::BulletText("Kenward-Roger (improved): "); ImGui::SameLine();
                    if (checker.updater.summary.power_kr2_2 == 909) {
                        ImGui::Text("ERR");
                    }
                    else {
                        ImGui::Text("%.2f", checker.updater.summary.ci_width_kr2_2);
                    }
                }

                ImGui::Text("Treatment Interaction");
                ImGui::BulletText("Mixed model: "); ImGui::SameLine();
                ImGui::Text("%.2f", checker.updater.summary.ci_width_12);
                ImGui::BulletText("Mixed model t-test (between-within): "); ImGui::SameLine();
                ImGui::Text("%.2f", checker.updater.summary.ci_width_bw_12);
                ImGui::BulletText("Satterthwaite: "); ImGui::SameLine();
                if (checker.updater.summary.power_sat_12 == 909) {
                    ImGui::Text("ERR");
                }
                else {
                    ImGui::Text("%.2f", checker.updater.summary.ci_width_sat_12);
                }
                ImGui::BulletText("Kenward-Roger: "); ImGui::SameLine();
                if (checker.updater.summary.power_kr_12 == 909) {
                    ImGui::Text("ERR");
                }
                else {
                    ImGui::Text("%.2f", checker.updater.summary.ci_width_kr_12);
                }
                ImGui::BulletText("GEE Independence working correlation, robust: "); ImGui::SameLine();
                if (checker.updater.summary.power_gee_indep_12 == 909) {
                    ImGui::Text("ERR");
                }
                else {
                    ImGui::Text("%.2f", checker.updater.summary.ci_width_gee_indep_12);
                }
                ImGui::BulletText("GEE Independence working correlation, robust, t-test: "); ImGui::SameLine();
                if (checker.updater.summary.power_gee_indep_12 == 909) {
                    ImGui::Text("ERR");
                }
                else {
                    ImGui::Text("%.2f", checker.updater.summary.ci_width_gee_indep_bw_12);
                }
                if (checker.updater.model.covariance != Covariance::exchangeable && checker.updater.model.covariance != Covariance::nested_exchangeable){
                    ImGui::BulletText("Kenward-Roger (improved): "); ImGui::SameLine();
                    if (checker.updater.summary.power_kr2_12 == 909) {
                        ImGui::Text("ERR");
                    }
                    else {
                        ImGui::Text("%.2f", checker.updater.summary.ci_width_kr2_12);
                    }
                }
            }
        }

        if(ImGui::CollapsingHeader("Standard errors")){

            ImGui::BulletText("Mixed model: "); ImGui::SameLine();
            ImGui::Text("%.3f", checker.updater.summary.se);
            ImGui::BulletText("Mixed model t-test (between-within): "); ImGui::SameLine();
            ImGui::Text("%.3f", checker.updater.summary.se);
            ImGui::BulletText("Satterthwaite: "); ImGui::SameLine();
            if (checker.updater.summary.power_sat == 909) {
                ImGui::Text("ERR");
            }
            else {
                ImGui::Text("%.3f", checker.updater.summary.se);
            }
            ImGui::BulletText("Kenward-Roger: "); ImGui::SameLine();
            if (checker.updater.summary.power_kr == 909) {
                ImGui::Text("ERR");
            }
            else {
                ImGui::Text("%.3f", checker.updater.summary.se_kr);
            }
            ImGui::BulletText("GEE Independence working correlation, robust: "); ImGui::SameLine();
            if (checker.updater.summary.power_gee_indep == 909) {
                ImGui::Text("ERR");
            }
            else {
                ImGui::Text("%.3f", checker.updater.summary.se_gee_indep);
            }
            if (checker.updater.model.covariance != Covariance::exchangeable && checker.updater.model.covariance != Covariance::nested_exchangeable){
                ImGui::BulletText("Kenward-Roger (improved): "); ImGui::SameLine();
                if (checker.updater.summary.power_kr2 == 909) {
                    ImGui::Text("ERR");
                }
                else {
                    ImGui::Text("%.3f", checker.updater.summary.se_kr2);
                }
            }

            if(option.two_treatments){
                ImGui::Text("Treatment 2");
                ImGui::BulletText("Mixed model: "); ImGui::SameLine();
                ImGui::Text("%.3f", checker.updater.summary.se_2);
                ImGui::BulletText("Mixed model t-test (between-within): "); ImGui::SameLine();
                ImGui::Text("%.3f", checker.updater.summary.se_2);
                ImGui::BulletText("Satterthwaite: "); ImGui::SameLine();
                if (checker.updater.summary.power_sat_2 == 909) {
                    ImGui::Text("ERR");
                }
                else {
                    ImGui::Text("%.3f", checker.updater.summary.se_2);
                }
                ImGui::BulletText("Kenward-Roger: "); ImGui::SameLine();
                if (checker.updater.summary.power_kr_2 == 909) {
                    ImGui::Text("ERR");
                }
                else {
                    ImGui::Text("%.3f", checker.updater.summary.se_kr_2);
                }
                ImGui::BulletText("GEE Independence working correlation, robust: "); ImGui::SameLine();
                if (checker.updater.summary.power_gee_indep_2 == 909) {
                    ImGui::Text("ERR");
                }
                else {
                    ImGui::Text("%.3f", checker.updater.summary.se_gee_indep_2);
                }
                if (checker.updater.model.covariance != Covariance::exchangeable && checker.updater.model.covariance != Covariance::nested_exchangeable){
                    ImGui::BulletText("Kenward-Roger (improved): "); ImGui::SameLine();
                    if (checker.updater.summary.power_kr2_2 == 909) {
                        ImGui::Text("ERR");
                    }
                    else {
                        ImGui::Text("%.3f", checker.updater.summary.se_kr2_2);
                    }
                }

                ImGui::Text("Treatment Interaction");
                ImGui::BulletText("Mixed model: "); ImGui::SameLine();
                ImGui::Text("%.3f", checker.updater.summary.se_12);
                ImGui::BulletText("Mixed model t-test (between-within): "); ImGui::SameLine();
                ImGui::Text("%.3f", checker.updater.summary.se_12);
                ImGui::BulletText("Satterthwaite: "); ImGui::SameLine();
                if (checker.updater.summary.power_sat_12 == 909) {
                    ImGui::Text("ERR");
                }
                else {
                    ImGui::Text("%.3f", checker.updater.summary.se_12);
                }
                ImGui::BulletText("Kenward-Roger: "); ImGui::SameLine();
                if (checker.updater.summary.power_kr_12 == 909) {
                    ImGui::Text("ERR");
                }
                else {
                    ImGui::Text("%.3f", checker.updater.summary.se_kr_12);
                }
                ImGui::BulletText("GEE Independence working correlation, robust: "); ImGui::SameLine();
                if (checker.updater.summary.power_gee_indep_12 == 909) {
                    ImGui::Text("ERR");
                }
                else {
                    ImGui::Text("%.3f", checker.updater.summary.se_gee_indep_12);
                }
                if (checker.updater.model.covariance != Covariance::exchangeable && checker.updater.model.covariance != Covariance::nested_exchangeable){
                    ImGui::BulletText("Kenward-Roger (improved): "); ImGui::SameLine();
                    if (checker.updater.summary.power_kr2_12 == 909) {
                        ImGui::Text("ERR");
                    }
                    else {
                        ImGui::Text("%.3f", checker.updater.summary.se_kr2_12);
                    }
                }
            }
        }

        if(ImGui::CollapsingHeader("Degrees of freedom")){

            ImGui::BulletText("Mixed model: "); ImGui::SameLine();
            ImGui::Text("%.1f", checker.updater.summary.dof);
            ImGui::BulletText("Mixed model t-test (between-within): "); ImGui::SameLine();
            ImGui::Text("%.1f", checker.updater.summary.dof_bw);
            ImGui::BulletText("Satterthwaite: "); ImGui::SameLine();
            if (checker.updater.summary.power_sat == 909) {
                ImGui::Text("ERR");
            }
            else {
                ImGui::Text("%.1f", checker.updater.summary.dof_kr);
            }
            ImGui::BulletText("Kenward-Roger: "); ImGui::SameLine();
            if (checker.updater.summary.power_kr == 909) {
                ImGui::Text("ERR");
            }
            else {
                ImGui::Text("%.1f", checker.updater.summary.dof_kr);
            }
            ImGui::BulletText("GEE Independence working correlation, robust: "); ImGui::SameLine();
            if (checker.updater.summary.power_gee_indep == 909) {
                ImGui::Text("ERR");
            }
            else {
                ImGui::Text("%.1f", checker.updater.summary.dof);
            }
            if (checker.updater.model.covariance != Covariance::exchangeable && checker.updater.model.covariance != Covariance::nested_exchangeable){
                ImGui::BulletText("Kenward-Roger (improved): "); ImGui::SameLine();
                if (checker.updater.summary.power_kr2 == 909) {
                    ImGui::Text("ERR");
                }
                else {
                    ImGui::Text("%.1f", checker.updater.summary.dof_kr);
                }
            }

            if(option.two_treatments){
                ImGui::Text("Treatment 2");
                ImGui::BulletText("Mixed model: "); ImGui::SameLine();
                ImGui::Text("%.1f", checker.updater.summary.dof_2);
                ImGui::BulletText("Mixed model t-test (between-within): "); ImGui::SameLine();
                ImGui::Text("%.1f", checker.updater.summary.dof_bw_2);
                ImGui::BulletText("Satterthwaite: "); ImGui::SameLine();
                if (checker.updater.summary.power_sat_2 == 909) {
                    ImGui::Text("ERR");
                }
                else {
                    ImGui::Text("%.1f", checker.updater.summary.dof_kr_2);
                }
                ImGui::BulletText("Kenward-Roger: "); ImGui::SameLine();
                if (checker.updater.summary.power_kr_2 == 909) {
                    ImGui::Text("ERR");
                }
                else {
                    ImGui::Text("%.1f", checker.updater.summary.dof_kr_2);
                }
                ImGui::BulletText("GEE Independence working correlation, robust: "); ImGui::SameLine();
                if (checker.updater.summary.power_gee_indep_2 == 909) {
                    ImGui::Text("ERR");
                }
                else {
                    ImGui::Text("%.1f", checker.updater.summary.dof_2);
                }
                if (checker.updater.model.covariance != Covariance::exchangeable && checker.updater.model.covariance != Covariance::nested_exchangeable){
                    ImGui::BulletText("Kenward-Roger (improved): "); ImGui::SameLine();
                    if (checker.updater.summary.power_kr2_2 == 909) {
                        ImGui::Text("ERR");
                    }
                    else {
                        ImGui::Text("%.1f", checker.updater.summary.dof_kr_2);
                    }
                }

                ImGui::Text("Treatment Interaction");
                ImGui::BulletText("Mixed model: "); ImGui::SameLine();
                ImGui::Text("%.1f", checker.updater.summary.dof_12);
                ImGui::BulletText("Mixed model t-test (between-within): "); ImGui::SameLine();
                ImGui::Text("%.1f", checker.updater.summary.dof_12);
                ImGui::BulletText("Satterthwaite: "); ImGui::SameLine();
                if (checker.updater.summary.power_sat_12 == 909) {
                    ImGui::Text("ERR");
                }
                else {
                    ImGui::Text("%.1f", checker.updater.summary.dof_kr_12);
                }
                ImGui::BulletText("Kenward-Roger: "); ImGui::SameLine();
                if (checker.updater.summary.power_kr_12 == 909) {
                    ImGui::Text("ERR");
                }
                else {
                    ImGui::Text("%.1f", checker.updater.summary.dof_kr_12);
                }
                ImGui::BulletText("GEE Independence working correlation, robust: "); ImGui::SameLine();
                if (checker.updater.summary.power_gee_indep_12 == 909) {
                    ImGui::Text("ERR");
                }
                else {
                    ImGui::Text("%.1f", checker.updater.summary.dof_12);
                }
                if (checker.updater.model.covariance != Covariance::exchangeable && checker.updater.model.covariance != Covariance::nested_exchangeable){
                    ImGui::BulletText("Kenward-Roger (improved): "); ImGui::SameLine();
                    if (checker.updater.summary.power_kr2_12 == 909) {
                        ImGui::Text("ERR");
                    }
                    else {
                        ImGui::Text("%.1f", checker.updater.summary.dof_kr_12);
                    }
                }
            }
        }

        enum Dmethod {
            Dmethod_glm,
            Dmethod_chisq
        };
        static int de_mode = 0;
        static int de_mode_curr_state = 0;

        if (!option.two_treatments && (checker.updater.model.covariance == ClusterApp::Covariance::exchangeable || checker.updater.model.covariance == ClusterApp::Covariance::nested_exchangeable)) {
            if (ImGui::CollapsingHeader("Design Effects")) {
                ImGui::Text("Design Effect Analysis"); ImGui::SameLine(); HelpMarker(
                    "For comparison we report the power and confidence interval half-width for exchangeable and nested exchangeable models using an adapted approach of Hooper et al (2016). Use of ICC, CAC, and IAC values for non - Gaussian - identity models uses the mean individual - level variance, which is approximated using the GLM weights, to convert to covariance parameter values.");


                ImGui::Text("Design effect: "); ImGui::SameLine();
                ImGui::Text("%.3f", checker.updater.summary.design_effect);
                ImGui::Text("Total observations: "); ImGui::SameLine();
                ImGui::Text("%.0f", checker.updater.summary.individual_n); ImGui::SameLine(); HelpMarker(
                    "The sample size for the single period, individual-level RCT. This value is corrected for the allocation ratio between treatment and control conditions.");
                if (checker.updater.model.family == ClusterApp::Family::binomial) {
                    if (ImGui::RadioButton("GLM variance", de_mode == Dmethod_glm)) { de_mode = Dmethod_glm; } ImGui::SameLine();
                    if (ImGui::RadioButton("Chi-squared test", de_mode == Dmethod_chisq)) { de_mode = Dmethod_chisq; }
                }
                else {
                    de_mode = Dmethod_glm;
                }
                checker.updater.de_mode = de_mode;

                if (de_mode_curr_state != de_mode) {
                    checker.updater.glmm.power_de(checker.updater.summary, de_mode);
                    de_mode_curr_state = de_mode;
                }
                if (checker.updater.summary.power_de == 909) {
                    ImGui::TextWrapped("Error. Check parameter values. For example, if the model produces probabilities outside [0,1] or if there are cells with zero probability for 2x2 test then this error will appear.");
                    checker.updater.log.AddLog("[%05d] [%s] Design effect power error. Int group mean: %.3f, Control group mean: %.3f, ratio: %.3f \n", ImGui::GetFrameCount(), checker.updater.log.cat[2], checker.updater.summary.se_de, checker.updater.summary.ci_width_de, checker.updater.summary.individual_se);
                    /*if (option.debug_info) {
                        ImGui::Text("Control group mean: "); ImGui::SameLine();
                        ImGui::Text("%.3f", checker.updater.summary.se_de);
                        ImGui::Text("Intervention group mean: "); ImGui::SameLine();
                        ImGui::Text("%.3f", checker.updater.summary.ci_width_de);
                        ImGui::Text("Ratio: "); ImGui::SameLine();
                        ImGui::Text("%.3f", checker.updater.summary.individual_se);
                    }*/
                }
                else {
                    if (de_mode == Dmethod_glm) {
                        ImGui::Text("Power: "); ImGui::SameLine();
                        ImGui::Text("%.1f", checker.updater.summary.power_de);
                        ImGui::Text("95%% Confidence-interval half-width: "); ImGui::SameLine();
                        ImGui::Text("%.3f", checker.updater.summary.ci_width_de);
                        ImGui::Text("Standard error: "); ImGui::SameLine();
                        ImGui::Text("%.3f", checker.updater.summary.se_de);
                        ImGui::Text("Assumed observation level variance: "); ImGui::SameLine();
                        ImGui::Text("%.3f", checker.updater.summary.individual_var);
                    }
                    else if (de_mode == Dmethod_chisq) {
                        ImGui::Text("Power: "); ImGui::SameLine();
                        ImGui::Text("%.1f", checker.updater.summary.power_de);
                    }
                }

            }

        }



        

        
        }

        




        ImGui::End();
    }

}

