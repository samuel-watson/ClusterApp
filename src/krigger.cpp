#include "clusterapp.h"

namespace ClusterApp {

    void RenderKriging(ClusterApp::krigingData& krig, ClusterApp::options& option) {
        ImGui::Begin("Kriging");
        ImGui::Text("Estimates the power across the sample size space");
        ImGui::TextWrapped("This feature is experimental. The power surface is estimated using a kriging approach and estimates will improve by sampling more points. The number of clusters is divided \
according to the proportions given by the main trial design. Individuals are the number of individuals per cluster-period.");
        ImGui::TextWrapped("The power plot may not meet exactly the specified plot limits due to rounding of the increments in the plot and currently only provides power for GLS estimators.");
        //const char* power_items[] = { "GLS", "Between-within", "Satterthwaite", "Kenward-Roger", "Design effect"};
        //static int power_item_current = 0;
        // plot of sampled points
        const short  s16_zero = 0, s16_one = 1;
        static int new_sample_size = 25;
        colourPicker colours;

        if(!krig.start)krig.update(false);
        //ImGui::SetNextItemWidth(200);
        //ImGui::Combo("Power: ", &power_item_current, power_items, IM_ARRAYSIZE(power_items));
        ImGui::SetNextItemWidth(150);
        ImGui::InputScalar("Sample size", ImGuiDataType_S16, &new_sample_size, &s16_one, NULL, "%d");
        ImGui::SetNextItemWidth(150);
        ImGui::InputScalar("Clusters min.", ImGuiDataType_S16, &krig.lower_int[0], &s16_one, NULL, "%d"); ImGui::SameLine();
        ImGui::SetNextItemWidth(150);
        ImGui::InputScalar("Clusters max.", ImGuiDataType_S16, &krig.upper_int[0], &s16_one, NULL, "%d");
        ImGui::SetNextItemWidth(150);
        ImGui::InputScalar("Individual min.", ImGuiDataType_S16, &krig.lower_int[1], &s16_one, NULL, "%d"); ImGui::SameLine();
        ImGui::SetNextItemWidth(150);
        ImGui::InputScalar("Individual max.", ImGuiDataType_S16, &krig.upper_int[1], &s16_one, NULL, "%d");
        if (ImGui::Button("Generate new sample")) {
            //krig.set_power_type(static_cast<PowerType>(power_item_current));
            krig.generate_grid();
            krig.new_sample(new_sample_size);
        }

        ImGui::Text("Sample size = %lu", krig.n_ind.size());
        if(option.debug_info)ImGui::Text("Mu = %.3f", krig.mu);
;       ImGui::Text("Sampled points");

        static ImPlotAxisFlags xflags = ImPlotAxisFlags_AutoFit; 
        static ImPlotAxisFlags yflags = ImPlotAxisFlags_AutoFit;

        if (ImPlot::BeginPlot("Scatter Plot", ImVec2(250, 250), ImPlotFlags_NoLegend | ImPlotFlags_NoMouseText)) {
            ImPlot::SetupAxis(ImAxis_X1, "Clusters", xflags);
            ImPlot::SetupAxis(ImAxis_Y1, "Individuals", yflags);
            ImPlot::PushStyleColor(ImPlotCol_Line, colours.red());
            ImPlot::PlotScatter("Sampled points", krig.n_cl.data(), krig.n_ind.data(), krig.n_ind.size());
            ImPlot::EndPlot();
        }
        if (option.debug_info) {
            ImGui::Text("Last sampled point: (%i", krig.n_ind[krig.n_ind.size() - 1]); ImGui::SameLine();
            ImGui::Text(" , %i )", krig.n_cl[krig.n_cl.size() - 1]);
        }
        ImGui::SetNextItemWidth(200);
        ImGui::DragFloat("Bandwidth", &krig.bandwidth, 0.01f, 0.0f, 2.0f, "%.2f", ImGuiSliderFlags_None); ImGui::SameLine(); HelpMarker(
            "The bandwidth controls the level of smoothing for the power surface. Sample size values are scaled to [0,1] x [0,1].");
        
        ImGui::SetNextItemWidth(200);
        ImGui::DragFloat("Power threshold", &krig.threshold_power, 0.01f, 0.0f, 1.0f, "%0.2f", ImGuiSliderFlags_None); ImGui::SameLine(); HelpMarker(
            "The points are sampled to minimised the uncertainty about where the contour for this level of power is.");
 
        if (ImGui::Button("Update heat map")) {
            krig.update(false);
        }
        ImGui::SameLine();
        if (ImGui::Button("Sample new point")) {
            krig.update(true);
        }
        
        //static const char* xlabels[] = { "10","15","20","25","30","40"};
        //static const char* ylabels[] = { "R1","R2","R3","R4","R5","R6","R7" };

        if (krig.surface_initialised) {
            static ImPlotColormap map = ImPlotColormap_Viridis;
            if (ImPlot::ColormapButton(ImPlot::GetColormapName(map), ImVec2(225, 0), map)) {
                map = (map + 1) % ImPlot::GetColormapCount();
                ImPlot::BustColorCache("Power plot");
            }

            ImPlot::PushColormap(map);

            static ImPlotAxisFlags axes_flags = ImPlotAxisFlags_Lock | ImPlotAxisFlags_NoGridLines;
            if (ImPlot::BeginPlot("Power plot", ImVec2(800,500), ImPlotFlags_NoLegend | ImPlotFlags_NoMouseText)) {
                ImPlot::SetupAxes("Clusters", "Individuals", axes_flags, axes_flags);
                ImPlot::SetupAxisTicks(ImAxis_X1, 0 + 1.0 / 40.0, 1 - 1.0 / 40.0, 20, krig.n_cl_grid_label);
                ImPlot::SetupAxisTicks(ImAxis_Y1, 0 + 1.0 / 40.0, 1 - 1.0 / 40.0, 20, krig.n_ind_grid_label);
                ImPlot::PlotHeatmap("Power", krig.surface, 20,20,0,0,"%.0f");
                ImPlot::EndPlot();
            }
        }

        ImGui::End();
    }
}