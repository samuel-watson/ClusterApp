#include "clusterapp.h"

namespace ClusterApp {

	void RenderDataSim(ClusterApp::glmmModel& glmm, ClusterApp::options& option) {
		ImGui::Begin("Data simulator");

		ImGui::TextWrapped("Simulate data from the statistical model. This function may be useful for checking whether the assumed parameters are realistic.");

		static std::vector<double> data;

		if (ImGui::Button("Simulate")) {
			data.clear();
			data = glmm.sim_data();
		}

		if (data.size() > 0) {
			if (ImPlot::BeginPlot("Data histogram")) {
				ImPlot::SetupAxes("Cluster means", nullptr, ImPlotAxisFlags_AutoFit, ImPlotAxisFlags_AutoFit);
				ImPlot::SetNextFillStyle(IMPLOT_AUTO_COL, 0.5f);
				ImPlot::PlotHistogram("Empirical", data.data(), data.size());
				ImPlot::EndPlot();
			}
		}

		ImGui::End();
	}

}