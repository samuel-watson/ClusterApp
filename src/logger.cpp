#include "clusterapp.h"


namespace ClusterApp {
    void RenderLog(ClusterApp::AppLog& log) {
        ImGui::SetNextWindowSize(ImVec2(500, 400), ImGuiCond_FirstUseEver);
        ImGui::Begin("Log");
        ImGui::Text("Log");
        ImGui::End();

        log.Draw("Log");
    }
}
