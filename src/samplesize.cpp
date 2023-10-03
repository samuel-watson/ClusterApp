#include "clusterapp.h"

namespace ClusterApp {

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

 
}

