#include "clusterapp.h"

namespace ClusterApp {

    void RenderDesigner(ClusterApp::design& designs, ClusterApp::modelUpdater& updater, ClusterApp::options& option) {
        ImGui::Begin("Trial Designer");

        // Variables and other setup for the designer
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
        enum Mode
        {
            Mode_Copy,
            Mode_Move,
            Mode_Swap
        };
        static int mode = 2;

        // View options for the designer
        if (ImGui::TreeNode("Drag and drop mode")) {            
            ImGui::Text("Drag and drop mode: "); ImGui::SameLine();
            if (ImGui::RadioButton("Copy", mode == Mode_Copy)) { mode = Mode_Copy; } ImGui::SameLine();
            if (ImGui::RadioButton("Move", mode == Mode_Move)) { mode = Mode_Move; } ImGui::SameLine();
            if (ImGui::RadioButton("Swap", mode == Mode_Swap)) { mode = Mode_Swap; }
            ImGui::TreePop();
        }

        ImGui::Spacing();

        // The designer 
        ImGui::Dummy(ImVec2(small_dim, small_dim)); ImGui::SameLine();
        ImGui::Dummy(ImVec2(small_dim, small_dim)); ImGui::SameLine();
        if (option.show_J_seq) {
            ImGui::Dummy(ImVec2(small_dim * 1.5, small_dim)); ImGui::SameLine();
        }
        ImGui::PushStyleColor(ImGuiCol_Button, colours.base1());
        ImGui::PushStyleColor(ImGuiCol_ButtonHovered, colours.base1(1));
        ImGui::PushStyleColor(ImGuiCol_ButtonActive, colours.base1(2));
        ImGui::PushStyleColor(ImGuiCol_Text, colours.base02(2));

        // First iterate over the time periods to add the + buttons to extend the columns
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
            ImGui::Dummy(ImVec2(small_dim * 1.5, small_dim)); ImGui::SameLine();
        }

        // now iterate again to create the column heads
        for (int t = 0; t < designs.time; t++) {
            ImGui::PushID(designs.time * designs.sequences + designs.time + 2 + designs.sequences + t);
            ImGui::Button(int_to_char(t + 1), ImVec2(large_dim, small_dim));
            if (ImGui::BeginPopupContextItem(NULL, ImGuiPopupFlags_MouseButtonLeft))
            {
                ImGui::Text("OPTIONS");
                if (ImGui::SmallButton("Activate all")) {
                    for (int s = 0; s < designs.sequences; s++)*designs.active(s, t) = true;
                }
                if (ImGui::SmallButton("De-activate all")) {
                    for (int s = 0; s < designs.sequences; s++)*designs.active(s, t) = false;
                }
                if (ImGui::SmallButton("Set all intervention")) {
                    for (int s = 0; s < designs.sequences; s++) {
                        *designs.intervention(s, t) = true;
                        *designs.intervention_2(s, t) = false;
                    }
                }
                if (option.dose_effect) {
                    static float dose_seq = 1.0;
                    ImGui::SetNextItemWidth(100);
                    ImGui::DragFloat("Dose", &dose_seq, 0.1f, 0.0, +FLT_MAX, "%.1f", ImGuiSliderFlags_None); ImGui::SameLine();
                    if (ImGui::SmallButton("Set dose")) {
                        for (int s = 0; s < designs.sequences; s++) {
                            *designs.dose(s, t) = dose_seq;
                        }
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

        // iterate over the sequences to produce the row heads and main buttons

        for (int n = 0; n < designs.sequences; n++)
        {
            if (option.show_J_seq) {
                std::string label_j = std::to_string(*designs.n_clusters(n));
                char* char_array_j = new char[label_j.length() + 1];
                strcpy(char_array_j, label_j.c_str());
                ImGui::PushID(20000 + n);
                ImGui::Button(char_array_j, ImVec2(small_dim * 1.5, small_dim));
                if (ImGui::BeginPopupContextItem(NULL, ImGuiPopupFlags_MouseButtonLeft)) {
                    static int n_clusters = 10;
                    ImGui::Text("Number of clusters");
                    ImGui::SetNextItemWidth(100);
                    ImGui::InputScalar("N", ImGuiDataType_S16, &n_clusters, &s16_one, NULL, "%d");
                    if (ImGui::Button("Set")) {
                        *(designs.n_clusters(n)) = n_clusters;
                    }
                    ImGui::EndPopup();
                }
                ImGui::PopID();
                ImGui::SameLine();
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
            // Row heads
            ImGui::PushID(designs.time * (designs.sequences + 1) + designs.time + designs.sequences + 2 + n);
            ImGui::Button(int_to_char(n + 1), ImVec2(small_dim, large_dim)); ImGui::SameLine();
            if (ImGui::BeginPopupContextItem(NULL, ImGuiPopupFlags_MouseButtonLeft))
            {
                static int seq_n = 10;
                static float dose = 1.0;
                ImGui::Text("Number of clusters");
                ImGui::SetNextItemWidth(100);
                ImGui::InputScalar("N", ImGuiDataType_S16, &seq_n, &s16_one, NULL, "%d"); ImGui::SameLine();
                if (ImGui::Button("Set")) {
                    *(designs.n_clusters(n)) = seq_n;
                }
                ImGui::Text("OPTIONS");
                if (ImGui::SmallButton("Activate all")) {
                    for (int s = 0; s < designs.time; s++)*designs.active(n, s) = true;
                }
                if (ImGui::SmallButton("De-activate all")) {
                    for (int s = 0; s < designs.time; s++)*designs.active(n, s) = false;
                }
                if (ImGui::SmallButton("Set all intervention")) {
                    for (int s = 0; s < designs.time; s++) {
                        *designs.intervention(n, s) = true;
                        *designs.intervention_2(n, s) = false;
                    }
                }
                if (option.dose_effect) {
                    ImGui::SetNextItemWidth(100);
                    ImGui::DragFloat("Dose", &dose, 0.1f, 0.0, +FLT_MAX, "%.1f", ImGuiSliderFlags_None); ImGui::SameLine();
                    if (ImGui::SmallButton("Set dose")) {
                        for (int s = 0; s < designs.time; s++) {
                            *designs.dose(n, s) = dose;
                        }
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
            
            // cluster-period buttons
            for (int t = 0; t < designs.time; t++) {
                int id = t + n * designs.time;
                ImGui::PushID(id);
                std::string label = "";
                if (*designs.active(n, t)) {
                    if (option.show_n_period) {
                        label += "n=";
                        label += std::to_string(*designs.n(n, t));
                        if (option.show_status_period) label += "\n"; //replace with | for single line
                    }
                    if (option.show_status_period) {                        
                        label += option.dose_effect ? std::format("{:.2f}", (*designs.intervention(n, t)) * (*designs.dose(n, t))) : std::to_string(*designs.intervention(n, t));
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
                        static int cell_n = 10;
                        static float dose_n = 1.0;
                        ImGui::Checkbox("Active", designs.active(n, t));
                        ImGui::Checkbox("Intervention", designs.intervention(n, t));
                        if (option.two_treatments)ImGui::Checkbox("Intervention 2", designs.intervention_2(n, t));
                        if (option.dose_effect && *designs.intervention(n, t)) {
                            ImGui::SetNextItemWidth(100);
                            ImGui::DragFloat("Dose", &dose_n, 0.1f, 0.0, +FLT_MAX, "%.1f", ImGuiSliderFlags_None);
                        }
                        ImGui::SetNextItemWidth(100);
                        ImGui::InputScalar("N", ImGuiDataType_S16, &cell_n, &s16_one, NULL, "%d"); ImGui::SameLine();
                        if (ImGui::Button("Set")) {
                            *(designs.n(n, t)) = cell_n;
                            if (option.dose_effect && *designs.intervention(n, t)) *designs.dose(n, t) = dose_n;
                        }
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

        /*if (option.debug_info) {
            ImGui::Text("Design Checksum"); ImGui::SameLine();
            ImGui::Text("%d", designs.crc_val);
        }*/

        ImGui::End();
    }

    
}

