#include "clusterapp.h"

namespace ClusterApp {

    void RenderPlotter(ClusterApp::plotData& plot, ClusterApp::options& option, ClusterApp::modelChecker& checker) {
        ImGui::Begin("Plotter");
        ImGui::Text("Plot settings");

        const char* xaxis_items_exchangeable[] = { "Clusters", "N per cluster-period", "ICC", "Treatment effect", "Baseline" };
        const char* xaxis_items_other[] = { "Clusters", "N per cluster-period", "ICC", "Treatment effect", "Baseline", "CAC" };
        static int xaxis_item_current = 2;
        static int series_item_current = 4;
        const char* yaxis_items[] = { "Power (GLS)", "CI width (GLS)", "Power (GLS-BW)", "CI width (GLS-BW)", "Power (Sat)", "CI width (Sat)", "Power (KR)", "CI width (KR)", "Power (Design effect)", "CI width (Design effect)" };
        static int yaxis_item_current = 0;
        static int series_clusters[] = { 1,2,3 };
        static int series_n[] = { 10,20,30 };
        static float series_icc[] = { 0.01,0.02,0.05 };
        static float series_te[] = { 0.0, 0.5, 1.0 };
        static float series_baseline[] = { 0.0, 0.5, 1.0 };
        static float series_cac[] = { 0.2,0.5,0.8 };
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
            plot.yaxis = ClusterApp::YAxis::power_sat;
            ylabel = "Power (%%)";
            break;
        case 5:
            plot.yaxis = ClusterApp::YAxis::ci_width_sat;
            ylabel = "95%% CI half-width";
            break;
        case 6:
            plot.yaxis = ClusterApp::YAxis::power_kr;
            ylabel = "Power (%%)";
            break;
        case 7:
            plot.yaxis = ClusterApp::YAxis::ci_width_kr;
            ylabel = "95%% CI half-width";
            break;
        case 8:
            plot.yaxis = ClusterApp::YAxis::power_de;
            ylabel = "Power (%%)";
            break;
        case 9:
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
            if (ImGui::Button(char_array, ImVec2(60, 40)))
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

        if (!option.auto_update && (checker.updater.plot_requires_update || checker.updater.requires_update)){
            ImGui::PushStyleColor(ImGuiCol_Button, colours.base1());
            ImGui::PushStyleColor(ImGuiCol_ButtonHovered, colours.base1(1));
            ImGui::PushStyleColor(ImGuiCol_ButtonActive, colours.base1(2));
            ImGui::PushStyleColor(ImGuiCol_Text, colours.base02(2));

            if (ImGui::Button("Refresh", ImVec2(80, 30))){
                checker.update();
            }
            ImGui::PopStyleColor(4);
        } else {
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
        }
        // if (option.debug_info) {
        //     ImGui::Text("X-axis limits:");
        //     ImGui::Text("%.3f", plot.x_axis_limits.first);
        //     ImGui::Text("%.3f", plot.x_axis_limits.second);
        // }


        ImGui::End();
    }
}

