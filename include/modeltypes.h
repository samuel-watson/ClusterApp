#pragma once

#include<string>
#include<cstring>
#include<boost/crc.hpp>
#include "imgui.h"
#include <boost/math/distributions/students_t.hpp>
#include <boost/math/distributions/non_central_chi_squared.hpp>
#include <boost/math/distributions/chi_squared.hpp>

namespace ClusterApp {

    inline char* int_to_char(const int& i) {
        std::string label = std::to_string(i);
        char* char_array = new char[label.length() + 1];
        strcpy(char_array, label.c_str());
        return char_array;
    }

    static void HelpMarker(const char* desc)
    {
        ImGui::TextDisabled("(?)");
        if (ImGui::BeginItemTooltip())
        {
            ImGui::PushTextWrapPos(ImGui::GetFontSize() * 35.0f);
            ImGui::TextUnformatted(desc);
            ImGui::PopTextWrapPos();
            ImGui::EndTooltip();
        }
    }

    struct options {
        bool sample_size = false;
        bool two_treatments = false;
        bool dose_effect = false;
        bool heterogeneous_te = false;
        bool time_averaged_te = false;
        bool model = true;
        bool results = true;
        bool optimiser = false;
        bool light_mode = true;
        bool dockspace = true;
        bool plotter = false;
        bool krigger = false;
        bool simulate = false;
        bool log = false;
        bool show_n_period = true;
        bool show_status_period = false;
        bool show_J_seq = true;
        bool use_icc_for_non_gaussian = true;
        bool debug_info = false;
        bool show_box = false;
        bool auto_update = false;
        options() {};
        ~options() = default;
    };

    enum class Family {
        gaussian = 1,
        binomial = 2,
        poisson = 3,
        beta = 4,
        gamma = 5,
        quantile = 6
    };

    enum class Link {
        identity = 1,
        log = 2,
        logit = 3,
        inverse = 4,
        probit = 5
    };

    enum class Covariance {
        exchangeable = 1,
        nested_exchangeable = 2,
        autoregressive = 3,
        exponential = 4,
        squared_exponential = 5
    };

    enum class IndividualCovariance {
        exchangeable = 1,
        autoregressive = 2,
        exponential = 3,
        squared_exponential = 4
    };

    enum class LinearPredictor {
        time_fixed_effects = 1,
        cluster_linear_trends = 2
    };

    enum class Sampling {
        cross_sectional = 1,
        cohort = 2,
        open_cohort = 3
    };

    // Need to add the other options below to the plotdata class
    enum class XAxis {
        clusters = 1,
        individual_n = 2,
        /*cov_par_1 = 3,
        cov_par_2 = 4,
        cov_par_3 = 5,
        cov_par_4 = 6,
        cov_par_5 = 7,*/
        treatment_effect = 8,
        baseline = 9,
        icc = 10,
        cac = 11
        //iac = 12,
        //total_n = 13
    };

    enum class YAxis {
        power = 0,
        ci_width = 1,
        power_bw = 2,
        ci_width_bw = 3,
        power_sat = 4,
        ci_width_sat = 5,
        power_kr = 6,
        ci_width_kr = 7,
        power_de = 8,
        ci_width_de = 9
    };

    struct CRC {
        boost::crc_32_type crc;
        template <typename T>
        void operator()(T const& i) {
            crc.process_bytes(&i, sizeof(i));
        }
        auto get() { return crc.checksum(); }
    };

    struct colourPicker {
        ImVec4 red(const int type = 0) {
            if (type == 0) {
                return (ImVec4)ImColor::HSV(0.002f, 0.79f, 0.86f);
            }
            else if (type == 1) {
                return (ImVec4)ImColor::HSV(0.002f, 0.79f, 0.96f);
            }
            else {
                return (ImVec4)ImColor::HSV(0.002f, 0.79f, 1.06f);
            }
        }
        ImVec4 blue(const int type = 0) {
            if (type == 0) {
                return (ImVec4)ImColor::HSV(0.569f, 0.82f, 0.82f);
            }
            else if (type == 1) {
                return (ImVec4)ImColor::HSV(0.569f, 0.82f, 0.92f);
            }
            else {
                return (ImVec4)ImColor::HSV(0.569f, 0.82f, 1.02f);
            }

        }
        ImVec4 base01(const int type = 0) {
            if (type == 0) {
                return (ImVec4)ImColor::HSV(0.54f, 0.25f, 0.46f);
            }
            else if (type == 1) {
                return (ImVec4)ImColor::HSV(0.54f, 0.25f, 0.56f);
            }
            else {
                return (ImVec4)ImColor::HSV(0.54f, 0.25f, 0.66f);
            }
        }
        ImVec4 base02(const int type = 0) {
            if (type == 0) {
                return (ImVec4)ImColor::HSV(0.534f, 0.894f, 0.259f);
            }
            else if (type == 1) {
                return (ImVec4)ImColor::HSV(0.534f, 0.894f, 0.359f);
            }
            else {
                return (ImVec4)ImColor::HSV(0.534f, 0.894f, 0.459f);
            }
        }
        ImVec4 base03(const int type = 0) {
            if (type == 0) {
                return (ImVec4)ImColor::HSV(0.534f, 1.000f, 0.212f);
            }
            else if (type == 1) {
                return (ImVec4)ImColor::HSV(0.534f, 1.000f, 0.312f);
            }
            else {
                return (ImVec4)ImColor::HSV(0.534f, 1.000f, 0.412f);
            }
        }
        ImVec4 base2(const int type = 0) {
            if (type == 0) {
                return (ImVec4)ImColor::HSV(0.127f, 0.105f, 0.933f);
            }
            else if (type == 1) {
                return (ImVec4)ImColor::HSV(0.127f, 0.105f, 1.033f);
            }
            else {
                return (ImVec4)ImColor::HSV(0.127f, 0.105f, 1.133f);
            }
        }
        ImVec4 base3(const int type = 0) {
            if (type == 0) {
                return (ImVec4)ImColor::HSV(0.123f, 0.103f, 0.992f);
            }
            else if (type == 1) {
                return (ImVec4)ImColor::HSV(0.123f, 0.103f, 1.002f);
            }
            else {
                return (ImVec4)ImColor::HSV(0.123f, 0.103f, 1.102f);
            }
        }
        ImVec4 base1(const int type = 0) {
            if (type == 0) {
                return (ImVec4)ImColor::HSV(0.5f, 0.087f, 0.631f);
            }
            else if (type == 1) {
                return (ImVec4)ImColor::HSV(0.5f, 0.087f, 0.731f);
            }
            else {
                return (ImVec4)ImColor::HSV(0.5f, 0.087f, 0.831f);
            }
        }
        ImVec4 yellow(const int type = 0) {
            if (type == 0) {
                return (ImVec4)ImColor::HSV(0.126f, 1.000f, 0.710f);
            }
            else if (type == 1) {
                return (ImVec4)ImColor::HSV(0.126f, 1.000f, 0.810f);
            }
            else {
                return (ImVec4)ImColor::HSV(0.126f, 1.000f, 0.910f);
            }
        }
        ImVec4 green(const int type = 0) {
            if (type == 0) {
                return (ImVec4)ImColor::HSV(0.189f, 1.000f, 0.600f);
            }
            else if (type == 1) {
                return (ImVec4)ImColor::HSV(0.189f, 1.000f, 0.700f);
            }
            else {
                return (ImVec4)ImColor::HSV(0.189f, 1.000f, 0.800f);
            }
        }
    };

}
