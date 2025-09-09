#pragma once

#define IMGUI_DEFINE_MATH_OPERATORS

#include "clusterclasses.h"
#include <iomanip>
#include <sstream>
#include <format>
#include "implot.h"

namespace ClusterApp {

    bool MultiColorButton(const char* desc_id, const ImVec4& col_upper_left, const ImVec4& col_upper_right, const ImVec4& col_btm_left, const ImVec4& col_btm_right,
        ImGuiColorEditFlags flags, ImGuiButtonFlags button_flags, const ImVec2& size_arg, const ImVec4& col_upper_left_hover, const ImVec4& col_upper_right_hover,
        const ImVec4& col_btm_left_hover, const ImVec4& col_btm_right_hover, const ImVec4& col_upper_left_active, const ImVec4& col_upper_right_active,
        const ImVec4& col_btm_left_active, const ImVec4& col_btm_right_active);
    
    void AppDockSpace(bool* p_open);

    void RenderSelector(std::array<bool,3>& isin, ClusterApp::options& option);
    
    void RenderMenuBar(ClusterApp::appModel& model, ClusterApp::options& option, bool* mopen);
    
    void RenderDesignSelector(ClusterApp::appModel& model, ClusterApp::options& option);

    void RenderModelSelector(ClusterApp::appModel& model, ClusterApp::options& option);

    void RenderLog(ClusterApp::AppLog& log);

}

