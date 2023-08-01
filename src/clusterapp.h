﻿#pragma once

#define IMGUI_DEFINE_MATH_OPERATORS

#include "imgui_internal.h"
#include "clusterclasses.h"

namespace ClusterApp {

    bool MultiColorButton(const char* desc_id, const ImVec4& col_upper_left, const ImVec4& col_upper_right, const ImVec4& col_btm_left, const ImVec4& col_btm_right,
        ImGuiColorEditFlags flags, ImGuiButtonFlags button_flags, const ImVec2& size_arg, const ImVec4& col_upper_left_hover, const ImVec4& col_upper_right_hover,
        const ImVec4& col_btm_left_hover, const ImVec4& col_btm_right_hover, const ImVec4& col_upper_left_active, const ImVec4& col_upper_right_active,
        const ImVec4& col_btm_left_active, const ImVec4& col_btm_right_active);

    void ShowMainMenu(ClusterApp::options& windows, ClusterApp::design& designs, ClusterApp::modelUpdater& updater, ClusterApp::modelSummary& summary);

    void AppDockSpace(bool* p_open);

    void RenderDesigner(ClusterApp::design& designs, ClusterApp::modelUpdater& updater, ClusterApp::options& option);

    void RenderSampleSize(ClusterApp::design& designs);

    void RenderModel(ClusterApp::design& design, ClusterApp::statisticalModel& model, ClusterApp::options& option);

    void RenderResults(ClusterApp::modelUpdater& updater, ClusterApp::options& option);

    void RenderOptimiser(ClusterApp::design& design, ClusterApp::modelUpdater& updater, ClusterApp::modelSummary& summary, ClusterApp::options& option);
}

