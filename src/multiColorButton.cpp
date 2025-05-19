#include "clusterapp.h"

bool ClusterApp::MultiColorButton(const char* desc_id, const ImVec4& col_upper_left, const ImVec4& col_upper_right, const ImVec4& col_btm_left, const ImVec4& col_btm_right,
    ImGuiColorEditFlags flags, ImGuiButtonFlags button_flags, const ImVec2& size_arg, const ImVec4& col_upper_left_hover, const ImVec4& col_upper_right_hover,
    const ImVec4& col_btm_left_hover, const ImVec4& col_btm_right_hover, const ImVec4& col_upper_left_active, const ImVec4& col_upper_right_active,
    const ImVec4& col_btm_left_active, const ImVec4& col_btm_right_active)
{
    using namespace ImGui;
    ImGuiWindow* window = ImGui::GetCurrentWindow();
    if (window->SkipItems)
        return false;

    ImGuiContext& g = *GImGui;
    const ImGuiStyle& style = g.Style;
    const ImVec2 label_size = CalcTextSize(desc_id, NULL, true);
    const ImGuiID id = window->GetID(desc_id);
    const float default_size = GetFrameHeight();
    const ImVec2 size(size_arg.x == 0.0f ? default_size : size_arg.x, size_arg.y == 0.0f ? default_size : size_arg.y);
    const ImRect bb(window->DC.CursorPos, window->DC.CursorPos + size);
    ImGui::ItemSize(bb, (size.y >= default_size) ? g.Style.FramePadding.y : 0.0f);
    if (!ImGui::ItemAdd(bb, id))
        return false;

    bool hovered, held;
    bool pressed = ButtonBehavior(bb, id, &hovered, &held, button_flags);

    if (flags & ImGuiColorEditFlags_NoAlpha)
        flags &= ~(ImGuiColorEditFlags_AlphaPreview | ImGuiColorEditFlags_AlphaPreviewHalf);
    ImVec4 col_rgb_upper_left;
    ImVec4 col_rgb_upper_right;
    ImVec4 col_rgb_btm_left;
    ImVec4 col_rgb_btm_right;
    if (hovered) {
        col_rgb_upper_left = col_upper_left_hover;
        col_rgb_upper_right = col_upper_right_hover;
        col_rgb_btm_left = col_btm_left_hover;
        col_rgb_btm_right = col_btm_right_hover;
    }
    else if (held) {
        col_rgb_upper_left = col_upper_left_active;
        col_rgb_upper_right = col_upper_right_active;
        col_rgb_btm_left = col_btm_left_active;
        col_rgb_btm_right = col_btm_right_active;
    }
    else {
        col_rgb_upper_left = col_upper_left;
        col_rgb_upper_right = col_upper_right;
        col_rgb_btm_left = col_btm_left;
        col_rgb_btm_right = col_btm_right;
    }

    if (flags & ImGuiColorEditFlags_InputHSV) {
        ColorConvertHSVtoRGB(col_rgb_upper_left.x, col_rgb_upper_left.y, col_rgb_upper_left.z, col_rgb_upper_left.x, col_rgb_upper_left.y, col_rgb_upper_left.z);
        ColorConvertHSVtoRGB(col_rgb_upper_right.x, col_rgb_upper_right.y, col_rgb_upper_right.z, col_rgb_upper_right.x, col_rgb_upper_right.y, col_rgb_upper_right.z);
        ColorConvertHSVtoRGB(col_rgb_btm_left.x, col_rgb_btm_left.y, col_rgb_btm_left.z, col_rgb_btm_left.x, col_rgb_btm_left.y, col_rgb_btm_left.z);
        ColorConvertHSVtoRGB(col_rgb_btm_right.x, col_rgb_btm_right.y, col_rgb_btm_right.z, col_rgb_btm_right.x, col_rgb_btm_right.y, col_rgb_btm_right.z);

    }

    float grid_step = ImMin(size.x, size.y) / 2.99f;
    float rounding = ImMin(g.Style.FrameRounding, grid_step * 0.5f);
    ImRect bb_inner = bb;
    float off = 0.0f;
    if ((flags & ImGuiColorEditFlags_NoBorder) == 0)
    {
        off = -0.75f; // The border (using Col_FrameBg) tends to look off when color is near-opaque and rounding is enabled. This offset seemed like a good middle ground to reduce those artifacts.
        bb_inner.Expand(off);
    }

    window->DrawList->AddRectFilledMultiColor(bb_inner.Min, bb_inner.Max, GetColorU32(col_rgb_upper_left), GetColorU32(col_rgb_upper_right),
        GetColorU32(col_rgb_btm_left), GetColorU32(col_rgb_btm_right));

    RenderNavHighlight(bb, id);
    //RenderFrame(bb.Min, bb.Max, col, true, style.FrameRounding);

    if (g.LogEnabled)
        LogSetNextTextDecoration("[", "]");
    RenderTextClipped(bb.Min + style.FramePadding, bb.Max - style.FramePadding, desc_id, NULL, &label_size, style.ButtonTextAlign, &bb);

    IMGUI_TEST_ENGINE_ITEM_INFO(id, desc_id, g.LastItemData.StatusFlags);

    return pressed;
}
/*
bool ClusterApp::BufferingBar(const char* label, float value,  const ImVec2& size_arg, const ImU32& bg_col, const ImU32& fg_col) {
    ImGuiWindow* window = ImGui::GetCurrentWindow();
    if (window->SkipItems)
        return false;
    
    ImGuiContext& g = *GImGui;
    const ImGuiStyle& style = g.Style;
    const ImGuiID id = window->GetID(label);

    ImVec2 pos = window->DC.CursorPos;
    ImVec2 size = size_arg;
    size.x -= style.FramePadding.x * 2;
    
    const ImRect bb(pos, ImVec2(pos.x + size.x, pos.y + size.y));
    ImGui::ItemSize(bb, style.FramePadding.y);
    if (!ImGui::ItemAdd(bb, id))
        return false;
    
    // Render
    const float circleStart = size.x * 0.7f;
    const float circleEnd = size.x;
    const float circleWidth = circleEnd - circleStart;
    
    window->DrawList->AddRectFilled(bb.Min, ImVec2(pos.x + circleStart, bb.Max.y), bg_col);
    window->DrawList->AddRectFilled(bb.Min, ImVec2(pos.x + circleStart*value, bb.Max.y), fg_col);
    
    const float t = g.Time;
    const float r = size.y / 2;
    const float speed = 1.5f;
    
    const float a = speed*0;
    const float b = speed*0.333f;
    const float c = speed*0.666f;
    
    const float o1 = (circleWidth+r) * (t+a - speed * (int)((t+a) / speed)) / speed;
    const float o2 = (circleWidth+r) * (t+b - speed * (int)((t+b) / speed)) / speed;
    const float o3 = (circleWidth+r) * (t+c - speed * (int)((t+c) / speed)) / speed;
    
    window->DrawList->AddCircleFilled(ImVec2(pos.x + circleEnd - o1, bb.Min.y + r), r, bg_col);
    window->DrawList->AddCircleFilled(ImVec2(pos.x + circleEnd - o2, bb.Min.y + r), r, bg_col);
    window->DrawList->AddCircleFilled(ImVec2(pos.x + circleEnd - o3, bb.Min.y + r), r, bg_col);
}

bool ClusterApp::Spinner(const char* label, float radius, int thickness, const ImU32& color) {
    ImGuiWindow* window = ImGui::GetCurrentWindow();
    if (window->SkipItems)
        return false;
    
    ImGuiContext& g = *GImGui;
    const ImGuiStyle& style = g.Style;
    const ImGuiID id = window->GetID(label);
    
    ImVec2 pos = window->DC.CursorPos;
    ImVec2 size((radius )*2, (radius + style.FramePadding.y)*2);
    
    const ImRect bb(pos, ImVec2(pos.x + size.x, pos.y + size.y));
    ImGui::ItemSize(bb, style.FramePadding.y);
    if (!ImGui::ItemAdd(bb, id))
        return false;
    
    // Render
    window->DrawList->PathClear();
    
    int num_segments = 30;
    int start = abs(ImSin(g.Time*1.8f)*(num_segments-5));
    
    const float a_min = IM_PI*2.0f * ((float)start) / (float)num_segments;
    const float a_max = IM_PI*2.0f * ((float)num_segments-3) / (float)num_segments;

    const ImVec2 centre = ImVec2(pos.x+radius, pos.y+radius+style.FramePadding.y);
    
    for (int i = 0; i < num_segments; i++) {
        const float a = a_min + ((float)i / (float)num_segments) * (a_max - a_min);
        window->DrawList->PathLineTo(ImVec2(centre.x + ImCos(a+g.Time*8) * radius,
                                            centre.y + ImSin(a+g.Time*8) * radius));
    }

    window->DrawList->PathStroke(color, false, thickness);
}*/