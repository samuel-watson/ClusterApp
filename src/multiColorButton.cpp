#include "clusterapp.h"

bool ClusterApp::MultiColorButton(const char* desc_id, const ImVec4& col_upper_left, const ImVec4& col_upper_right, const ImVec4& col_btm_left, const ImVec4& col_btm_right,
    ImGuiColorEditFlags flags, ImGuiButtonFlags button_flags, const ImVec2& size_arg, const ImVec4& col_upper_left_hover, const ImVec4& col_upper_right_hover,
    const ImVec4& col_btm_left_hover, const ImVec4& col_btm_right_hover, const ImVec4& col_upper_left_active, const ImVec4& col_upper_right_active,
    const ImVec4& col_btm_left_active, const ImVec4& col_btm_right_active)
{
    using namespace ImGui;
    ImGuiWindow* window = GetCurrentWindow();
    if (window->SkipItems)
        return false;

    ImGuiContext& g = *GImGui;
    const ImGuiStyle& style = g.Style;
    const ImVec2 label_size = CalcTextSize(desc_id, NULL, true);
    const ImGuiID id = window->GetID(desc_id);
    const float default_size = GetFrameHeight();
    const ImVec2 size(size_arg.x == 0.0f ? default_size : size_arg.x, size_arg.y == 0.0f ? default_size : size_arg.y);
    const ImRect bb(window->DC.CursorPos, window->DC.CursorPos + size);
    ItemSize(bb, (size.y >= default_size) ? g.Style.FramePadding.y : 0.0f);
    if (!ItemAdd(bb, id))
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
