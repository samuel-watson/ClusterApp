#include <stdio.h>

#ifdef __EMSCRIPTEN__
#include <emscripten.h>
#endif

#define GLFW_INCLUDE_ES3
#include <GLES3/gl3.h>
#include <GLFW/glfw3.h>

#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"
#include <iostream>
#include "src/clusterapp.h"


GLFWwindow* g_window;
ImVec4 clear_color = ImVec4(0.45f, 0.55f, 0.60f, 1.00f);
bool show_demo_window = true;
bool show_another_window = false;
int g_width;
int g_height;

// Function used by c++ to get the size of the html canvas
EM_JS(int, canvas_get_width, (), {
  return Module.canvas.width;
});

// Function used by c++ to get the size of the html canvas
EM_JS(int, canvas_get_height, (), {
  return Module.canvas.height;
});

// Function called by javascript
EM_JS(void, resizeCanvas, (), {
  js_resizeCanvas();
});

void on_size_changed()
{
  glfwSetWindowSize(g_window, g_width, g_height);

  ImGui::SetCurrentContext(ImGui::GetCurrentContext());
}

void loop()
{
  int width = canvas_get_width();
  int height = canvas_get_height();

  if (width != g_width || height != g_height)
  {
    g_width = width;
    g_height = height;
    on_size_changed();
  }

  glfwPollEvents();

  ImGui_ImplOpenGL3_NewFrame();
  ImGui_ImplGlfw_NewFrame();
  ImGui::NewFrame();

  ImGui::DockSpaceOverViewport(ImGui::GetMainViewport());

  static ClusterApp::design designs;
  static ClusterApp::options windows;
  static ClusterApp::statisticalModel model;
  static ClusterApp::colourPicker colour;
  static ClusterApp::modelSummary results(designs);
  static ClusterApp::glmmModel glmm(model, windows);
  static ClusterApp::modelUpdater updater(designs, model, results, glmm);
  static ClusterApp::modelChecker checker(designs, model, updater);

  ImGui::PushFont(main_font);
  if (windows.light_mode) {
      ImGui::PushStyleColor(ImGuiCol_WindowBg, colour.base3());
      ImGui::PushStyleColor(ImGuiCol_Text, colour.base03());
      ImGui::PushStyleColor(ImGuiCol_MenuBarBg, colour.base2());
      ImGui::PushStyleColor(ImGuiCol_TableHeaderBg, colour.base2());
      ImGui::PushStyleColor(ImGuiCol_PopupBg, colour.base2());
  }
  else {
      ImGui::PushStyleColor(ImGuiCol_WindowBg, colour.base03());
      ImGui::PushStyleColor(ImGuiCol_Text, colour.base1());
      ImGui::PushStyleColor(ImGuiCol_MenuBarBg, colour.base02());
      ImGui::PushStyleColor(ImGuiCol_TableHeaderBg, colour.base02());
      ImGui::PushStyleColor(ImGuiCol_PopupBg, colour.base02());
  }


  ClusterApp::ShowMainMenu(windows);
  ClusterApp::RenderDesigner(designs, updater, windows);
  if (windows.sample_size)
      ClusterApp::RenderSampleSize(designs);
  if (windows.model)
      ClusterApp::RenderModel(designs, model, windows, uni_font);
  if (windows.results)
      ClusterApp::RenderResults(updater, windows);
  if (windows.optimiser)
      ClusterApp::RenderOptimiser(designs, updater, results, windows);
  ImGui::PopFont();
  ImGui::PopStyleColor(5);
  // Rendering

  ImGui::Render();

  int display_w, display_h;
  glfwMakeContextCurrent(g_window);
  glfwGetFramebufferSize(g_window, &display_w, &display_h);
  glViewport(0, 0, display_w, display_h);
  glClearColor(clear_color.x, clear_color.y, clear_color.z, clear_color.w);
  glClear(GL_COLOR_BUFFER_BIT);

  ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
  glfwMakeContextCurrent(g_window);
}


int init_gl()
{
  if( !glfwInit() )
  {
      fprintf( stderr, "Failed to initialize GLFW\n" );
      return 1;
  }

  glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE); // We don't want the old OpenGL

  // Open a window and create its OpenGL context
  int canvasWidth = g_width;
  int canvasHeight = g_height;
  g_window = glfwCreateWindow(canvasWidth, canvasHeight, "Cluster Trials App", NULL, NULL);
  if( g_window == NULL )
  {
      fprintf( stderr, "Failed to open GLFW window.\n" );
      glfwTerminate();
      return -1;
  }
  glfwMakeContextCurrent(g_window); // Initialize GLEW

  return 0;
}


int init_imgui()
{
  // Setup Dear ImGui binding
  IMGUI_CHECKVERSION();
  ImGui::CreateContext();
  ImGui_ImplGlfw_InitForOpenGL(g_window, true);
  ImGui_ImplOpenGL3_Init();

  // Setup style
  ImGui::StyleColorsDark();

  ImGuiIO& io = ImGui::GetIO();

  // Load Fonts
  ImVector<ImWchar> ranges;
  ImFontGlyphRangesBuilder builder;                         // Add a specific character
  builder.AddRanges(io.Fonts->GetGlyphRangesGreek()); // Add one of the default ranges
  builder.AddChar(0x00B2);
  builder.BuildRanges(&ranges);

  ImFont* main_font = io.Fonts->AddFontFromFileTTF("data\twcen.ttf", 18.0f, nullptr, ranges.Data);
  ImFont* uni_font = io.Fonts->AddFontFromFileTTF("data\didact.ttf", 18.0f, nullptr, ranges.Data);
  io.Fonts->Build();

  resizeCanvas();

  return 0;
}


int init()
{
  init_gl();
  init_imgui();
  return 0;
}


void quit()
{
  glfwTerminate();
}


extern "C" int main(int argc, char** argv)
{
  g_width = canvas_get_width();
  g_height = canvas_get_height();
  if (init() != 0) return 1;

  #ifdef __EMSCRIPTEN__
  emscripten_set_main_loop(loop, 0, 1);
  #endif

  quit();

  return 0;
}
