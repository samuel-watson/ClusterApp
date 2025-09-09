#include <stdio.h>

#ifdef __EMSCRIPTEN__
#include <emscripten.h>
#endif

#define GLFW_INCLUDE_ES3
#include <GLES3/gl3.h>
#include <GLFW/glfw3.h>

#define IMGUI_DEFINE_MATH_OPERATORS

#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"
#include "imgui_spectrum.h"
#include <iostream>
#include "clusterapp.h"

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
  
  static ClusterApp::options option;
  static ClusterApp::AppLog logger;
  static std::array<bool, 3>  isin{true,false,false};
  static ClusterApp::appModel model1(option, logger, 0);
  static ClusterApp::appModel model2(option, logger, 1);
  static ClusterApp::appModel model3(option, logger, 2);

  static ClusterApp::colourPicker colour;
  ImGui::SetNextWindowSize(ImVec2(160, 500), ImGuiCond_FirstUseEver);
  ImGui::SetNextWindowPos(ImVec2(0, 0), ImGuiCond_FirstUseEver);
  ClusterApp::RenderSelector(isin, option);

  static float maxx = ImGui::GetWindowPos().x;
  static float xwidth = maxx > 700.0 ? maxx : 700.0;

  int counter = 0;
  for(int i = 0; i < 3; i++){
    if(isin[i]){
      
      ImGui::SetNextWindowPos(ImVec2((counter*50)+175, 50+(counter*50)), ImGuiCond_FirstUseEver);
      if(i == 0){
        model1.checker.check_time();
        if(!model1.first_design_click && !model1.first_model_click && isin[0]){
          ImGui::SetNextWindowSize(ImVec2(300,350), ImGuiCond_FirstUseEver);
          ClusterApp::RenderDesignSelector(model1, option);
        } else if(model1.first_design_click && !model1.first_model_click && isin[0]) {
          ImGui::SetNextWindowSize(ImVec2(500,600), ImGuiCond_FirstUseEver);
          ClusterApp::RenderModelSelector(model1, option);
        } else {
          ImGui::SetNextWindowSize(ImVec2((int)xwidth,800), ImGuiCond_FirstUseEver);
          ClusterApp::RenderMenuBar(model1, option, &isin[0]);
        }        
      } else if(i == 1){
        model2.checker.check_time();
        if(!model2.first_design_click && !model2.first_model_click && isin[1]){
          ImGui::SetNextWindowSize(ImVec2(300,350), ImGuiCond_FirstUseEver);
          ClusterApp::RenderDesignSelector(model2, option);
        } else if(model2.first_design_click && !model2.first_model_click && isin[1]) {
          ImGui::SetNextWindowSize(ImVec2(500,600), ImGuiCond_FirstUseEver);
          ClusterApp::RenderModelSelector(model2, option);
        } else {
          ImGui::SetNextWindowSize(ImVec2((int)xwidth,800), ImGuiCond_FirstUseEver);
          ClusterApp::RenderMenuBar(model2, option, &isin[1]);
        } 
      } else {
        model3.checker.check_time();
        if(!model3.first_design_click && !model3.first_model_click && isin[2]){
          ImGui::SetNextWindowSize(ImVec2(300,350), ImGuiCond_FirstUseEver);
          ClusterApp::RenderDesignSelector(model3, option);
        } else if(model3.first_design_click && !model3.first_model_click && isin[2]) {
          ImGui::SetNextWindowSize(ImVec2(500,600), ImGuiCond_FirstUseEver);
          ClusterApp::RenderModelSelector(model3, option);
        } else {
          ImGui::SetNextWindowSize(ImVec2((int)xwidth,800), ImGuiCond_FirstUseEver);
          ClusterApp::RenderMenuBar(model3, option, &isin[2]);
        } 
      }
      counter++;
    }    
  }

  if (option.log) {
      ClusterApp::RenderLog(logger);
  }
      

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
  g_window = glfwCreateWindow(canvasWidth, canvasHeight, "Cluster Trials Designer", NULL, NULL);
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
  ImPlot::CreateContext();
  ImGui_ImplGlfw_InitForOpenGL(g_window, true);
  ImGui_ImplOpenGL3_Init();

  // Setup style
  ImGui::StyleColorsDark();

  ImGuiIO& io = ImGui::GetIO();
  io.ConfigFlags |= ImGuiConfigFlags_DockingEnable;

  ImGui::Spectrum::StyleColorsSpectrum();
  ImGui::Spectrum::LoadFont(18);
  // Load Fonts
  /*
  ImVector<ImWchar> ranges;
  ImFontGlyphRangesBuilder builder;                      
  builder.AddRanges(io.Fonts->GetGlyphRangesGreek()); // Add one of the default ranges - doesn't actually have greek but leaving here in case I change the font
  builder.AddChar(0x00B2);
  builder.BuildRanges(&ranges);

  io.Fonts->AddFontFromFileTTF("data/twcen.ttf", 18.0f, nullptr, ranges.Data);
  //ImFont* uni_font = io.Fonts->AddFontFromFileTTF("data/didact.ttf", 18.0f, nullptr, ranges.Data); // this has greek but haven't figured out how to choose fonts
  io.Fonts->Build();*/

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
