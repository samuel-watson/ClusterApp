CXX = emcc
OUTPUT = imgui.html
IMGUI_DIR:=C:/cpp/imgui/imgui
EIGEN_DIR:=C:/cpp/eigen/eigen-3.4.0
BOOST_DIR:=C:/cpp/boost/boost_1_82_0

SOURCES = main.cpp
SOURCES += $(IMGUI_DIR)/backends/imgui_impl_glfw.cpp $(IMGUI_DIR)/backends/imgui_impl_opengl3.cpp
SOURCES += $(IMGUI_DIR)/imgui.cpp $(IMGUI_DIR)/imgui_draw.cpp $(IMGUI_DIR)/imgui_demo.cpp $(IMGUI_DIR)/imgui_widgets.cpp $(IMGUI_DIR)/imgui_tables.cpp

LIBS = -lGL
WEBGL_VER = -s USE_WEBGL2=1 -s USE_GLFW=3 -s FULL_ES3=1
#WEBGL_VER = USE_GLFW=2
USE_WASM = -s WASM=1

all: $(SOURCES) $(OUTPUT)

$(OUTPUT): $(SOURCES) 
	$(CXX)  $(SOURCES) -std=c++11 -o $(OUTPUT) $(LIBS) $(WEBGL_VER) -O2 --preload-file data $(USE_WASM) -I$(IMGUI_DIR) -I$(IMGUI_DIR)/backends -I$(EIGEN_DIR) -I$(BOOST_DIR)

clean:
	rm -f $(OUTPUT)
