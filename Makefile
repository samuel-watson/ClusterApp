CXX = emcc
OUTPUT = index.html
IMGUI_DIR:=C:/cpp/imgui/imgui
EIGEN_DIR:=C:/cpp/eigen/eigen-3.4.0
BOOST_DIR:=C:/cpp/boost/boost_1_82_0

SOURCES = main.cpp
SOURCES += $(IMGUI_DIR)/backends/imgui_impl_glfw.cpp $(IMGUI_DIR)/backends/imgui_impl_opengl3.cpp
SOURCES += $(IMGUI_DIR)/imgui.cpp $(IMGUI_DIR)/imgui_draw.cpp $(IMGUI_DIR)/imgui_demo.cpp $(IMGUI_DIR)/imgui_widgets.cpp $(IMGUI_DIR)/imgui_tables.cpp
SOURCES += ./src/clusterapp.cpp ./src/sequenceperiod.cpp ./src/multiColorButton.cpp  ./src/design.cpp
SOURCES += ./src/statisticalmodel.cpp ./src/glmmmodel.cpp ./src/modelupdater.cpp ./src/modelchecker.cpp

LIBS = -lGL
WEBGL_VER = -s USE_WEBGL2=1 -s USE_GLFW=3 -s FULL_ES3=1 
USE_WASM = -s WASM=1 -s ALLOW_MEMORY_GROWTH=1  -s NO_DISABLE_EXCEPTION_CATCHING
CPPFLAGS += -I$(BOOST_DIR) -I$(EIGEN_DIR) 
LDFLAGS += --shell-file shell_full.html

# for debigging add: -s NO_DISABLE_EXCEPTION_CATCHING and remove optimisation flag

all: $(SOURCES) $(OUTPUT)

$(OUTPUT): $(SOURCES) 
	$(CXX)  $(SOURCES) -std=c++17 -o $(OUTPUT) $(LIBS) $(WEBGL_VER) --preload-file data $(USE_WASM) -I$(IMGUI_DIR) -I$(IMGUI_DIR)/backends $(LDFLAGS) $(CPPFLAGS)
clean:
	rm -f $(OUTPUT)

