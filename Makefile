CXX = emcc
OUTPUT = index.html
IMGUI_DIR:=C:/cpp/imgui/imgui
IMPLOT_DIR:=C:/cpp/implot
EIGEN_DIR:=C:/cpp/eigen/eigen-3.4.0
BOOST_DIR:=C:/cpp/boost/boost_1_82_0
INCLUDE:=./include

SOURCES = main.cpp
SOURCES += $(IMGUI_DIR)/backends/imgui_impl_glfw.cpp $(IMGUI_DIR)/backends/imgui_impl_opengl3.cpp
SOURCES += $(IMGUI_DIR)/imgui.cpp $(IMGUI_DIR)/imgui_draw.cpp $(IMGUI_DIR)/imgui_widgets.cpp $(IMGUI_DIR)/imgui_tables.cpp
SOURCES += $(IMGUI_DIR)/imgui_spectrum.cpp 
SOURCES += ./src/sequenceperiod.cpp ./src/multiColorButton.cpp  ./src/design.cpp
SOURCES += ./src/statisticalmodel.cpp ./src/menubar.cpp ./src/glmmmodel.cpp ./src/modelupdater.cpp ./src/plotdata.cpp 
SOURCES += ./src/modelchecker.cpp ./src/krigingdata.cpp  
SOURCES += ./src/applog.cpp ./src/logger.cpp ./src/appmodel.cpp 
SOURCES += $(IMPLOT_DIR)/implot.cpp $(IMPLOT_DIR)/implot_items.cpp
# ./src/results.cpp ./src/krigger.cpp ./src/datasimulate.cpp 
LIBS = -lGL
WEBGL_VER = -s USE_WEBGL2=1 -s USE_GLFW=3 -s FULL_ES3=1 
USE_WASM = -s WASM=1 -s ALLOW_MEMORY_GROWTH=1 -s STACK_SIZE=15MB 
CPPFLAGS += -I$(BOOST_DIR) -I$(EIGEN_DIR) -I$(INCLUDE) -fexperimental-library
LDFLAGS += --shell-file shell_full.html

# for debugging add: -s NO_DISABLE_EXCEPTION_CATCHING -g and remove optimisation flags -Os -g2 -fno-math-errno -fno-math-errno -Os

all: $(SOURCES) $(OUTPUT)

$(OUTPUT): $(SOURCES) 
	$(CXX)  $(SOURCES) -std=c++20 -o $(OUTPUT) $(LIBS) $(WEBGL_VER) -s NO_DISABLE_EXCEPTION_CATCHING -Os -g2 --preload-file data $(USE_WASM) -I$(IMGUI_DIR) -I$(IMGUI_DIR)/backends $(LDFLAGS) $(CPPFLAGS) 
clean:
	rm -f $(OUTPUT)