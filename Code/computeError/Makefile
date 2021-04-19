
CUDAINC   = "/Developer/GPU Computing/shared/inc"
CUDAROOT  = /usr/local/cuda
CCFLAGS   = -I. -Iinclude -O3 -I/usr/local/cuda/include
GLLIB     = -lGL -lGLU -lglut -lGLEW
CUDAGEN   = -gencode=arch=compute_30,code=\"sm_30,compute_30\"
#--------------------------------------------------------------------------------------------------------------------------------

CXX  = g++
NVCC = $(CUDAROOT)/bin/nvcc

OBJECTS = main.o vis.o texture2D.o lodepng.o

all: skeleton

skeleton: $(OBJECTS)
	$(CXX) -o skeleton $(OBJECTS) $(CCFLAGS) $(GLLIB) -L$(CUDAROOT)/targets/x86_64-linux/lib -lcudart

main.o: main.cpp
	$(CXX) -c $< $(CCFLAGS) -FGLUT -FOpenGL

vis.o: vis.cpp
	$(CXX) -c $< $(CCFLAGS)

texture2D.o: texture2D.cpp
	$(CXX) -c $< $(CCFLAGS)
	
lodepng.o: lodepng.cpp
	$(CXX) -c $< $(CCFLAGS)
	
clean:
	-rm *.o skeleton