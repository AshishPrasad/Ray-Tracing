all: assign3.o Ray.o Vector4.o Matrix4.o
	g++ assign3.o Ray.o Vector4.o Matrix4.o -lglut -lGLU -lGL -o ashish.out

assign3: assign3.cpp ImageDataStructure.h
	g++ -c assign3.cpp -lglut -lGLU -lGL
Ray: Ray.cpp Ray.h
	g++ -c Ray.cpp -lglut -lGLU -lGL
Vector4: Vector4.cpp Vector4.h
	g++ -c Vector4.cpp -lglut -lGLU -lGL
Matrix4: Matrix4.cpp Matrix4.h
	g++ -c Matrix4.cpp -lglut -lGLU -lGL

clean:
	rm -rf *.o ashish.out
