lorenz:	src/drawing.cpp src/main.cpp src/sim.cpp src/solver.cpp src/vbo.cpp src/vector.cpp include/*.h
	clang++ -O3 -o lorenz src/drawing.cpp src/main.cpp src/sim.cpp src/solver.cpp src/vbo.cpp src/vector.cpp -w -I include -framework OpenGL -framework GLUT

clean:
	rm -rf lorenz
