# the executable file locates in src directory
CODE = main.cpp
EXE = Hsbg
LALIB = ./inc/

Default:
	g++ -I$(LALIB) ./src/$(CODE) -o ./src/$(EXE)

clean:
	rm -rf *.o $(EXE)
