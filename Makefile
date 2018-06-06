# the executable file locates in src directory
CODE = main.cpp
EXE = Hsbg
LALIB = ./lib/

Default:
	g++ -I$(LALIB) ./src/$(CODE) -o ./src/$(EXE)

clean:
	rm -rf *.o $(EXE)
