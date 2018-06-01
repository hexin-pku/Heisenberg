#
CODE = main.cc
EXE = scf

Default:
    g++ -I../inc/Eigen ./src/$(CODE) -o $(EXE)

clean:
    rm -rf *.o $(EXE)
