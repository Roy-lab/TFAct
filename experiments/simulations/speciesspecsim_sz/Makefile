SRC = Framework.C  SpeciesDistance.C Graph.C Vertex.C Distance.C Clique.C

CC=g++
CFLAGS = -g
INCPATH=/ahg/regev/users/sroy/thirdparty/gsl_install/include
LIBPATH=/ahg/regev/users/sroy/thirdparty/gsl_install/lib
simulateNw: $(SRC)
	$(CC) $(SRC) $(CFLAGS) -I $(INCPATH) -L $(LIBPATH) -lgsl -lgslcblas -o simulateNw

