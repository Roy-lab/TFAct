# INCLUDE_PATH = ./include
# LIBRARY_PATH = ./lib
LFLAG = -lgsl -lgslcblas -lrt -lpthread -lgfortran
CONDA_PREFIX = /mnt/dv/wid/projects7/Roy-singlecell2/yuda_work/envs/IBPG
INCLUDE_PATH = ${CONDA_PREFIX}/include
LIBRARY_PATH = ${CONDA_PREFIX}/lib
CC=g++
CFLAGS = -g -DHAVE_INLINE -O3 -fopenmp

BIN=tetraNMF
# SRC=EdgeList.C Error.C EvidenceManager.C Framework.C Graph.C main.C Matrix.C MemoryCheck.C NCA.C NCALearner.C Variable.C VariableManager.C VariableSelection.C

$(BIN): modules/NCAFramework/*.C modules/nmf/*.cpp modules/nmf/random_svd/*.o
	$(CC) $^ $(CFLAGS) -o $(BIN) -L${LIBRARY_PATH} $(LFLAG) -I${INCLUDE_PATH}

clean:
	rm -f $(BIN) *~
