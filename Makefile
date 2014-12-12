DIR_BUILD = ../Mie_eff_build
DIR_SRC = ./
DIR_OBJ = ${DIR_BUILD}/obj
DIR_BIN = ${DIR_BUILD}/

SRC = $(wildcard ${DIR_SRC}/*.for)  
OBJ = $(patsubst %.for,${DIR_OBJ}/%.o,$(notdir ${SRC})) 

TARGET = Mie_sphere

BIN_TARGET = ${DIR_BIN}/${TARGET}

FC = gfortran
FCFLAGS = -O3 -g -Wall 

all:buildrepo ${BIN_TARGET}

${BIN_TARGET}:${OBJ}
	$(FC) $(OBJ)  -o $@
    
${DIR_OBJ}/%.o:${DIR_SRC}/%.for
	$(FC) $(FCFLAGS) -c  $< -o $@
.PHONY:clean
clean:
	find ${DIR_OBJ} -name *.o *.mod -exec rm -rf {}

buildrepo:
	mkdir -p ${DIR_BUILD} ${DIR_BIN} ${DIR_OBJ}
