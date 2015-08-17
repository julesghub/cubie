### IMPORTANT ###
# PUT YOUR METIS PATH BELOW
PATH_METIS=/usr/local/metis-5.1.0

CC=gcc
LIB_PATH=${PATH_METIS}/lib/
INC_PATH=${PATH_METIS}/include/

all: cubedsphere.c
	${CC} -g cubedsphere.c -o go -I${INC_PATH} -Wl,-rpath=${LIB_PATH} -L${LIB_PATH} -lmetis -lm
