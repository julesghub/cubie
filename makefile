### IMPORTANT ###
# PUT YOUR METIS PATH BELOW or LEAVE COMMENTED
#metis_path=/usr/local/metis-5.1.0

ifneq ($(metis_path),) 
	LIB_PATH=${metis_path}/lib/
	INC_PATH=${metis_path}/include/
	ARGS=-I${INC_PATH} -Wl,-rpath=${LIB_PATH} -L${LIB_PATH} -lmetis
	MACRO_DEFS=-DMETIS_ENABLED=0
else
	ARGS=
	MACRO_DEFS=
endif

CC=gcc

all: cubedsphere.c
	${CC} ${MACRO_DEFS} -g cubedsphere.c -o go -lm ${ARGS}
