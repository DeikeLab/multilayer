#SET FILES OF INTEREST 
SRC_FILES = $(addsuffix .c,$(OUTPUT))
MPI_FILES = $(addprefix _,$(SRC_FILES))
OBJ_FILES = $(addsuffix .o,$(basename $(SRC_FILES)))
OUT_FOLDER = $(addprefix ./f,$(OUTPUT))

#SET C FLAGS AND BASILISK COMPILER, LINKING FLAGS
CC = gcc
CFLAGS += -Wall -O2 -DLAYERS=1 #-fopenmp
BAS = qcc
LDFLAGS = -L${BASILISK}/ppr -lppr -lgfortran -lm 
OPENGLIBS = -L/usr/lib/gcc/x86_64-linux-gnu/7 -L${BASILISK}/gl -lfb_osmesa -lglutils -lGLU -lOSMesa
#LDFLAGS = -L/usr/local/MATLAB/R2018b/sys/os/glnxa64 -L/usr/share/doc -L/usr/lib/gcc/x86_64-linux-gnu/7 -L/usr/lib/x86_64-linux-gnu -lgfortran -lppr -lm -L${BASILISK}/ppr 
#OPENGLIBS = -L${BASILISK}/gl -lfb_osmesa -lglutils -lGLU -lOSMesa
OPENGLINC = ${BASILISK}

#MAKE COMMAND: TAKES ARGUMENT AND MAKES EXECUTABLE IN DESIRED FOLDER
$(OUTPUT): 
	$(BAS) $(CFLAGS) $(SRC_FILES) -I$(OPENGLINC) -o $(OUT_FOLDER)/$(OUTPUT) $(OPENGLIBS) $(LDFLAGS)
