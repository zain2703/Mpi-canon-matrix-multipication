DEBUG = 0
CC = mpicc
LD = $(CC)
CFLAGS = -Wall # -std=c99
OPT_FLAGS = -O3 -fopt-info-vec #-ffast-math -fopt-info-vec # -march=native -g # -fopt-info-vec-missed 
#OPT_FLAGS = -Ofast -march=haswell -ffast-math -fopt-info-vec  #-fopt-info-vec-missed -g
INCLUDES=-I/usr/lib/x86_64-linux-gnu/openmpi/include
LIBS = -lm
RM = /bin/rm -f
OBJS = V4.o utils.o
EXEC = V4

ifeq ($(DEBUG), 1)
	DEFINES += -DDEBUG
	CFLAGS  += -g
else
	DEFINES += -DNDEBUG
endif

$(EXEC): $(OBJS)
	$(LD) $(LDFLAGS) $(OBJS) $(LIBS) -o $(EXEC)

utils.o: utils.c utils.h
	$(CC) $(CFLAGS) $(OPT_FLAGS) $(INCLUDES) $(DEFINES) -c utils.c

V4.o: V4.c utils.h
	$(CC) $(CFLAGS) $(OPT_FLAGS) $(INCLUDES) $(DEFINES) -c V4.c

clean:
	$(RM) $(EXEC) $(OBJS)
