ifdef OPTIM
    OPTIMFLAGS+=-O3
endif

ifdef PROFILE
    PROFILEFLAGS+=-pg
endif

CC=gcc-7
LD=${CC}
CFLAGS+=-Wall -Wextra -std=c99 -fopenmp $(OPTIMFLAGS) $(PROFILEFLAGS)
LDFLAGS+= -lm

EXE=ice_model

all: ice_model
ice_model: main.o
	$(CC) $(CFLAGS) -o $(EXE) main.o $(LDFLAGS)
main.o: main.c main.h
	$(CC) $(CFLAGS) -c main.c

clean:
	rm -f $(EXE) *.o *~
