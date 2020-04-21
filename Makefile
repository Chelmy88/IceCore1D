ifdef OPTIM
    OPTIMFLAGS+=-Ofast -march=skylake
endif
ifdef PROFILE
    PROFILEFLAGS+=-p
endif

ifeq ($(PREFIX),)
    PREFIX := /usr/local/bin
endif

OPTIMFLAGS?=-O0 -Wall -Wextra

CC= gcc-9
LD=${CC}
CFLAGS+= -std=c99 -fopenmp $(OPTIMFLAGS) $(PROFILEFLAGS)
LDFLAGS+= -lm

SOURCEDIR=src
BUILDIR=build
SOURCES = $(wildcard $(SOURCEDIR)/*.c)
OBJECTS = $(patsubst %.c, %.o, $(SOURCES))
OBJECTS:= $(subst $(SOURCEDIR), $(BUILDIR), $(OBJECTS))

BINARY=ice_model
BINARYDIR=bin

all: directories $(BINARY)

directories: $(BUILDIR) $(BINARYDIR)

$(BUILDIR) $(BINARYDIR):
	mkdir -p $@

$(BINARY): $(OBJECTS)
	$(CC) $(CFLAGS) $(OBJECTS) $(LDFLAGS) -o $(BINARYDIR)/$(BINARY)

$(BUILDIR)/%.o: $(SOURCEDIR)/%.c
	$(CC) $(CFLAGS) -c $< -o $@

install:
	install -m 775 $(BINARYDIR)/$(BINARY) $(PREFIX)

clean:
	rm -f $(BINARYDIR)/* $(BUILDIR)/*
