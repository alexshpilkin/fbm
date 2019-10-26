.PHONY: all clean

TARGETS = fbm rusage
OBJECTS = fbm.o rusage.o
CFLAGS += -g -O2 -Wall -Wextra
LDLIBS += -lfftw3q -lfftw3l -lfftw3 -lgsl -lgslcblas -lblas -lquadmath -lm

all: $(TARGETS)

clean:
	rm -f $(TARGETS) $(OBJECTS)
