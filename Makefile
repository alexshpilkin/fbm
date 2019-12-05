.PHONY: all clean

TARGETS = fpt fpt-pb max max-pb rusage
OBJECTS = fpt.o fpt-pb.o max.o max-pb.o rusage.o
CFLAGS += -g -O2 -Wall -Wextra
LDLIBS += -lfftw3q -lfftw3l -lfftw3 -lgsl -lgslcblas -lblas -lquadmath -lm

all: $(TARGETS)

clean:
	rm -f $(TARGETS) $(OBJECTS)

fpt.o fpt-pb.o max.o max-pb.o: fbm.c
