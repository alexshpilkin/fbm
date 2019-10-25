.PHONY: clean

TARGET = fbm
OBJECTS = fbm.o
CFLAGS += -g -O2 -Wall -Wextra
LDLIBS += -lfftw3q -lfftw3l -lfftw3 -lgsl -lgslcblas -lblas -lquadmath -lm

$(TARGET): $(OBJECTS)

clean:
	rm -f $(TARGET) $(OBJECTS)
