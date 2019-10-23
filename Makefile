.PHONY: clean

TARGET = fbm
OBJECTS = fbm.o
CFLAGS += -g -O2 -Wall -Wextra
LDLIBS += -lfftw3 -lgsl -lgslcblas -lblas -lm

$(TARGET): $(OBJECTS)

clean:
	rm -f $(TARGET) $(OBJECTS)
