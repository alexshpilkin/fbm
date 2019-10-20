.PHONY: clean

TARGET = dh
OBJECTS = dh.o
CFLAGS += -g -O2
LDLIBS += -lfftw3 -lgsl -lgslcblas -lblas -lm

$(TARGET): $(OBJECTS)

clean:
	rm -f $(TARGET) $(OBJECTS)
