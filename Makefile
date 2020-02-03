.PHONY: all clean

TARGETS = fpt fptpb fptl fptq max maxpb maxl maxq rusage
CFLAGS += -g -O2 -Wall -Wextra
LDLIBS  = -lfftw3  -lgsl -lgslcblas -lm
LDLIBSL = -lfftw3l -lgsl -lgslcblas -lm
LDLIBSQ = -lfftw3q -lgsl -lgslcblas -lquadmath -lm

all: $(TARGETS)

clean:
	rm -f $(TARGETS) *.o

fpt.o: fbm.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -DDO_FPT -c -o $@ $<
fptpb.o: fbm.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -DDO_FPT -DDO_PHONEBOOK -c -o $@ $<
fptl.o: fbm.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -DDO_FPT -DUSE_LDBL -c -o $@ $<
fptq.o: fbm.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -DDO_FPT -DUSE_QUAD -c -o $@ $<

max.o: fbm.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -DDO_MAX -c -o $@ $<
maxpb.o: fbm.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -DDO_MAX -DDO_PHONEBOOK -c -o $@ $<
maxl.o: fbm.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -DDO_MAX -DUSE_LDBL -c -o $@ $<
maxq.o: fbm.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -DDO_MAX -DUSE_QUAD -c -o $@ $<

fpt: fpt.o
	$(CC) $(LDFLAGS) -o $@ $< $(LDLIBS)
fptpb: fptpb.o
	$(CC) $(LDFLAGS) -o $@ $< $(LDLIBS)
fptl: fptl.o
	$(CC) $(LDFLAGS) -o $@ $< $(LDLIBSL)
fptq: fptq.o
	$(CC) $(LDFLAGS) -o $@ $< $(LDLIBSQ)

max: max.o
	$(CC) $(LDFLAGS) -o $@ $< $(LDLIBS)
maxpb: maxpb.o
	$(CC) $(LDFLAGS) -o $@ $< $(LDLIBS)
maxl: maxl.o
	$(CC) $(LDFLAGS) -o $@ $< $(LDLIBSL)
maxq: maxq.o
	$(CC) $(LDFLAGS) -o $@ $< $(LDLIBSQ)

rusage: rusage.o
	$(CC) $(LDFLAGS) -o $@ $<
