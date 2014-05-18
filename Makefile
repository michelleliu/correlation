#
# Makefile for correlation.c
#

CC = g++
LIBS = -lm
CFLAGS = -std=c++11

TARGETS = correlation correlation-window correlation-tinker test

all: $(TARGETS)

correlation-tinker: correlation-tinker.o
	$(CC) -o $@ $(LIBS) $<
correlation-tinker.o: correlation-tinker.cpp correlation.h
	$(CC) -c $(CFLAGS) correlation-tinker.cpp

correlation-window: correlation-window.o
	$(CC) -o $@ $(LIBS) $<
correlation-window.o: correlation-window.cpp correlation.h
	$(CC) -c $(CFLAGS) correlation-window.cpp

correlation: correlation.o
	$(CC) -o $@ $(LIBS) $<
correlation.o: correlation.cpp correlation.h
	$(CC) -c $(CFLAGS) correlation.cpp

test: test.o
	$(CC) -o $@ $(LIBS) $<
test.o: test.cpp
	$(CC) -c $(CFLAGS) test.cpp


clean:
	rm -f *.o $(TARGETS)
