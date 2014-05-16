#
# Makefile for correlation.c
#

CC = g++
LIBS = -lm

TARGETS = correlation

all: $(TARGETS)

correlation: correlation.o
	$(CC) -o $@ $(LIBS) $<

correlation.o: correlation.cpp correlation.h
	$(CC) -c $(CFLAGS) correlation.cpp

clean:
	rm -f *.o $(TARGETS)
