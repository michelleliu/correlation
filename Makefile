#
# Makefile for correlation.c
#

CC = g++
LIBS = -lm
CFLAGS = -std=c++11

TARGETS = correlation

all: $(TARGETS)

correlation: correlation.o
	$(CC) -o $@ $(LIBS) $<

correlation.o: correlation.cpp correlation.h
	$(CC) -c $(CFLAGS) correlation.cpp

clean:
	rm -f *.o $(TARGETS)
