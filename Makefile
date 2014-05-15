#
# Makefile for correlation.c
#

CC = gcc
LIBS = -lm

TARGETS = correlation

all: $(TARGETS)

correlation: correlation.o
	$(CC) -o $@ $(LIBS) correlation.o

correlation.o: correlation.c correlation.h
	$(CC) -c $(CFLAGS) correlation.c

clean:
	rm -f *.o $(TARGETS)
