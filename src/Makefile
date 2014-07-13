#
# Makefile for correlation.c
#

CC = g++
LIBS = -lm
CFLAGS = -std=c++11 -pedantic -Wall

TARGETS = correlation dipole hbonds #correlation-r0 # test correlation-window correlation-tinker

all: $(TARGETS)

correlation-tinker: correlation-tinker.o
	$(CC) -o $@ $(LIBS) $<
	mv $@ ../
correlation-tinker.o: correlation-tinker.cpp
	$(CC) -c $(CFLAGS) correlation-tinker.cpp

correlation-window: correlation-window.o
	$(CC) -o $@ $(LIBS) $<
	mv $@ ../
correlation-window.o: correlation-window.cpp
	$(CC) -c $(CFLAGS) correlation-window.cpp

hbonds: hbonds.o common.o
	$(CC) -o $@ $(LIBS) $< common.o
	mv $@ ../
hbonds.o: hbonds.cpp common.h
	$(CC) -c $(CFLAGS) hbonds.cpp

map: map.o
	$(CC) -o $@ $(LIBS) $< common.o
	mv $@ ../
map.o: map.cpp common.h
	$(CC) -c $(CFLAGS) map.cpp

bin: bin.o
	$(CC) -o $@ $(LIBS) $< common.o
	mv $@ ../
bin.o: bin.cpp common.h
	$(CC) -c $(CFLAGS) bin.cpp

dipole: dipole.o
	$(CC) -o $@ $(LIBS) $< common.o
	mv $@ ../
dipole.o: dipole.cpp common.h
	$(CC) -c $(CFLAGS) dipole.cpp

correlation: correlation.o common.o
	$(CC) -o $@ $(LIBS) $< common.o
	mv $@ ../
correlation.o: correlation.cpp common.h
	$(CC) -c $(CFLAGS) correlation.cpp

correlation-r0: correlation-r0.o
	$(CC) -o $@ $(LIBS) $<
	mv $@ ../
correlation-r0.o: correlation-r0.cpp
	$(CC) -c $(CFLAGS) correlation-r0.cpp

test: test.o
	$(CC) -o $@ $(LIBS) $< common.o
	mv $@ ../
test.o: test.cpp common.h
	$(CC) -c $(CFLAGS) test.cpp

common.o: common.cpp common.h
	$(CC) -c $(CFLAGS) common.cpp

clean:
	rm -f *.o
	cd ../
	rm -f $(TARGETS)
