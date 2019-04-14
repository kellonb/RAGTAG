# -*- Makefile -*-

CFLAGS_DEBUG = -DDEBUG=3 -DTHRUST_DEBUG -G -g
NVCC = nvcc
CC = g++

CFLAGS= -Wall -std=c++11
sm30flags=-gencode arch=compute_30,code=sm_30
LIBS = -lcurand -dlink


objects = genA.o parse.o load.o

debug: $(objects)
	$(NVCC)  $(sm30flags) $(LIBS) -std=c++11 $(CFLAGS_DEBUG) $(objects) -o ../bin/genA

all: $(objects)
	$(NVCC) $(sm30flags) $(LIBS) -std=c++11 $(objects) -o ../bin/genA

parse.o: parse.cpp
	$(CC) $(CFLAGS) -o parse.o -c $<

load.o: load.cpp
	$(CC) $(CFLAGS) -o load.o -c $<

genA.o: genA.cu 
	$(NVCC) -arch=sm_30 $(LIBS) -std=c++11 -o genA.o -c $<

clean: 
	rm -f *.o 

#.PHONY: clean all
