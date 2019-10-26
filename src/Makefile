# -*- Makefile -*-

CFLAGS_DEBUG = -DDEBUG=3 -DTHRUST_DEBUG -G -g
NVCC = nvcc
CC = g++

CFLAGS= -Wall -std=c++11
sm30flags=-gencode arch=compute_30,code=sm_30
sm35flags=-gencode arch=compute_35,code=sm_35
sm37flags=-gencode arch=compute_37,code=sm_37
sm50flags=-gencode arch=compute_50,code=sm_50
sm52flags=-gencode arch=compute_52,code=sm_52
sm53flags=-gencode arch=compute_53,code=sm_53
sm60flags=-gencode arch=compute_60,code=sm_60
sm61flags=-gencode arch=compute_61,code=sm_510
sm70flags=-gencode arch=compute_70,code=sm_70

LIBS = -lcurand -dlink


objects = genA.o parse.o load.o

debug: $(objects)
	$(NVCC)  $(sm30flags) $(sm35flags) $(sm37flags) $(sm50flags) $(sm52flags) $(sm53flags) $(sm60flags) $(sm61flags) $(sm70flags) $(LIBS) -std=c++11 -Wno-deprecated-declarations $(CFLAGS_DEBUG) $(objects) -o ../bin/genA

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
