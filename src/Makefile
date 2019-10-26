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
sm61flags=-gencode arch=compute_61,code=sm_61
sm70flags=-gencode arch=compute_70,code=sm_70

LIBS = -lcurand -std=c++11 


#nvcc -arch=sm_30 $([[ "$1" == "-d" ]] && echo -n " -DDEBUG=3 -DTHRUST_DEBUG -G -g") -o genA{,.cu} -lcurand -std=c++11
all:
	$(NVCC) genA.cu $(sm30flags) $(sm30flags) $(sm35flags) $(sm37flags) $(sm50flags) $(sm52flags) $(sm53flags) $(sm60flags) $(sm61flags) $(sm70flags) -o ../bin/genA $(LIBS)

clean: 
	rm -f *.o 

#.PHONY: clean all
