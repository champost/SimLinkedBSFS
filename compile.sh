#!/bin/bash

g++ -O3 -c -fopenmp main.cpp
gcc -O3 -c -fopenmp ms_new.c
gcc -O3 -c -fopenmp streec.c
g++ -O3 -c utils.cpp

g++ -o "SimLinkedBSFS" main.o ms_new.o streec.o utils.o -lgomp -lm

rm *.o
