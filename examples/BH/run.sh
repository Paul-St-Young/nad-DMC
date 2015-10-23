#!/bin/bash

export OMP_NUM_THREADS=1
ROOT=../..
BIN=$ROOT/build/bin/qmcapp

mpiexec -np 1 $BIN noBO.xml
#mpiexec -np 2 $BIN BH.s002.cont.xml
