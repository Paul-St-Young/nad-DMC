#!/bin/bash

export OMP_NUM_THREADS=1
ROOT=../..
BIN=$ROOT/branch/atom/build/bin/qmcapp

$BIN noBO.xml
