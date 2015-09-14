#!/bin/bash

export OMP_NUM_THREADS=1
ROOT=../..
BIN=$ROOT/build/bin/qmcapp

$BIN noBO.xml
