#!/bin/bash

## compile
g++ -O3 -std=c++11 -Wno-deprecated -o bin/mainDynamic src/main.cpp

## execute
EXEC="./bin/mainDynamic"

VF3DYNAMIC_TARGET="example-graphs/directed-graphs/target-graph-example"
NB_TIMESTAMPS_TARGET=8

VF3DYNAMIC_PATTERN="example-graphs/directed-graphs/pattern-graph-example"
NB_TIMESTAMPS_PATTERN=5

WINDOW_SIZE_INCREMENT=0

$EXEC $VF3DYNAMIC_TARGET $NB_TIMESTAMPS_TARGET $VF3DYNAMIC_PATTERN $NB_TIMESTAMPS_PATTERN $WINDOW_SIZE_INCREMENT

echo "Done!"
