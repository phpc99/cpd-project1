#!/bin/bash

mkdir -p ./results
mkdir -p ./bin

rm -f ./bin/programm
g++ -fopenmp -O2 matrixproduct.cpp -o ./bin/programm -lpapi

echo "[INFO]  Running OnMult 600 - 3000"
for SIZE in $(seq 600 400 3000); do
    ./bin/programm $SIZE $SIZE $SIZE 1 "results/cpp.csv"
done

echo "[INFO]  Running OnMultLine 600 - 3000"
for SIZE in $(seq 600 400 3000); do
    ./bin/programm $SIZE $SIZE $SIZE 2 "results/cpp.csv"
done

echo "[INFO]  Running OnMultLine 4096 - 10240"
for SIZE in $(seq 4096 2048 10240); do
    ./bin/programm $SIZE $SIZE $SIZE 2 "results/cpp.csv"
done

echo "[INFO]  Running OnMultBlock 4096 - 10240"
for SIZE in $(seq 4096 2048 10240); do
    for BLOCK_SIZE in 128 256 512; do
        ./bin/programm $SIZE $SIZE $SIZE 3 "results/cpp.csv" $BLOCK_SIZE
    done
done

echo "[INFO]  C++ tests finished"

echo "[INFO]  Running Python OnMult 600 - 3000"
for SIZE in $(seq 600 400 3000); do
    python3 "matrixproduct.py" $SIZE $SIZE $SIZE 1 "results/python.csv"
done

echo "[INFO]  Running Python OnMultLine 600 - 3000"
for SIZE in $(seq 600 400 3000); do
    python3 "matrixproduct.py" $SIZE $SIZE $SIZE 2 "results/python.csv"
done

echo "[INFO]  Python tests finished"

echo "[INFO]  Running OnMultLineParallelV1 600 - 3000"
for SIZE in $(seq 4096 2048 10240); do
    ./bin/programm $SIZE $SIZE $SIZE 4 "results/cpp.csv"
done

echo "[INFO]  Running OnMultLineParallelV2 600 - 3000"
for SIZE in $(seq 4096 2048 10240); do
    ./bin/programm $SIZE $SIZE $SIZE 5 "results/cpp.csv"
done

echo "[INFO]  Finish Parallel tests"
