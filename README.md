# 🧬 CRISPR Genome Editing Simulator (OpenMPI)

A simple parallelized CRISPR simulator in C using OpenMPI. This project mimics gene editing by identifying and modifying DNA sequences, leveraging distributed computation.

## 🚀 Features

- DNA sequence parsing from file  
- Parallelized editing using OpenMPI  
- Segment-based CRISPR logic  
- Example input and shell utilities

## 🛠 Usage

### 1. Compile

```bash
mpicc -o crispr_mpi crispr_mpi.c crispr_seq.c
```

### 2. Run

```bash
mpirun -np 4 ./crispr_mpi
```

## 📄 Input Format

- File: `strand.txt`  
- Example content:

```
ATGCGTACGTAGCTAGCTAGCTAGCTAGCTA
```

## 📘 Reference

Inspired by [this article](https://medium.freecodecamp.org/programming-the-genome-with-crispr-bd567a214e2a)
