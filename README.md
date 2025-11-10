# Inverting-GF2-matrices

![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)  
![C++17](https://img.shields.io/badge/C%2B%2B-17-blue.svg)

**Author:** Ibrahim Mammadov  
**Contact:** ibrahim.22@intl.zju.edu.cn

A suite of C++ tools for performing advanced matrix operations over the finite field GF(2). This toolkit provides algorithms for fast matrix inversion, guaranteed matrix repair, and the generation of various types of test matrices.

---

## Project Status and Philosophy

Please note that the code in this repository is provided as a proof‑of‑concept and a practical demonstration of the described algorithms. The primary goal is to serve as a clear, functional reference for these advanced matrix operations.

While the algorithms themselves are robust, the code is structured as simple, self‑contained programs for ease of compilation and understanding. A production‑grade implementation would involve a more modular, library‑based architecture (with separate `.hpp`/`.cpp` files), the use of a build system like CMake, and potentially parallelization for enhanced performance.

This repository successfully demonstrates the functionality and correctness of the core ideas. It is an excellent starting point for academic exploration, testing, and further development. Feel free to contact me if you have ideas that may improve the code.

---

## Project Structure

```text
.
├── Binary_Matrix_Inversion/
│   ├── main.cpp
│   └── strassen_inversion.hpp
├── Generating_Invertible_Matrices/
│   ├── main.cpp
│   └── generate_all_invertible_matrices.hpp
├── Repairing_Singular_Matrices/
│   ├── main.cpp
│   └── matrix_repair.hpp
└── README.md

```
---

# Building the Code

The tools are single‑file programs. You can compile them directly from the command line with a C++17‑compliant compiler (e.g., `g++`).

## Binary Matrix Inversion
cd Binary_Matrix_Inversion/

g++ -std=c++17 -O2 -o invert main.cpp

## Generating Invertible Matrices
cd ../Generating_Invertible_Matrices/

g++ -std=c++17 -O2 -o gen_all main.cpp

## Repairing Singular Matrices
cd ../Repairing_Singular_Matrices/

g++ -std=c++17 -O2 -o gen_sing main.cpp


## Workflows and Usage
# 1. Binary Matrix Inversion
This workflow demonstrates generating an invertible matrix and then finding its inverse.

cd Binary_Matrix_Inversion/

./gen_inv            # Enter the size of the matrix: 7

Creates matrix.txt containing a random, invertible 7×7 matrix

./invert             # Reads matrix.txt, computes its inverse, saves to answer.txt

version 0.9.0: Beta version might have given late responses in edge cases

version 1.0.0: Base version where it can give consistent result by switching strassen and gauss-jordan accoding to matrix

version 1.0.1: Multithread (openMP) application which speeds up process up to 4 times.
# 2. Generating All Invertible Matrices
This tool exhaustively generates all unique invertible matrices for a given (small) size.

cd Generating_Invertible_Matrices/
./gen_all            # Enter matrix size n (recommended ≤ 4): 3

version 1.0.0: Base code with stable input and output interaction

version 1.1.0: On-fly generation: instead of saving all generated matrixes, save just last matrix and generate out of that last matrix

version 1.1.1: Multithread (openMP) application which speeds up process up to 4 times.

# Computes all 168 unique invertible 3×3 matrices and saves to unique_GF2_matrices.txt

# 3. Repairing Singular Matrices
This workflow demonstrates generating a singular matrix and then repairing it to make it invertible.

version 1.0.0: Base code with stable input and output interaction

version 1.0.1: Multithread (openMP) application which extremely fast computational power up to the size n=2048.

cd Repairing_Singular_Matrices/

./gen_sing           # Enter matrix size n: 7

Creates matrix.txt containing a random, singular 7×7 matrix

./repair             # Applies the guaranteed repair algorithm, saves to answer.txt








| File                                   | Version | Purpose                                                       | Input              | Output                    |
| -------------------------------------- | ------- | ------------------------------------------------------------- | ------------------ | ------------------------- |
| `strassen_inversion.cpp`       | 1.0.1   | Inverts a matrix using a pivoted Strassen‑based algorithm     | User prompt or `matrix.txt` | `answer.txt`              |
| `generate_all_invertible_matrices.cpp` | 1.1.1   | Exhaustively generates all unique invertible matrices for n   | User prompt (size) | `unique_GF2_matrices.txt` |
| `matrix_repair.cpp`         | 1.0.1   | Repairs a singular matrix using a provably correct algorithm  | `matrix.txt` or User prompt       | `answer.txt`              |

