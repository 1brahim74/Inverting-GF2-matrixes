# Inverting-GF2-matrices

![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)  
![C++17](https://img.shields.io/badge/C%2B%2B-17-blue.svg)

**Author:** Ibrahim Mammadov  
**Contact:** ibrahim.22@intl.zju.edu.cn

A suite of C++ tools for performing advanced matrix operations over the finite field GF(2). This toolkit provides algorithms for fast matrix inversion, guaranteed matrix repair, and the generation of various types of test matrices.

---

## Project Status and Philosophy

Please note that the code in this repository is provided as a proof‑of‑concept and a practical demonstration of the described algorithms. The primary goal is to serve as a clear, functional reference for these advanced matrix operations.

While the algorithms themselves are robust, the code is structured as simple, self‑contained programs for ease of compilation and understanding. A production‑grade implementation would involve a more modular, library‑based architecture (with separate `.h`/`.cpp` files), the use of a build system like CMake, and potentially parallelization for enhanced performance.

This repository successfully demonstrates the functionality and correctness of the core ideas. It is an excellent starting point for academic exploration, testing, and further development. Feel free to contact me if you have ideas that may improve the code.

---

## Project Structure

.
├── Binary_Matrix_Inversion/
│ ├── Invertible_matrix_generator.cpp
│ └── strassen_inversion_pivoted.cpp
│
├── Generating_Invertible_Matrices/
│ └── generate_all_invertible_matrices.cpp
│
├── Repairing_Singular_Matrices/
│ ├── singular_matrix_generator.cpp
│ └── matrix_repair_guaranteed.cpp
│
└── README.md


---

## Building the Code

The tools are single‑file programs. You can compile them directly from the command line with a C++17‑compliant compiler (e.g., `g++`).

# Binary Matrix Inversion
cd Binary_Matrix_Inversion/

g++ -std=c++17 -O2 -o gen_inv Invertible_matrix_generator.cpp

g++ -std=c++17 -O2 -o invert strassen_inversion_pivoted.cpp

# Generating Invertible Matrices
cd ../Generating_Invertible_Matrices/

g++ -std=c++17 -O2 -o gen_all generate_all_invertible_matrices.cpp

# Repairing Singular Matrices
cd ../Repairing_Singular_Matrices/

g++ -std=c++17 -O2 -o gen_sing singular_matrix_generator.cpp

g++ -std=c++17 -O2 -o repair matrix_repair_guaranteed.cpp


## Workflows and Usage
# 1. Binary Matrix Inversion
This workflow demonstrates generating an invertible matrix and then finding its inverse.

cd Binary_Matrix_Inversion/

./gen_inv            # Enter the size of the matrix (must be a power of 2): 8

Creates matrix.txt containing a random, invertible 8×8 matrix

./invert             # Reads matrix.txt, computes its inverse, saves to answer.txt

# 2. Repairing Singular Matrices
This workflow demonstrates generating a singular matrix and then repairing it to make it invertible.


cd Repairing_Singular_Matrices/

./gen_sing           # Enter matrix size n: 8

Creates matrix.txt containing a random, singular 8×8 matrix

./repair             # Applies the guaranteed repair algorithm, saves to answer.txt

# 3. Generating All Invertible Matrices
This tool exhaustively generates all unique invertible matrices for a given (small) size.

cd Generating_Invertible_Matrices/
./gen_all            # Enter matrix size n (recommended ≤ 4): 3
# Computes all 168 unique invertible 3×3 matrices and saves to unique_GF2_matrices.txt






| File                                   | Version | Purpose                                                       | Input              | Output                    |
| -------------------------------------- | ------- | ------------------------------------------------------------- | ------------------ | ------------------------- |
| `strassen_inversion_pivoted.cpp`       | 0.9.0   | Inverts a matrix using a pivoted Strassen‑based algorithm     | `matrix.txt`       | `answer.txt`              |
| `Invertible_matrix_generator.cpp`      | 1.0.0   | Generates a guaranteed‑invertible matrix via LU decomposition | User prompt (size) | `matrix.txt`              |
| `singular_matrix_generator.cpp`        | 1.0.0   | Generates a guaranteed‑singular matrix by checking the rank   | User prompt (size) | `matrix.txt`              |
| `matrix_repair_guaranteed.cpp`         | 1.0.0   | Repairs a singular matrix using a provably correct algorithm  | `matrix.txt`       | `answer.txt`              |
| `generate_all_invertible_matrices.cpp` | 1.0.0   | Exhaustively generates all unique invertible matrices for n   | User prompt (size) | `unique_GF2_matrices.txt` |
