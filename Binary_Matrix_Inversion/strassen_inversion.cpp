//
// Copyright (c) 2023, Ibrahim Mammadov
// Refactored by an AI assistant to support arbitrary matrix sizes.
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
//

/**
 * @file strassen_inversion_dynamic.cpp
 * @author Ibrahim Mammadov (Original), AI Assistant (Refactoring)
 * @contact ibrahim.22@intl.zju.edu.cn
 * @brief Implements a pivoted, recursive block matrix inversion for matrices of
 *        ANY size over GF(2), using Strassen's algorithm for multiplication.
 * @version 1.0.0
 * This program reads a square binary matrix of any dimension from "matrix.txt".
 * If the dimension is not a power of two, it pads the matrix to the next
 * power of two. It then computes the inverse using a recursive method and
 * Strassen's algorithm. The final result (un-padded) is written to "answer.txt".
 */

#include <iostream>
#include <vector>
#include <optional>
#include <numeric>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>

// --- Type Definitions and Constants ---

// Use std::vector<bool> which is a space-optimized dynamic bitset
using MatrixRow = std::vector<bool>;
using Matrix = std::vector<MatrixRow>;

// The recursion cutoff size. Below this, use Gauss-Jordan.
const size_t STRASSEN_CUTOFF = 16; 

// --- Forward Declarations ---

std::optional<Matrix> strassenInverse(const Matrix& A);
Matrix strassenMultiply(const Matrix& A, const Matrix& B);

// --- Utility and Helper Functions ---

void printMatrix(const Matrix& A, std::ostream& out = std::cout) {
    if (A.empty()) {
        out << "[Empty Matrix]\n";
        return;
    }
    for (const auto& row : A) {
        for (bool val : row) {
            out << val << " ";
        }
        out << "\n";
    }
}

Matrix add(const Matrix& A, const Matrix& B) {
    size_t n = A.size();
    Matrix C(n, MatrixRow(n));
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            C[i][j] = A[i][j] ^ B[i][j];
        }
    }
    return C;
}

// In GF(2), subtraction is the same as addition (XOR)
Matrix subtract(const Matrix& A, const Matrix& B) {
    return add(A, B);
}

void split(const Matrix& A, Matrix& A11, Matrix& A12, Matrix& A21, Matrix& A22) {
    size_t k = A.size() / 2;
    for (size_t i = 0; i < k; ++i) {
        for (size_t j = 0; j < k; ++j) {
            A11[i][j] = A[i][j];
            A12[i][j] = A[i][j + k];
            A21[i][j] = A[i + k][j];
            A22[i][j] = A[i + k][j + k];
        }
    }
}

Matrix combine(const Matrix& A11, const Matrix& A12, const Matrix& A21, const Matrix& A22) {
    size_t k = A11.size();
    size_t n = 2 * k;
    Matrix A(n, MatrixRow(n));
    for (size_t i = 0; i < k; ++i) {
        for (size_t j = 0; j < k; ++j) {
            A[i][j] = A11[i][j];
            A[i][j + k] = A12[i][j];
            A[i + k][j] = A21[i][j];
            A[i + k][j + k] = A22[i][j];
        }
    }
    return A;
}

// --- Core Algorithms ---

std::optional<Matrix> gaussJordanInverse(Matrix A) {
    size_t n = A.size();
    if (n == 0) return Matrix{};
    
    Matrix I(n, MatrixRow(n, false));
    for (size_t i = 0; i < n; ++i) I[i][i] = true;

    for (size_t j = 0; j < n; ++j) {
        size_t pivot_row = j;
        while (pivot_row < n && !A[pivot_row][j]) {
            pivot_row++;
        }
        if (pivot_row == n) return std::nullopt; // Singular
        std::swap(A[j], A[pivot_row]);
        std::swap(I[j], I[pivot_row]);
        for (size_t i = 0; i < n; ++i) {
            if (i != j && A[i][j]) {
                for (size_t k = 0; k < n; ++k) {
                    A[i][k] = A[i][k] ^ A[j][k];
                    I[i][k] = I[i][k] ^ I[j][k];
                }
            }
        }
    }
    return I;
}

Matrix strassenMultiply(const Matrix& A, const Matrix& B) {
    size_t n = A.size();
    if (n <= STRASSEN_CUTOFF) { // Base case: standard multiplication
        Matrix C(n, MatrixRow(n, false));
        for (size_t i = 0; i < n; ++i) {
            for (size_t j = 0; j < n; ++j) {
                bool sum = false;
                for (size_t k = 0; k < n; ++k) {
                    sum = sum ^ (A[i][k] & B[k][j]);
                }
                C[i][j] = sum;
            }
        }
        return C;
    }

    size_t k = n / 2;
    Matrix A11(k, MatrixRow(k)), A12(k, MatrixRow(k)), A21(k, MatrixRow(k)), A22(k, MatrixRow(k));
    Matrix B11(k, MatrixRow(k)), B12(k, MatrixRow(k)), B21(k, MatrixRow(k)), B22(k, MatrixRow(k));

    split(A, A11, A12, A21, A22);
    split(B, B11, B12, B21, B22);

    Matrix P1 = strassenMultiply(add(A11, A22), add(B11, B22));
    Matrix P2 = strassenMultiply(add(A21, A22), B11);
    Matrix P3 = strassenMultiply(A11, subtract(B12, B22));
    Matrix P4 = strassenMultiply(A22, subtract(B21, B11));
    Matrix P5 = strassenMultiply(add(A11, A12), B22);
    Matrix P6 = strassenMultiply(subtract(A21, A11), add(B11, B12));
    Matrix P7 = strassenMultiply(subtract(A12, A22), add(B21, B22));

    Matrix C11 = add(subtract(add(P1, P4), P5), P7);
    Matrix C12 = add(P3, P5);
    Matrix C21 = add(P2, P4);
    Matrix C22 = add(subtract(add(P1, P3), P2), P6);

    return combine(C11, C12, C21, C22);
}

std::optional<Matrix> strassenInverse(const Matrix& A) {
    size_t n = A.size();
    if (n == 0) return Matrix{};
    if (n <= STRASSEN_CUTOFF) {
        return gaussJordanInverse(A);
    }
    // This function now assumes N is a power of 2, handled by the wrapper.

    size_t k = n / 2;
    Matrix A11(k, MatrixRow(k)), A12(k, MatrixRow(k)), A21(k, MatrixRow(k)), A22(k, MatrixRow(k));
    split(A, A11, A12, A21, A22);

    // Strategy 1: Pivot on A11
    auto A11_inv_opt = strassenInverse(A11);
    if (A11_inv_opt) {
        Matrix A11_inv = *A11_inv_opt;
        Matrix A21_A11_inv = strassenMultiply(A21, A11_inv);
        Matrix S = subtract(A22, strassenMultiply(A21_A11_inv, A12));

        auto S_inv_opt = strassenInverse(S);
        if (!S_inv_opt) return std::nullopt; // Singular
        Matrix S_inv = *S_inv_opt;

        Matrix A11_inv_A12 = strassenMultiply(A11_inv, A12);
        Matrix C12 = strassenMultiply(A11_inv_A12, S_inv);
        Matrix C21 = strassenMultiply(S_inv, A21_A11_inv);
        Matrix C11 = add(A11_inv, strassenMultiply(C12, A21_A11_inv));
        Matrix C22 = S_inv;

        return combine(C11, C12, C21, C22);
    }

    // Strategy 2: Pivot on A22
    auto A22_inv_opt = strassenInverse(A22);
    if (A22_inv_opt) {
        Matrix A22_inv = *A22_inv_opt;
        Matrix A12_A22_inv = strassenMultiply(A12, A22_inv);
        Matrix S_prime = subtract(A11, strassenMultiply(A12_A22_inv, A21));

        auto S_prime_inv_opt = strassenInverse(S_prime);
        if (!S_prime_inv_opt) return std::nullopt; // Singular
        Matrix S_prime_inv = *S_prime_inv_opt;

        Matrix A22_inv_A21 = strassenMultiply(A22_inv, A21);
        Matrix C11 = S_prime_inv;
        Matrix C12 = strassenMultiply(S_prime_inv, A12_A22_inv);
        Matrix C21 = strassenMultiply(A22_inv_A21, S_prime_inv);
        Matrix C22 = add(A22_inv, strassenMultiply(A22_inv_A21, C12));

        return combine(C11, C12, C21, C22);
    }

    return std::nullopt;
}

// --- Top-Level "Driver" Function with Padding ---

// Calculates the next power of two for a given number n
size_t nextPowerOfTwo(size_t n) {
    if (n == 0) return 1;
    if ((n & (n - 1)) == 0) return n; // Already a power of two
    size_t p = 1;
    while (p < n) {
        p <<= 1;
    }
    return p;
}

std::optional<Matrix> invert(const Matrix& A) {
    size_t original_size = A.size();
    if (original_size == 0) {
        return Matrix{};
    }

    size_t padded_size = nextPowerOfTwo(original_size);

    if (original_size == padded_size) {
        // No padding needed, invert directly
        return strassenInverse(A);
    }

    // Create the padded matrix
    Matrix padded_A(padded_size, MatrixRow(padded_size, false));

    // Copy original matrix A to the top-left corner
    for (size_t i = 0; i < original_size; ++i) {
        for (size_t j = 0; j < original_size; ++j) {
            padded_A[i][j] = A[i][j];
        }
    }
    // Place identity matrix in the bottom-right corner
    for (size_t i = original_size; i < padded_size; ++i) {
        padded_A[i][i] = true;
    }

    std::cout << "\nOriginal " << original_size << "x" << original_size 
              << " matrix was padded to " << padded_size << "x" << padded_size << ".\n";
              
    // Invert the padded matrix
    auto padded_inv_opt = strassenInverse(padded_A);

    if (!padded_inv_opt) {
        return std::nullopt; // The padded matrix was singular (implies original was too)
    }

    // Extract the top-left sub-matrix which is the answer
    Matrix result(original_size, MatrixRow(original_size));
    for (size_t i = 0; i < original_size; ++i) {
        for (size_t j = 0; j < original_size; ++j) {
            result[i][j] = (*padded_inv_opt)[i][j];
        }
    }

    return result;
}

// --- Main Program ---

int main() {
    std::ifstream inFile("matrix.txt");
    if (!inFile.is_open()) {
        std::cerr << "Error opening input file: matrix.txt" << std::endl;
        std::cerr << "Creating a sample 5x5 matrix file for you..." << std::endl;
        std::ofstream sampleFile("matrix.txt");
        sampleFile << "1 0 1 1 0\n";
        sampleFile << "0 1 1 0 1\n";
        sampleFile << "1 1 0 1 0\n";
        sampleFile << "1 0 1 0 1\n";
        sampleFile << "0 1 0 1 1\n";
        sampleFile.close();
        inFile.open("matrix.txt");
    }

    Matrix A;
    std::string line;
    while (std::getline(inFile, line)) {
        std::stringstream ss(line);
        int val;
        MatrixRow row;
        while (ss >> val) {
            row.push_back(val != 0);
        }
        if (!row.empty()) {
            A.push_back(row);
        }
    }
    inFile.close();

    if (A.empty() || A.size() != A[0].size()) {
        std::cerr << "Error: Matrix is not square or is empty." << std::endl;
        return 1;
    }

    std::cout << "Reading a " << A.size() << "x" << A.size() << " matrix from matrix.txt..." << std::endl;

    std::cout << "\nOriginal Matrix:\n";
    printMatrix(A);

    auto A_inv_opt = invert(A); // Use the new driver function

    std::ofstream outFile("answer.txt");
    if (!outFile.is_open()) {
        std::cerr << "Error opening output file: answer.txt" << std::endl;
        return 1;
    }

    if (A_inv_opt) {
        Matrix A_inv = *A_inv_opt;
        std::cout << "\nComputed Inverse:\n";
        printMatrix(A_inv);
        printMatrix(A_inv, outFile);
        std::cout << "\nSuccessfully wrote inverse to answer.txt." << std::endl;
        
        // Verification (only works if size is a power of 2, due to strassenMultiply)
        if (A.size() == nextPowerOfTwo(A.size())) {
            std::cout << "\nVerifying... (Original * Inverse) should be Identity Matrix:\n";
            Matrix I = strassenMultiply(A, A_inv);
            printMatrix(I);
        } else {
             std::cout << "\nVerification via Strassen multiply skipped (size is not power of 2)." << std::endl;
        }

    } else {
        std::cerr << "\nThe matrix is singular, or the pivoting strategy failed." << std::endl;
        outFile << "Matrix is singular or could not be inverted with this algorithm.\n";
    }

    outFile.close();
    return 0;
}