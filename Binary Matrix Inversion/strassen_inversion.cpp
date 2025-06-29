//
// Copyright (c) 2023, Ibrahim Mammadov
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
 * @file strassen_inversion.cpp
 * @author Ibrahim Mammadov
 * @contact ibrahim.22@intl.zju.edu.cn
 * @brief Implements a pivoted, recursive block matrix inversion algorithm over
 *        the finite field GF(2), using Strassen's algorithm for multiplication.
 * @version 0.9.0
 * This program reads a square binary matrix of a power-of-2 dimension from
 * a file named "matrix.txt". It then computes the inverse of this matrix
 * using a recursive method that can pivot between the A11 and A22 blocks
 * to handle cases where the top-left sub-matrix is singular. The sub-problem
 * multiplications are performed using Strassen's algorithm for efficiency.
 * The result is written to "answer.txt".
 */

#include <iostream>
#include <vector>
#include <bitset>
#include <optional>
#include <numeric>
#include <fstream>
#include <iomanip>
#include <cmath>

// --- Type Definitions and Constants ---

template <size_t N>
using Matrix = std::vector<std::bitset<N>>;

const size_t MATRIX_SIZE = 4;
const size_t STRASSEN_CUTOFF = 2;

// --- Forward Declarations ---

template <size_t N>
std::optional<Matrix<N>> strassenInverse(const Matrix<N>& A);

template <size_t N>
Matrix<N> strassenMultiply(const Matrix<N>& A, const Matrix<N>& B);

// --- Utility and Helper Functions ---

template <size_t N>
void printMatrix(const Matrix<N>& A, std::ostream& out = std::cout) {
    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < N; ++j) {
            out << A[i][j] << " ";
        }
        out << "\n";
    }
}

template <size_t N>
Matrix<N> add(const Matrix<N>& A, const Matrix<N>& B) {
    Matrix<N> C(N);
    for (size_t i = 0; i < N; ++i) {
        C[i] = A[i] ^ B[i];
    }
    return C;
}

template <size_t N>
Matrix<N> subtract(const Matrix<N>& A, const Matrix<N>& B) {
    return add(A, B);
}

template <size_t N>
void split(const Matrix<N>& A, Matrix<N / 2>& A11, Matrix<N / 2>& A12, Matrix<N / 2>& A21, Matrix<N / 2>& A22) {
    const size_t k = N / 2;
    for (size_t i = 0; i < k; ++i) {
        for (size_t j = 0; j < k; ++j) {
            A11[i][j] = A[i][j];
            A12[i][j] = A[i][j + k];
            A21[i][j] = A[i + k][j];
            A22[i][j] = A[i + k][j + k];
        }
    }
}

template <size_t N>
Matrix<N> combine(const Matrix<N / 2>& A11, const Matrix<N / 2>& A12, const Matrix<N / 2>& A21, const Matrix<N / 2>& A22) {
    Matrix<N> A(N);
    const size_t k = N / 2;
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

template <size_t N>
std::optional<Matrix<N>> gaussJordanInverse(Matrix<N> A) {
    Matrix<N> I(N);
    for (size_t i = 0; i < N; ++i) I[i][i] = 1;

    for (size_t j = 0; j < N; ++j) {
        size_t pivot_row = j;
        while (pivot_row < N && !A[pivot_row][j]) {
            pivot_row++;
        }
        if (pivot_row == N) return std::nullopt;
        std::swap(A[j], A[pivot_row]);
        std::swap(I[j], I[pivot_row]);
        for (size_t i = 0; i < N; ++i) {
            if (i != j && A[i][j]) {
                A[i] ^= A[j];
                I[i] ^= I[j];
            }
        }
    }
    return I;
}

template <size_t N>
Matrix<N> strassenMultiply(const Matrix<N>& A, const Matrix<N>& B) {
    if (N <= STRASSEN_CUTOFF) {
        Matrix<N> C(N);
        Matrix<N> B_T(N);
        for (size_t i = 0; i < N; ++i) {
            for (size_t j = 0; j < N; ++j) {
                if (B[i][j]) B_T[j][i] = 1;
            }
        }
        for (size_t i = 0; i < N; i++) {
            for (size_t j = 0; j < N; j++) {
                if ((A[i] & B_T[j]).count() % 2 != 0) {
                    C[i][j] = 1;
                }
            }
        }
        return C;
    }

    const size_t k = N / 2;
    Matrix<k> A11(k), A12(k), A21(k), A22(k);
    Matrix<k> B11(k), B12(k), B21(k), B22(k);

    split(A, A11, A12, A21, A22);
    split(B, B11, B12, B21, B22);

    Matrix<k> P1 = strassenMultiply(add(A11, A22), add(B11, B22));
    Matrix<k> P2 = strassenMultiply(add(A21, A22), B11);
    Matrix<k> P3 = strassenMultiply(A11, subtract(B12, B22));
    Matrix<k> P4 = strassenMultiply(A22, subtract(B21, B11));
    Matrix<k> P5 = strassenMultiply(add(A11, A12), B22);
    Matrix<k> P6 = strassenMultiply(subtract(A21, A11), add(B11, B12));
    Matrix<k> P7 = strassenMultiply(subtract(A12, A22), add(B21, B22));

    Matrix<k> C11 = add(subtract(add(P1, P4), P5), P7);
    Matrix<k> C12 = add(P3, P5);
    Matrix<k> C21 = add(P2, P4);
    Matrix<k> C22 = add(subtract(add(P1, P3), P2), P6);

    return combine<N>(C11, C12, C21, C22);
}

/**
 * @brief Computes the inverse of a matrix using a pivoted recursive blockwise method.
 *
 * @details This function first attempts to use the standard block inversion
 *          formula, which requires the top-left sub-matrix (A11) to be invertible.
 *          If A11 is singular, it then attempts an alternative block inversion
 *          formula that requires the bottom-right sub-matrix (A22) to be
 *          invertible. This pivoting strategy makes the algorithm significantly
 *          more robust.
 *
 * @warning This improved algorithm is more robust but may still fail if both
 *          diagonal blocks (A11 and A22) are singular, even if the full matrix
 *          is invertible. A perfect implementation would require more complex
 *          pivoting strategies (e.g., LUP decomposition).
 *
 * @tparam N The dimension of the matrix.
 * @param A The matrix to invert.
 * @return An std::optional containing the inverse, or std::nullopt if the
 *         matrix is singular or the pivoting strategy fails.
 */
template <size_t N>
std::optional<Matrix<N>> strassenInverse(const Matrix<N>& A) {
    if (N <= STRASSEN_CUTOFF) {
        return gaussJordanInverse(A);
    }

    const size_t k = N / 2;
    Matrix<k> A11(k), A12(k), A21(k), A22(k);
    split(A, A11, A12, A21, A22);

    // --- Strategy 1: Try to pivot on A11 (standard case) ---
    auto A11_inv_opt = strassenInverse(A11);
    if (A11_inv_opt) {
        Matrix<k> A11_inv = *A11_inv_opt;
        Matrix<k> A21_A11_inv = strassenMultiply(A21, A11_inv);
        Matrix<k> S = subtract(A22, strassenMultiply(A21_A11_inv, A12));

        auto S_inv_opt = strassenInverse(S);
        if (!S_inv_opt) return std::nullopt; // Singular matrix
        Matrix<k> S_inv = *S_inv_opt;

        Matrix<k> A11_inv_A12 = strassenMultiply(A11_inv, A12);
        Matrix<k> C12 = strassenMultiply(A11_inv_A12, S_inv);
        Matrix<k> C21 = strassenMultiply(S_inv, A21_A11_inv);
        Matrix<k> C11 = add(A11_inv, strassenMultiply(C12, A21_A11_inv));
        Matrix<k> C22 = S_inv;

        return combine<N>(C11, C12, C21, C22);
    }

    // --- Strategy 2: A11 was singular, try to pivot on A22 ---
    auto A22_inv_opt = strassenInverse(A22);
    if (A22_inv_opt) {
        Matrix<k> A22_inv = *A22_inv_opt;
        Matrix<k> A12_A22_inv = strassenMultiply(A12, A22_inv);
        Matrix<k> S_prime = subtract(A11, strassenMultiply(A12_A22_inv, A21));

        auto S_prime_inv_opt = strassenInverse(S_prime);
        if (!S_prime_inv_opt) return std::nullopt; // Singular matrix
        Matrix<k> S_prime_inv = *S_prime_inv_opt;

        Matrix<k> A22_inv_A21 = strassenMultiply(A22_inv, A21);
        Matrix<k> C11 = S_prime_inv;
        Matrix<k> C12 = strassenMultiply(S_prime_inv, A12_A22_inv);
        Matrix<k> C21 = strassenMultiply(A22_inv_A21, S_prime_inv);
        Matrix<k> C22 = add(A22_inv, strassenMultiply(A22_inv_A21, C12));

        return combine<N>(C11, C12, C21, C22);
    }

    // --- Failure: Both A11 and A22 are singular, this algorithm cannot proceed.
    return std::nullopt;
}

int main() {
    static_assert((MATRIX_SIZE > 0) && ((MATRIX_SIZE & (MATRIX_SIZE - 1)) == 0),
                  "MATRIX_SIZE must be a power of 2 for this Strassen implementation.");

    std::ifstream inFile("matrix.txt");
    if (!inFile.is_open()) {
        std::cerr << "Error opening input file: matrix.txt" << std::endl;
        std::cerr << "Creating a sample invertible 4x4 matrix for you..." << std::endl;
        std::ofstream sampleFile("matrix.txt");
        sampleFile << "1 0 1 1\n";
        sampleFile << "0 1 1 0\n";
        sampleFile << "1 1 0 1\n";
        sampleFile << "1 0 1 0\n";
        sampleFile.close();
        inFile.open("matrix.txt");
    }

    Matrix<MATRIX_SIZE> A(MATRIX_SIZE);
    std::cout << "Reading a " << MATRIX_SIZE << "x" << MATRIX_SIZE << " matrix from matrix.txt..." << std::endl;
    for (size_t i = 0; i < MATRIX_SIZE; ++i) {
        for (size_t j = 0; j < MATRIX_SIZE; ++j) {
            int val;
            if (!(inFile >> val)) {
                std::cerr << "Error reading matrix from file. Ensure it is "
                          << MATRIX_SIZE << "x" << MATRIX_SIZE << "." << std::endl;
                return 1;
            }
            if (val % 2 != 0) {
                A[i][j] = 1;
            }
        }
    }
    inFile.close();

    std::cout << "\nOriginal Matrix:\n";
    printMatrix(A);

    auto A_inv_opt = strassenInverse(A);

    std::ofstream outFile("answer.txt");
    if (!outFile.is_open()) {
        std::cerr << "Error opening output file: answer.txt" << std::endl;
        return 1;
    }

    if (A_inv_opt) {
        Matrix<MATRIX_SIZE> A_inv = *A_inv_opt;
        std::cout << "\nComputed Inverse:\n";
        printMatrix(A_inv);
        printMatrix(A_inv, outFile);
        std::cout << "\nSuccessfully wrote inverse to answer.txt." << std::endl;
        std::cout << "\nVerifying... (Original * Inverse) should be Identity Matrix:\n";
        Matrix<MATRIX_SIZE> I = strassenMultiply(A, A_inv);
        printMatrix(I);
    } else {
        std::cerr << "\nThe matrix is singular, or the pivoting strategy failed." << std::endl;
        outFile << "Matrix is singular or could not be inverted with this algorithm.\n";
    }

    outFile.close();
    return 0;
}