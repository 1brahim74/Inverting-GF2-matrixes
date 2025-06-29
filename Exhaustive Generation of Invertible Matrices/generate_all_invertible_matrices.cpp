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
 * @file generate_all_invertible_matrices.cpp
 * @author Ibrahim Mammadov
 * @contact ibrahim.22@intl.zju.edu.cn
 * @brief Exhaustively generates all unique invertible n x n matrices over GF(2).
 * @version 1.0
 * This program leverages the LUP decomposition theorem, which states that any
 * invertible matrix A can be expressed as A = P * L * U, where:
 *   - P is a permutation matrix.
 *   - L is a unit lower-triangular matrix (1s on the diagonal).
 *   - U is a unit upper-triangular matrix (1s on the diagonal).
 *
 * The program systematically generates every possible P, L, and U matrix for a
 * given size n, computes their product, and stores each unique resulting matrix.
 * The final set of unique matrices is written to a file.
 *
 * @warning This is an exhaustive generator. The number of matrices grows
 *          extremely rapidly with n. It is only feasible for very small n (e.g., n <= 4).
 */

#include <iostream>
#include <vector>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <set>
#include <numeric> // For std::iota

/**
 * @brief Defines a matrix as a 2D vector of integers.
 */
using Matrix = std::vector<std::vector<int>>;

/**
 * @brief Creates an n x m matrix initialized with zeros.
 * @param n The number of rows.
 * @param m The number of columns.
 * @return An n x m matrix filled with zeros.
 */
Matrix createMatrix(int n, int m) {
    return Matrix(n, std::vector<int>(m, 0));
}

/**
 * @brief Multiplies two matrices over the finite field GF(2).
 * @param A The left-hand side matrix.
 * @param B The right-hand side matrix.
 * @return The resulting matrix C = (A * B) mod 2.
 */
Matrix multiplyMatrix(const Matrix& A, const Matrix& B) {
    int n = A.size();
    int m = B[0].size();
    int p = A[0].size();
    Matrix C = createMatrix(n, m);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            for (int k = 0; k < p; ++k) {
                // In GF(2): multiplication is AND, addition is XOR.
                C[i][j] ^= (A[i][k] & B[k][j]);
            }
        }
    }
    return C;
}

/**
 * @brief Generates a permutation matrix from a permutation vector.
 * @param perm A vector where perm[i] = j means row i maps to column j.
 * @return The corresponding n x n permutation matrix.
 */
Matrix generatePermutationMatrix(const std::vector<int>& perm) {
    int n = perm.size();
    Matrix P = createMatrix(n, n);
    for (int i = 0; i < n; ++i) {
        P[i][perm[i]] = 1;
    }
    return P;
}

/**
 * @brief Generates all possible n x n unit lower-triangular matrices.
 * @details A unit lower-triangular matrix has 1s on its diagonal. There are
 *          n*(n-1)/2 entries below the diagonal, each can be 0 or 1. This
 *          function iterates through all 2^(n*(n-1)/2) combinations.
 * @param n The dimension of the matrix.
 * @return A vector containing all possible unit lower-triangular matrices.
 */
std::vector<Matrix> generateAllUnitLower(int n) {
    int bits = n * (n - 1) / 2;
    long long num_matrices = 1LL << bits;
    std::vector<Matrix> result;
    result.reserve(num_matrices);

    for (long long mask = 0; mask < num_matrices; ++mask) {
        Matrix L = createMatrix(n, n);
        int b_idx = 0;
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < i; ++j) {
                L[i][j] = (mask >> b_idx++) & 1;
            }
            L[i][i] = 1;
        }
        result.push_back(L);
    }
    return result;
}

/**
 * @brief Generates all possible n x n unit upper-triangular matrices.
 * @param n The dimension of the matrix.
 * @return A vector containing all possible unit upper-triangular matrices.
 */
std::vector<Matrix> generateAllUnitUpper(int n) {
    int bits = n * (n - 1) / 2;
    long long num_matrices = 1LL << bits;
    std::vector<Matrix> result;
    result.reserve(num_matrices);

    for (long long mask = 0; mask < num_matrices; ++mask) {
        Matrix U = createMatrix(n, n);
        int b_idx = 0;
        for (int i = 0; i < n; ++i) {
            U[i][i] = 1;
            for (int j = i + 1; j < n; ++j) {
                U[i][j] = (mask >> b_idx++) & 1;
            }
        }
        result.push_back(U);
    }
    return result;
}

/**
 * @brief Writes a matrix to the specified output file stream.
 * @param A The matrix to write.
 * @param out The output file stream.
 */
void writeMatrixToFile(const Matrix& A, std::ofstream& out) {
    for (size_t i = 0; i < A.size(); ++i) {
        for (size_t j = 0; j < A[i].size(); ++j) {
            out << A[i][j] << (j == A[i].size() - 1 ? "" : " ");
        }
        out << "\n";
    }
    out << "\n"; // Add a blank line between matrices for readability
}

/**
 * @brief Main function to drive the matrix generation.
 */
int main() {
    int n;
    std::cout << "Enter matrix size n (recommended <= 4): ";
    std::cin >> n;

    if (!std::cin || n <= 0) {
        std::cerr << "Error: Size must be a positive integer." << std::endl;
        return 1;
    }
    if (n > 4) {
        std::cout << "Warning: n > 4 will take a very long time and consume significant memory.\nProceeding...\n";
    }

    // --- 1. Generate all component matrices (P, L, U) ---
    std::cout << "Generating component matrices..." << std::endl;
    auto lowers = generateAllUnitLower(n);
    auto uppers = generateAllUnitUpper(n);

    std::vector<int> p_vec(n);
    std::iota(p_vec.begin(), p_vec.end(), 0); // Initialize p_vec to {0, 1, 2, ...}

    // This set will store unique matrices by flattening them into a single vector.
    std::set<std::vector<int>> seen_matrices;
    int count = 0;

    std::ofstream outFile("unique_GF2_matrices.txt");
    if (!outFile) {
        std::cerr << "Error opening output file." << std::endl;
        return 1;
    }

    // --- 2. Iterate through all P, L, U combinations ---
    std::cout << "Generating and checking all P*L*U combinations..." << std::endl;
    do {
        Matrix P = generatePermutationMatrix(p_vec);
        for (const auto& L : lowers) {
            Matrix PL = multiplyMatrix(P, L); // Pre-calculate P*L
            for (const auto& U : uppers) {
                Matrix A = multiplyMatrix(PL, U);

                // Flatten the matrix to a single vector to use as a key in the set.
                std::vector<int> flat_matrix;
                flat_matrix.reserve(n * n);
                for (const auto& row : A) {
                    flat_matrix.insert(flat_matrix.end(), row.begin(), row.end());
                }

                // If insert succeeds (.second is true), it's a new unique matrix.
                if (seen_matrices.insert(flat_matrix).second) {
                    writeMatrixToFile(A, outFile);
                    ++count;
                }
            }
        }
    } while (std::next_permutation(p_vec.begin(), p_vec.end()));

    outFile.close();

    // --- 3. Report results ---
    // Calculate the theoretical number of invertible matrices over GF(2)
    unsigned long long theoretical_count = 1;
    for (int i = 0; i < n; ++i) {
        unsigned long long term = (1ULL << n) - (1ULL << i);
        theoretical_count *= term;
    }

    std::cout << "\nDone!" << std::endl;
    std::cout << "Found " << count << " unique invertible " << n << "x" << n << " matrices over GF(2)." << std::endl;
    std::cout << "Theoretical count: " << theoretical_count << std::endl;
    std::cout << "Results saved to unique_GF2_matrices.txt" << std::endl;

    return 0;
}