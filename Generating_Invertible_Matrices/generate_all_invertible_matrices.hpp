#ifndef GENERATE_ALL_INVERTIBLE_MATRICES_HPP
#define GENERATE_ALL_INVERTIBLE_MATRICES_HPP

//
// Copyright (c) 2023, Ibrahim Mammadov
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, to the following conditions:
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

#include <vector>
#include <string>
#include <fstream>
#include <numeric>   // For std::iota
#include <algorithm> // For std::next_permutation
#include <utility>   // For std::pair

/**
 * @brief A namespace for generating invertible matrices using Bruhat decomposition.
 */
namespace InvertibleMatrixGenerator {

using Matrix = std::vector<std::vector<int>>;

// --- Internal Helper Functions (Anonymous Namespace) ---
namespace {

    Matrix createMatrix(int n, int m) {
        return Matrix(n, std::vector<int>(m, 0));
    }

    Matrix multiplyMatrix(const Matrix& A, const Matrix& B) {
        int n = A.size();
        int m = B[0].size();
        int p = A[0].size();
        Matrix C = createMatrix(n, m);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < m; ++j) {
                for (int k = 0; k < p; ++k) {
                    C[i][j] ^= (A[i][k] & B[k][j]); // XOR for GF(2) addition
                }
            }
        }
        return C;
    }
    
    Matrix transposeMatrix(const Matrix& M) {
        int n = M.size();
        Matrix T = createMatrix(n, n);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                T[j][i] = M[i][j];
            }
        }
        return T;
    }

    Matrix generatePermutationMatrix(const std::vector<int>& perm) {
        int n = perm.size();
        Matrix P = createMatrix(n, n);
        for (int i = 0; i < n; ++i) {
            P[i][perm[i]] = 1;
        }
        return P;
    }

    // Finds positions (i, j) where i > j that correspond to inversions in p.
    // An inversion is a pair (row_a, row_b) where row_a < row_b but p[row_a] > p[row_b].
    // These are the positions where L can have a '1'.
    std::vector<std::pair<int, int>> getInversionPositionsForL(const std::vector<int>& perm) {
        int n = perm.size();
        std::vector<std::pair<int, int>> positions;
        for (int i = 1; i < n; ++i) {
            for (int j = 0; j < i; ++j) {
                if (perm[j] > perm[i]) {
                    positions.push_back({i, j});
                }
            }
        }
        return positions;
    }

    // Generates a Unit Lower triangular matrix constrained by allowed '1' positions.
    Matrix generateConstrainedUnitLower(int n, const std::vector<std::pair<int, int>>& allowed_pos, long long mask) {
        Matrix L = createMatrix(n, n);
        for (int i = 0; i < n; ++i) {
            L[i][i] = 1;
        }
        for (size_t k = 0; k < allowed_pos.size(); ++k) {
            if ((mask >> k) & 1) {
                L[allowed_pos[k].first][allowed_pos[k].second] = 1;
            }
        }
        return L;
    }

    // Generates a Unit Upper triangular matrix from a bitmask.
    Matrix generateUnitUpperFromMask(int n, long long mask) {
        Matrix U = createMatrix(n, n);
        int b_idx = 0;
        for (int i = 0; i < n; ++i) {
            U[i][i] = 1;
            for (int j = i + 1; j < n; ++j) {
                U[i][j] = (mask >> b_idx++) & 1;
            }
        }
        return U;
    }

} // end anonymous namespace

// --- Public API Functions ---

/**
 * @brief Generates all unique invertible n x n matrices over GF(2) using Bruhat decomposition.
 * @param n The dimension of the matrices. Feasible for n <= 5.
 * @return A vector of all unique invertible matrices.
 */
std::vector<Matrix> generateAll(int n) {
    if (n <= 0) return {};
    int bits = n * (n - 1) / 2;
    // Safety check for upper matrix generation. 2^15 is already large.
    if (bits > 15 && n > 5) return {};

    std::vector<Matrix> unique_results;
    
    // The theoretical count can be very large, reserve memory if reasonable.
    // For n=4, 20160. For n=5, 9999360.
    if (n == 4) unique_results.reserve(20160);

    std::vector<int> p_vec(n);
    std::iota(p_vec.begin(), p_vec.end(), 0);

    long long num_U_matrices = 1LL << bits;

    // --- Main Loop: Iterate through each Bruhat cell, defined by a permutation P ---
    do {
        Matrix P = generatePermutationMatrix(p_vec);
        Matrix PT = transposeMatrix(P); // P^-1 = P^T

        // 1. Determine the structure of L matrices for this cell
        auto inversion_pos = getInversionPositionsForL(p_vec);
        long long num_L_matrices = 1LL << inversion_pos.size();

        // 2. Generate all matrices in this cell: A = P^T * L * U
        for (long long l_mask = 0; l_mask < num_L_matrices; ++l_mask) {
            Matrix L = generateConstrainedUnitLower(n, inversion_pos, l_mask);
            Matrix PTL = multiplyMatrix(PT, L);

            for (long long u_mask = 0; u_mask < num_U_matrices; ++u_mask) {
                Matrix U = generateUnitUpperFromMask(n, u_mask);
                Matrix A = multiplyMatrix(PTL, U);
                unique_results.push_back(A);
            }
        }
    } while (std::next_permutation(p_vec.begin(), p_vec.end()));

    return unique_results;
}

/**
 * @brief Saves a collection of matrices to a file.
 * @param matrices The vector of matrices to save.
 * @param filename The name of the file to save to.
 * @return true on success, false on file opening error.
 */
bool saveMatricesToFile(const std::vector<Matrix>& matrices, const std::string& filename) {
    std::ofstream outFile(filename);
    if (!outFile.is_open()) {
        return false;
    }
    for (const auto& matrix : matrices) {
        for (const auto& row : matrix) {
            for (size_t j = 0; j < row.size(); ++j) {
                outFile << row[j] << (j == row.size() - 1 ? "" : " ");
            }
            outFile << "\n";
        }
        outFile << "\n"; // Blank line separator
    }
    return true;
}

} // namespace InvertibleMatrixGenerator

#endif // GENERATE_ALL_INVERTIBLE_MATRICES_HPP