#ifndef MATRIX_REPAIR_HPP // CHANGED
#define MATRIX_REPAIR_HPP // CHANGED

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

#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <optional>
#include <bitset>
#include <numeric>   // For std::iota
#include <utility>   // For std::swap and std::pair
#include <algorithm> // For std::find
#include <cstdlib>
#include <ctime>

namespace MatrixRepair {

// --- Type Definitions ---
using IntMatrix = std::vector<std::vector<int>>;

/// Define a maximum matrix size for the bitset-based repair algorithm.
const int MAXN = 2048;
using BitsetMatrix = std::vector<std::bitset<MAXN>>;

/**
 * @brief A structure to hold the results of a matrix repair operation.
 */
struct RepairResult {
    BitsetMatrix repaired_matrix;                ///< The final, repaired (or original) matrix.
    int initial_rank;                          ///< The rank of the matrix before repair.
    int deficiency;                            ///< The rank deficiency (n - rank).
    std::vector<std::pair<int, int>> flips;      ///< A list of (row, col) coordinates that were flipped.
    bool was_already_invertible;               ///< True if no repair was needed.
};

// --- Internal Helper Functions (Anonymous Namespace) ---
namespace {

    // --- Helpers for Singular Matrix Generation ---
    IntMatrix generateRandomMatrix(int n) {
        IntMatrix A(n, std::vector<int>(n, 0));
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                A[i][j] = std::rand() % 2;
            }
        }
        return A;
    }

    int computeRank(const IntMatrix &A) {
        int n = A.size();
        if (n == 0) return 0;
        IntMatrix M = A;
        int rank = 0;
        for (int col = 0; col < n && rank < n; ++col) {
            int pivotRow = rank;
            while (pivotRow < n && M[pivotRow][col] == 0) {
                pivotRow++;
            }
            if (pivotRow < n) {
                std::swap(M[pivotRow], M[rank]);
                for (int row = 0; row < n; ++row) {
                    if (row != rank && M[row][col] == 1) {
                        for (int k = col; k < n; k++) {
                            M[row][k] ^= M[rank][k];
                        }
                    }
                }
                rank++;
            }
        }
        return rank;
    }

    bool isInvertible(const IntMatrix &A) {
        if (A.empty()) return false;
        return computeRank(A) == (int)A.size();
    }
} // end anonymous namespace

// --- Public API Functions ---

/**
 * @brief Generates a random, singular square matrix over GF(2).
 * @param n The dimension of the matrix.
 * @return An optional containing the matrix. Returns std::nullopt if n is invalid.
 */
std::optional<IntMatrix> generateSingularMatrix(int n) {
    if (n <= 0) return std::nullopt;
    std::srand(static_cast<unsigned int>(std::time(0)));

    IntMatrix A;
    // Generate random matrices until a non-invertible one is found.
    do {
        A = generateRandomMatrix(n);
    } while (isInvertible(A));
    
    return A;
}

/**
 * @brief Repairs a singular binary matrix to make it invertible.
 * @param singular_matrix The input matrix to repair.
 * @return A RepairResult struct containing the repaired matrix and details of the operation.
 */
RepairResult repairMatrix(const BitsetMatrix& singular_matrix) {
    const int N = singular_matrix.size();
    if (N == 0) {
        return { {}, 0, 0, {}, false };
    }

    auto mat = singular_matrix;
    std::vector<int> original_row_indices(N);
    std::iota(original_row_indices.begin(), original_row_indices.end(), 0);

    std::vector<int> pivot_cols;
    int rank = 0;
    for (int j = 0; j < N && rank < N; ++j) {
        int p = rank;
        while (p < N && !mat[p].test(j)) ++p;
        if (p < N) {
            std::swap(mat[p], mat[rank]);
            std::swap(original_row_indices[p], original_row_indices[rank]);
            pivot_cols.push_back(j);
            for (int i = 0; i < N; ++i) {
                if (i != rank && mat[i].test(j)) mat[i] ^= mat[rank];
            }
            rank++;
        }
    }

    int deficiency = N - rank;
    if (deficiency == 0) {
        return { singular_matrix, rank, 0, {}, true };
    }

    // Identify dependent rows (in their original positions) and free columns
    std::vector<int> dependent_rows;
    for (int i = rank; i < N; ++i) {
        dependent_rows.push_back(original_row_indices[i]);
    }

    std::vector<int> free_cols;
    for (int j = 0; j < N; ++j) {
        if (std::find(pivot_cols.begin(), pivot_cols.end(), j) == pivot_cols.end()) {
            free_cols.push_back(j);
        }
    }

    // Apply flips to a fresh copy of the original matrix
    auto repaired_mat = singular_matrix;
    std::vector<std::pair<int, int>> flips;
    for (int k = 0; k < deficiency; ++k) {
        int row_to_fix = dependent_rows[k];
        int col_to_flip = free_cols[k];
        repaired_mat[row_to_fix].flip(col_to_flip);
        flips.emplace_back(row_to_fix, col_to_flip);
    }

    return { repaired_mat, rank, deficiency, flips, false };
}


// --- File I/O and Conversion Utilities ---

BitsetMatrix convertToBitsetMatrix(const IntMatrix& int_matrix) {
    int n = int_matrix.size();
    if (n == 0 || n > MAXN) return {};
    BitsetMatrix bs_matrix(n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (int_matrix[i][j] == 1) {
                bs_matrix[i].set(j);
            }
        }
    }
    return bs_matrix;
}

bool saveMatrixToFile(const IntMatrix& matrix, const std::string& filename) {
    std::ofstream fout(filename);
    if (!fout) return false;
    for (const auto& row : matrix) {
        for (size_t j = 0; j < row.size(); ++j) {
            fout << row[j] << (j == row.size() - 1 ? "" : " ");
        }
        fout << "\n";
    }
    return true;
}

bool saveMatrixToFile(const BitsetMatrix& matrix, int n, const std::string& filename) {
    std::ofstream fout(filename);
    if (!fout) return false;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            fout << matrix[i].test(j) << (j == n - 1 ? "" : " ");
        }
        fout << "\n";
    }
    return true;
}

} // namespace MatrixRepair

#endif // MATRIX_REPAIR_HPP // CHANGED