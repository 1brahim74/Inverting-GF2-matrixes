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
#include <set>
#include <numeric>   // For std::iota
#include <algorithm> // For std::next_permutation

/**
 * @brief A namespace for the exhaustive invertible matrix generator.
 */
namespace ExhaustiveGenerator {

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

    Matrix generatePermutationMatrix(const std::vector<int>& perm) {
        int n = perm.size();
        Matrix P = createMatrix(n, n);
        for (int i = 0; i < n; ++i) {
            P[i][perm[i]] = 1;
        }
        return P;
    }

    std::vector<Matrix> generateAllUnitLower(int n) {
        int bits = n * (n - 1) / 2;
        // Safety check: 2^20 is already over a million matrices.
        if (bits > 20) return {};
        
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

    std::vector<Matrix> generateAllUnitUpper(int n) {
        int bits = n * (n - 1) / 2;
        if (bits > 20) return {};

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

} // end anonymous namespace

// --- Public API Functions ---

/**
 * @brief Exhaustively generates all unique invertible n x n matrices over GF(2).
 * @param n The dimension of the matrices. Only feasible for small n (<= 4).
 * @return A vector of matrices. Returns an empty vector if n is invalid or too large.
 */
std::vector<Matrix> generateAll(int n) {
    if (n <= 0) return {};

    auto lowers = generateAllUnitLower(n);
    auto uppers = generateAllUnitUpper(n);
    if (lowers.empty() || uppers.empty()) {
        // This check catches cases where n is too large for the helpers.
        return {};
    }

    std::vector<int> p_vec(n);
    std::iota(p_vec.begin(), p_vec.end(), 0);

    std::set<std::vector<int>> seen_matrices;
    std::vector<Matrix> unique_results;

    do {
        Matrix P = generatePermutationMatrix(p_vec);
        for (const auto& L : lowers) {
            Matrix PL = multiplyMatrix(P, L);
            for (const auto& U : uppers) {
                Matrix A = multiplyMatrix(PL, U);

                std::vector<int> flat_matrix;
                flat_matrix.reserve(n * n);
                for (const auto& row : A) {
                    flat_matrix.insert(flat_matrix.end(), row.begin(), row.end());
                }

                if (seen_matrices.insert(flat_matrix).second) {
                    unique_results.push_back(A);
                }
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

} // namespace ExhaustiveGenerator

#endif // GENERATE_ALL_INVERTIBLE_MATRICES_HPP