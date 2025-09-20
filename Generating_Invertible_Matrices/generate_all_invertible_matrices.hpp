#ifndef GENERATE_ALL_INVERTIBLE_MATRICES_HPP
#define GENERATE_ALL_INVERTIBLE_MATRICES_HPP

//
// Copyright (c) 2023, Ibrahim Mammadov
// MIT License (same as before)
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
        int n = (int)A.size();
        int m = (int)B[0].size();
        int p = (int)A[0].size();
        Matrix C = createMatrix(n, m);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < m; ++j) {
                int acc = 0;
                for (int k = 0; k < p; ++k) {
                    acc ^= (A[i][k] & B[k][j]); // XOR for GF(2) addition
                }
                C[i][j] = acc;
            }
        }
        return C;
    }

    Matrix generatePermutationMatrix(const std::vector<int>& perm) {
        int n = (int)perm.size();
        Matrix P = createMatrix(n, n);
        for (int i = 0; i < n; ++i) {
            P[i][perm[i]] = 1;
        }
        return P;
    }

    // Old helpers (kept) ------------------------------------------------------

    std::vector<Matrix> generateAllUnitLower(int n) {
        int bits = n * (n - 1) / 2;
        // Safety check: 2^20 is already over a million matrices.
        if (bits > 20) return {};
        
        long long num_matrices = 1LL << bits;
        std::vector<Matrix> result;
        result.reserve((size_t)num_matrices);

        for (long long mask = 0; mask < num_matrices; ++mask) {
            Matrix L = createMatrix(n, n);
            int b_idx = 0;
            for (int i = 0; i < n; ++i) {
                for (int j = 0; j < i; ++j) {
                    L[i][j] = (int)((mask >> b_idx++) & 1LL);
                }
                L[i][i] = 1;
            }
            result.push_back(std::move(L));
        }
        return result;
    }

    std::vector<Matrix> generateAllUnitUpper(int n) {
        int bits = n * (n - 1) / 2;
        if (bits > 20) return {};

        long long num_matrices = 1LL << bits;
        std::vector<Matrix> result;
        result.reserve((size_t)num_matrices);

        for (long long mask = 0; mask < num_matrices; ++mask) {
            Matrix U = createMatrix(n, n);
            int b_idx = 0;
            for (int i = 0; i < n; ++i) {
                U[i][i] = 1;
                for (int j = i + 1; j < n; ++j) {
                    U[i][j] = (int)((mask >> b_idx++) & 1LL);
                }
            }
            result.push_back(std::move(U));
        }
        return result;
    }

    // New helpers for Bruhat-aware generation --------------------------------

    // Build the list of subdiagonal positions (i,j), i>j, that are FREE for L
    // under the permutation p, i.e., those that are NOT inversions of p.
    // Constraint: for any inversion (i>j with p[i] < p[j]) => L_{i,j} MUST be 0.
    // So free positions are pairs (i>j) with p[i] > p[j].
    inline std::vector<std::pair<int,int>>
    freeLowerPositionsForL(const std::vector<int>& p) {
        const int n = (int)p.size();
        std::vector<std::pair<int,int>> freePos;
        freePos.reserve(n*(n-1)/2);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < i; ++j) {
                if (p[i] > p[j]) { // non-inversion relative to p
                    freePos.emplace_back(i, j);
                }
            }
        }
        return freePos;
    }

    // Generate ALL unit lower-triangular L respecting the zero-pattern given by p.
    // Only positions in freePos can be nonzero; all others below diagonal are fixed 0.
    std::vector<Matrix> generateAllUnitLowerBruhat(int n, const std::vector<int>& p) {
        auto freePos = freeLowerPositionsForL(p);
        const int f = (int)freePos.size();
        if (f > 20) {
            // Safety cap like before; you can lift this if you want to go bigger.
            return {};
        }
        const long long total = 1LL << f;
        std::vector<Matrix> result;
        result.reserve((size_t)total);

        for (long long mask = 0; mask < total; ++mask) {
            Matrix L = createMatrix(n, n);
            // Unit diagonal:
            for (int i = 0; i < n; ++i) L[i][i] = 1;
            // Fill ONLY free positions from mask:
            for (int b = 0; b < f; ++b) {
                if ((mask >> b) & 1LL) {
                    auto [i,j] = freePos[b];
                    L[i][j] = 1;
                }
            }
            // All other subdiagonal entries (inversions) remain 0 by construction
            result.push_back(std::move(L));
        }
        return result;
    }

} // end anonymous namespace

// --- Public API Functions ---

/**
 * @brief OLD brute-force enumeration via A = P * L * U. Kept for reference.
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
                    unique_results.push_back(std::move(A));
                }
            }
        }
    } while (std::next_permutation(p_vec.begin(), p_vec.end()));

    return unique_results;
}

/**
 * @brief NEW Bruhat-aware enumeration: A = P * L * U with L constrained by permutation inversions.
 *
 * For permutation p, subdiagonal entries (i>j) with p[i] < p[j] are FORCED zeros in L.
 * Only (i>j) with p[i] > p[j] are free. U is any unit upper-triangular matrix.
 *
 * @param n The dimension (practical n <= 6–7 with current safety caps).
 * @return A vector of unique invertible matrices; empty if n too large for caps.
 */
std::vector<Matrix> generateAllBruhat(int n) {
    if (n <= 0) return {};

    // Pre-enumerate all unit upper-triangular U (unconstrained)
    auto uppers = generateAllUnitUpper(n);
    if (uppers.empty()) return {};

    // Permutations
    std::vector<int> p_vec(n);
    std::iota(p_vec.begin(), p_vec.end(), 0);

    std::set<std::vector<int>> seen;       // Safety: ensure uniqueness
    std::vector<Matrix> results;
    do {
        // For this permutation, enumerate only L that respect the Bruhat zero pattern
        auto lowersForP = generateAllUnitLowerBruhat(n, p_vec);
        if (lowersForP.empty()) {
            // If cap hit (too many free bits), abort to avoid explosion.
            return {};
        }

        Matrix P = generatePermutationMatrix(p_vec);
        for (const auto& L : lowersForP) {
            Matrix PL = multiplyMatrix(P, L);
            for (const auto& U : uppers) {
                Matrix A = multiplyMatrix(PL, U);

                // Flatten to deduplicate:
                std::vector<int> flat;
                flat.reserve(n*n);
                for (const auto& r : A) flat.insert(flat.end(), r.begin(), r.end());

                if (seen.insert(flat).second) {
                    results.push_back(std::move(A));
                }
            }
        }
    } while (std::next_permutation(p_vec.begin(), p_vec.end()));

    return results;
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
                outFile << row[j] << (j + 1 == row.size() ? "" : " ");
            }
            outFile << "\n";
        }
        outFile << "\n"; // Blank line separator
    }
    return true;
}

} // namespace ExhaustiveGenerator

#endif // GENERATE_ALL_INVERTIBLE_MATRICES_HPP
