#ifndef GENERATE_ALL_INVERTIBLE_MATRICES_HPP
#define GENERATE_ALL_INVERTIBLE_MATRICES_HPP

#include <vector>
#include <string>
#include <numeric>   // For std::iota
#include <algorithm> // For std::next_permutation
#include <utility>   // For std::pair
#include <iterator>  // For iterator traits
#include <omp.h>     // For OpenMP

/**
 * @brief A namespace for generating invertible matrices using Bruhat decomposition.
 */
namespace InvertibleMatrixGenerator {

using Matrix = std::vector<std::vector<int>>;

// --- Internal Helper Functions (Anonymous Namespace) ---
namespace {
    // These functions are implementation details and not part of the public API.

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

    Matrix generatePermutationMatrix(int n, const std::vector<int>& perm) {
        Matrix P = createMatrix(n, n);
        for (int i = 0; i < n; ++i) {
            P[i][perm[i]] = 1;
        }
        return P;
    }

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

    Matrix generateConstrainedUnitLower(int n, const std::vector<std::pair<int, int>>& allowed_pos, long long mask) {
        Matrix L = createMatrix(n, n);
        for (int i = 0; i < n; ++i) L[i][i] = 1;
        for (size_t k = 0; k < allowed_pos.size(); ++k) {
            if ((mask >> k) & 1) {
                L[allowed_pos[k].first][allowed_pos[k].second] = 1;
            }
        }
        return L;
    }

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


/**
 * @class BruhatGenerator
 * @brief A C++17 range-based generator for iterating through all invertible n x n matrices over GF(2).
 *
 * This class allows iterating through the entire General Linear Group GL(n, 2) without
 * storing all matrices in memory. It is highly memory-efficient.
 *
 * Usage:
 *   BruhatGenerator generator(n);
 *   for (const auto& matrix : generator) {
 *       // process matrix
 *   }
 */
class BruhatGenerator {
public:
    // --- The Iterator Class ---
    class iterator {
    public:
        // C++ Iterator Traits
        using iterator_category = std::input_iterator_tag;
        using value_type = Matrix;
        using difference_type = std::ptrdiff_t;
        using pointer = const Matrix*;
        using reference = const Matrix&;

        // Dereference operator
        reference operator*() const { return current_matrix_; }
        pointer operator->() const { return &current_matrix_; }

        // Pre-increment operator (advances the state)
        iterator& operator++() {
            advance();
            return *this;
        }

        // Equality comparison
        bool operator!=(const iterator& other) const {
            return is_done_ != other.is_done_;
        }
        
        // Default constructor for the "end" iterator
        iterator() : is_done_(true) {}

        // Constructor for the "begin" iterator
        explicit iterator(int n) : n_(n), is_done_(n <= 0) {
            if (!is_done_) {
                p_vec_.resize(n);
                std::iota(p_vec_.begin(), p_vec_.end(), 0);

                num_U_matrices_ = 1LL << (n * (n - 1) / 2);
                u_mask_ = 0;
                l_mask_ = 0;

                update_permutation_dependent_state();
                compute_current_matrix();
            }
        }

    private:
        void update_permutation_dependent_state() {
            inversion_pos_ = getInversionPositionsForL(p_vec_);
            num_L_matrices_ = 1LL << inversion_pos_.size();
            PT_ = transposeMatrix(generatePermutationMatrix(n_, p_vec_));
        }

        void compute_current_matrix() {
            Matrix L = generateConstrainedUnitLower(n_, inversion_pos_, l_mask_);
            Matrix U = generateUnitUpperFromMask(n_, u_mask_);
            current_matrix_ = multiplyMatrix(multiplyMatrix(PT_, L), U);
        }

        void advance() {
            u_mask_++;
            if (u_mask_ < num_U_matrices_) {
                // Common case: just advance U
            } else {
                u_mask_ = 0;
                l_mask_++;
                if (l_mask_ < num_L_matrices_) {
                    // Advance L
                } else {
                    l_mask_ = 0;
                    if (std::next_permutation(p_vec_.begin(), p_vec_.end())) {
                        // Advance to next permutation
                        update_permutation_dependent_state();
                    } else {
                        // All permutations are done
                        is_done_ = true;
                        return;
                    }
                }
            }
            compute_current_matrix();
        }

        int n_ = 0;
        bool is_done_ = false;

        // State variables
        std::vector<int> p_vec_;
        std::vector<std::pair<int, int>> inversion_pos_;
        long long u_mask_ = 0, l_mask_ = 0;
        long long num_U_matrices_ = 0, num_L_matrices_ = 0;
        Matrix PT_; // Cached transposed permutation matrix

        Matrix current_matrix_;
    };

    // --- BruhatGenerator public methods ---
    explicit BruhatGenerator(int n) : n_(n) {}

    iterator begin() const { return iterator(n_); }
    iterator end() const { return iterator(); }

private:
    int n_;
};


/**
 * @brief Generates all unique invertible n x n matrices in parallel and processes them with a callback.
 * @param n The dimension of the matrices.
 * @param callback A function (e.g., a lambda) that will be called for each generated matrix.
 *                 The callback must be thread-safe if it modifies shared data.
 *
 * This function is the fastest way to process all matrices for n >= 4, as it uses
 * OpenMP to distribute the work across all available CPU cores.
 *
 * Example:
 *   std::atomic<long long> count = 0;
 *   generateAllParallel(n, [&](const Matrix& m) {
 *       count++;
 *   });
 */
template<typename Func>
void generateAllParallel(int n, Func process_matrix) {
    if (n <= 0) return;

    // 1. Pre-generate all permutations
    std::vector<std::vector<int>> all_perms;
    std::vector<int> p_vec(n);
    std::iota(p_vec.begin(), p_vec.end(), 0);
    do {
        all_perms.push_back(p_vec);
    } while (std::next_permutation(p_vec.begin(), p_vec.end()));

    const long long num_U_matrices = 1LL << (n * (n - 1) / 2);

    // 2. Parallel loop over the permutations (this is the main work distribution)
    #pragma omp parallel for schedule(dynamic)
    for (size_t i = 0; i < all_perms.size(); ++i) {
        const auto& p = all_perms[i];

        // Each thread computes its own P^T and L structure
        Matrix PT = transposeMatrix(generatePermutationMatrix(n, p));
        auto inversion_pos = getInversionPositionsForL(p);
        long long num_L_matrices = 1LL << inversion_pos.size();

        // Inner loops are executed by each thread for its assigned permutation
        for (long long l_mask = 0; l_mask < num_L_matrices; ++l_mask) {
            Matrix L = generateConstrainedUnitLower(n, inversion_pos, l_mask);
            Matrix PTL = multiplyMatrix(PT, L); // Pre-compute for the inner-most loop

            for (long long u_mask = 0; u_mask < num_U_matrices; ++u_mask) {
                Matrix U = generateUnitUpperFromMask(n, u_mask);
                Matrix A = multiplyMatrix(PTL, U);
                
                // Call the user-provided function to process the matrix
                process_matrix(A);
            }
        }
    }
}


} // namespace InvertibleMatrixGenerator

#endif // GENERATE_ALL_INVERTIBLE_MATRICES_HPP