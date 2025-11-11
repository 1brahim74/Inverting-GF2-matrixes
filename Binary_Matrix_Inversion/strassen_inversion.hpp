#ifndef STRASSEN_INVERSION_HPP
#define STRASSEN_INVERSION_HPP

#include <vector>
#include <optional>
#include <string>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <algorithm> // For std::swap

namespace MatrixOps {

    // --- Public Type Definitions ---
    using BoolMatrix = std::vector<std::vector<bool>>;

    // --- Public Function Declarations ---
    inline std::optional<BoolMatrix> invertMatrixRobust(const BoolMatrix& A);
    inline BoolMatrix strassenMultiply(const BoolMatrix& A, const BoolMatrix& B);
    inline void printMatrix(const BoolMatrix& A, std::ostream& out, const std::string& title = "");

    // --- Internal Implementation Details ---
    namespace {

        // CUTOFF: Below this size, Strassen's has too much overhead.
        const size_t STRASSEN_CUTOFF = 16; 

        struct LUPResult {
            BoolMatrix P, L, U;
        };
        
        // --- Utility Functions ---
        
        BoolMatrix identity(size_t n) {
            BoolMatrix I(n, std::vector<bool>(n, false));
            for (size_t i = 0; i < n; ++i) I[i][i] = true;
            return I;
        }

        BoolMatrix add(const BoolMatrix& A, const BoolMatrix& B) {
            size_t n = A.size();
            BoolMatrix C(n, std::vector<bool>(n));
            for (size_t i = 0; i < n; ++i) {
                for (size_t j = 0; j < n; ++j) {
                    C[i][j] = A[i][j] ^ B[i][j];
                }
            }
            return C;
        }
        
        BoolMatrix transpose(const BoolMatrix& A) {
            if (A.empty()) return {};
            size_t n = A.size();
            BoolMatrix A_T(n, std::vector<bool>(n));
            for(size_t i = 0; i < n; ++i) {
                for(size_t j = 0; j < n; ++j) {
                    A_T[j][i] = A[i][j];
                }
            }
            return A_T;
        }

        // --- THE ITERATIVE BUNCH-HOPCROFT ALGORITHM ---

        // This is a direct, iterative implementation of LUP factorization with pivoting.
        // It forms the basis of a robust solution.
        std::optional<LUPResult> iterativeLUP(const BoolMatrix& A) {
            size_t n = A.size();
            if (n == 0) return LUPResult{{}, {}, {}};

            BoolMatrix LU = A;
            std::vector<size_t> p_vec(n);
            std::iota(p_vec.begin(), p_vec.end(), 0);

            for (size_t i = 0; i < n; ++i) {
                // --- Pivoting Step ---
                size_t pivot_row = i;
                while (pivot_row < n && !LU[pivot_row][i]) {
                    pivot_row++;
                }
                if (pivot_row == n) return std::nullopt; // Singular matrix

                if (pivot_row != i) {
                    std::swap(LU[i], LU[pivot_row]);
                    std::swap(p_vec[i], p_vec[pivot_row]);
                }

                // --- Elimination Step ---
                for (size_t j = i + 1; j < n; ++j) {
                    if (LU[j][i]) { // If the element is 1
                        LU[j][i] = true; // Store the multiplier Lji in place
                        // Update the rest of the row (j)
                        for (size_t k = i + 1; k < n; ++k) {
                            LU[j][k] = LU[j][k] ^ LU[i][k];
                        }
                    }
                }
            }

            // --- Deconstruct LU into separate L and U matrices ---
            BoolMatrix L = identity(n);
            BoolMatrix U = BoolMatrix(n, std::vector<bool>(n, false));
            for (size_t i = 0; i < n; i++) {
                for (size_t j = 0; j < n; j++) {
                    if (i > j) {
                        L[i][j] = LU[i][j];
                    } else {
                        U[i][j] = LU[i][j];
                    }
                }
            }
            
            // --- Construct Permutation Matrix P from p_vec ---
            BoolMatrix P(n, std::vector<bool>(n, false));
            for(size_t i = 0; i < n; ++i) P[i][p_vec[i]] = true;

            return LUPResult{P, L, U};
        }


        // --- Inversion Helpers ---
        
        std::optional<BoolMatrix> invertLowerTriangular(const BoolMatrix& L) {
            size_t n = L.size(); if (n == 0) return BoolMatrix{};
            BoolMatrix L_inv = identity(n);
            for (size_t j = 0; j < n; ++j) {
                if (!L[j][j]) return std::nullopt;
                for (size_t i = j + 1; i < n; ++i) {
                    bool sum = false; for (size_t k = j; k < i; ++k) sum ^= (L[i][k] & L_inv[k][j]);
                    L_inv[i][j] = sum;
                }
            }
            return L_inv;
        }

        std::optional<BoolMatrix> invertUpperTriangular(const BoolMatrix& U) {
            size_t n = U.size(); if (n == 0) return BoolMatrix{};
            BoolMatrix U_inv = identity(n);
            for (int j = n-1; j >= 0; --j) {
                 if (!U[j][j]) return std::nullopt;
                for (int i = j-1; i >= 0; --i) {
                    bool sum = false; for (size_t k = i+1; k <= (size_t)j; ++k) sum ^= (U[i][k] & U_inv[k][j]);
                    U_inv[i][j] = sum;
                }
            }
            return U_inv;
        }

    } // end anonymous namespace

    // --- Public Function Implementations ---
    
    inline BoolMatrix strassenMultiply(const BoolMatrix& A, const BoolMatrix& B) {
        size_t n = A.size();
        if (n == 0 || B.empty() || A[0].size() != B.size()) return {};

        if (n <= STRASSEN_CUTOFF || A.size() != A[0].size() || B.size() != B[0].size() || A.size() != B.size() || n % 2 != 0) { // Standard multiply
            size_t m = B[0].size();
            size_t p = B.size();
            BoolMatrix C(n, std::vector<bool>(m, false));
            for (size_t i = 0; i < n; ++i) {
                for (size_t j = 0; j < m; ++j) {
                    for (size_t l = 0; l < p; ++l) {
                        C[i][j] = C[i][j] ^ (A[i][l] & B[l][j]);
                    }
                }
            }
            return C;
        }

        size_t k = n / 2;
        BoolMatrix A11(k,std::vector<bool>(k)), A12(k,std::vector<bool>(k)), A21(k,std::vector<bool>(k)), A22(k,std::vector<bool>(k));
        BoolMatrix B11(k,std::vector<bool>(k)), B12(k,std::vector<bool>(k)), B21(k,std::vector<bool>(k)), B22(k,std::vector<bool>(k));
        // Simple split for square matrices
        for (size_t i = 0; i < k; ++i) for (size_t j = 0; j < k; ++j) { A11[i][j] = A[i][j]; A12[i][j] = A[i][j + k]; A21[i][j] = A[i + k][j]; A22[i][j] = A[i + k][j + k]; }
        for (size_t i = 0; i < k; ++i) for (size_t j = 0; j < k; ++j) { B11[i][j] = B[i][j]; B12[i][j] = B[i][j + k]; B21[i][j] = B[i + k][j]; B22[i][j] = B[i + k][j + k]; }

        BoolMatrix P1 = strassenMultiply(add(A11, A22), add(B11, B22));
        BoolMatrix P2 = strassenMultiply(add(A21, A22), B11);
        BoolMatrix P3 = strassenMultiply(A11, add(B12, B22));
        BoolMatrix P4 = strassenMultiply(A22, add(B21, B11));
        BoolMatrix P5 = strassenMultiply(add(A11, A12), B22);
        BoolMatrix P6 = strassenMultiply(add(A21, A11), add(B11, B12));
        BoolMatrix P7 = strassenMultiply(add(A12, A22), add(B21, B22));
        BoolMatrix C11 = add(add(add(P1, P4), P7), P5);
        BoolMatrix C12 = add(P3, P5);
        BoolMatrix C21 = add(P2, P4);
        BoolMatrix C22 = add(add(add(P1, P3), P6), P2);
        
        // Combine for square matrices
        BoolMatrix C(n, std::vector<bool>(n));
        for (size_t i = 0; i < k; ++i) for (size_t j = 0; j < k; ++j) { C[i][j] = C11[i][j]; C[i][j + k] = C12[i][j]; C[i + k][j] = C21[i][j]; C[i + k][j + k] = C22[i][j]; }
        return C;
    }
    
    inline std::optional<BoolMatrix> invertMatrixRobust(const BoolMatrix& A) {
        if (A.empty()) return BoolMatrix{};
        if (A.size() != A[0].size()) throw std::runtime_error("Matrix must be square.");
        
        // A^-1 = (P^-1 * L * U)^-1 = U^-1 * L^-1 * P
        // Use the robust iterative LUP decomposition
        auto lup_opt = iterativeLUP(A); 
        if (!lup_opt) return std::nullopt;
        
        auto& [P, L, U] = *lup_opt;
        
        auto L_inv_opt = invertLowerTriangular(L);
        auto U_inv_opt = invertUpperTriangular(U);
        if (!L_inv_opt || !U_inv_opt) return std::nullopt;
        
        // Use Strassen for the final multiplications if matrices are large enough
        auto temp = strassenMultiply(*U_inv_opt, *L_inv_opt);
        return strassenMultiply(temp, transpose(P)); // P^-1 is P^T
    }

    inline void printMatrix(const BoolMatrix& A, std::ostream& out, const std::string& title) {
        if (!title.empty()) { out << title << " (" << A.size() << "x" << (A.empty() ? 0 : A[0].size()) << "):\n"; }
        if (A.empty()) { out << "[Empty Matrix]\n"; return; }
        for (const auto& row : A) {
            for (size_t j = 0; j < row.size(); ++j) {
                out << row[j] << (j == row.size() - 1 ? "" : " ");
            }
            out << "\n";
        }
        out << std::endl;
    }

} // namespace MatrixOps

#endif // STRASSEN_INVERSION_HPP