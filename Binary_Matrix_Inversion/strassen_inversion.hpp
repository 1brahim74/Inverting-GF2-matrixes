#ifndef ROBUST_INVERSION_HPP
#define ROBUST_INVERSION_HPP

#include <vector>
#include <optional>
#include <string>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <algorithm>

namespace MatrixOps {

    // --- Public Type Definitions ---
    using BoolMatrix = std::vector<std::vector<bool>>;

    // --- Public Function Declarations ---
    inline std::optional<BoolMatrix> invertMatrixRobust(const BoolMatrix& A);
    inline BoolMatrix strassenMultiply(const BoolMatrix& A, const BoolMatrix& B);
    inline void printMatrix(const BoolMatrix& A, std::ostream& out, const std::string& title = "");
    inline bool isIdentity(const BoolMatrix& A);

    // --- Internal Implementation Details ---
    namespace {

        const size_t STRASSEN_CUTOFF = 16;
        
        // Forward declare the main recursive inversion function
        std::optional<BoolMatrix> fastInvert(const BoolMatrix& A);

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

        void split(const BoolMatrix& A, BoolMatrix& A11, BoolMatrix& A12, BoolMatrix& A21, BoolMatrix& A22) {
            size_t k = A.size() / 2;
            for (size_t i = 0; i < k; ++i) {
                for (size_t j = 0; j < k; ++j) {
                    A11[i][j] = A[i][j]; A12[i][j] = A[i][j + k];
                    A21[i][j] = A[i + k][j]; A22[i][j] = A[i + k][j + k];
                }
            }
        }

        BoolMatrix combine(const BoolMatrix& C11, const BoolMatrix& C12, const BoolMatrix& C21, const BoolMatrix& C22) {
            size_t k = C11.size();
            if (k == 0) return {};
            size_t n = 2 * k;
            BoolMatrix C(n, std::vector<bool>(n));
            for (size_t i = 0; i < k; ++i) {
                for (size_t j = 0; j < k; ++j) {
                    C[i][j] = C11[i][j]; C[i][j + k] = C12[i][j];
                    C[i + k][j] = C21[i][j]; C[i + k][j + k] = C22[i][j];
                }
            }
            return C;
        }

        BoolMatrix standardMultiply(const BoolMatrix& A, const BoolMatrix& B) {
            size_t n = A.size();
            if (n == 0) return {};
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

        // Base case for inversion: Gaussian-Jordan elimination (O(n^3))
        std::optional<BoolMatrix> gaussJordanInvert(const BoolMatrix& A) {
            size_t n = A.size();
            if (n == 0) return BoolMatrix{};
            
            BoolMatrix M = A;
            BoolMatrix I = identity(n);

            for (size_t j = 0; j < n; ++j) {
                size_t pivot_row = j;
                while (pivot_row < n && !M[pivot_row][j]) {
                    pivot_row++;
                }
                if (pivot_row == n) return std::nullopt; // Singular

                std::swap(M[j], M[pivot_row]);
                std::swap(I[j], I[pivot_row]);
                
                for (size_t i = 0; i < n; ++i) {
                    if (i != j && M[i][j]) {
                        for (size_t k = 0; k < n; ++k) {
                            M[i][k] = M[i][k] ^ M[j][k];
                            I[i][k] = I[i][k] ^ I[j][k];
                        }
                    }
                }
            }
            return I;
        }
        
        // Fully robust O(n^log2(7)) recursive matrix inversion.
        std::optional<BoolMatrix> fastInvert(const BoolMatrix& A) {
            size_t n = A.size();
            if (n % 2 != 0 || n <= STRASSEN_CUTOFF) {
                return gaussJordanInvert(A);
            }

            size_t k = n / 2;
            BoolMatrix A11(k,std::vector<bool>(k)), A12(k,std::vector<bool>(k)), 
                         A21(k,std::vector<bool>(k)), A22(k,std::vector<bool>(k));
            split(A, A11, A12, A21, A22);

            // --- Path 1: Assume A11 is invertible and use the standard block inversion formula ---
            auto A11_inv_opt = fastInvert(A11);
            if (A11_inv_opt) {
                BoolMatrix A11_inv = *A11_inv_opt;
                BoolMatrix S12 = strassenMultiply(A21, A11_inv);
                BoolMatrix S21 = strassenMultiply(A11_inv, A12);
                BoolMatrix S = add(A22, strassenMultiply(S12, A12)); // S = A22 - A21*A11_inv*A12
                
                auto S_inv_opt = fastInvert(S);
                if (!S_inv_opt) return std::nullopt; // Schur complement is singular, so A is singular
                BoolMatrix S_inv = *S_inv_opt;

                BoolMatrix C12 = strassenMultiply(S21, S_inv);
                BoolMatrix C21 = strassenMultiply(S_inv, S12);
                BoolMatrix C11 = add(A11_inv, strassenMultiply(C12, S12));
                
                // In GF(2), -X = X, so the formula simplifies
                return combine(C11, C12, C21, S_inv);
            }

            // --- Path 2: A11 is singular. Try assuming A22 is invertible and use the alternate formula ---
            auto A22_inv_opt = fastInvert(A22);
            if (A22_inv_opt) {
                BoolMatrix A22_inv = *A22_inv_opt;
                // T is the Schur complement of A22
                BoolMatrix T = add(A11, strassenMultiply(strassenMultiply(A12, A22_inv), A21)); // T = A11 - A12*A22_inv*A21

                auto T_inv_opt = fastInvert(T);
                if(!T_inv_opt) return std::nullopt; // Schur complement is singular, so A is singular
                BoolMatrix T_inv = *T_inv_opt;

                BoolMatrix C12 = strassenMultiply(strassenMultiply(T_inv, A12), A22_inv);
                BoolMatrix C21 = strassenMultiply(strassenMultiply(A22_inv, A21), T_inv);
                BoolMatrix C22 = add(A22_inv, strassenMultiply(C21, A12));

                return combine(T_inv, C12, C21, C22);
            }

            // --- Path 3: Both A11 and A22 are singular. Perform a permutation pivot. ---
            // We form a new matrix A' by swapping the block rows: A' = [[A21, A22], [A11, A12]]
            // If A is invertible, then A21 must be invertible in this case. We pivot on A21.
            auto A21_inv_opt = fastInvert(A21);
            if (A21_inv_opt) {
                BoolMatrix A21_inv = *A21_inv_opt;
                
                // Apply the standard block inversion formula to A', where:
                // A'_11=A21, A'_12=A22, A'_21=A11, A'_22=A12
                BoolMatrix S12_p = strassenMultiply(A11, A21_inv);
                BoolMatrix S21_p = strassenMultiply(A21_inv, A22);
                BoolMatrix S_p   = add(A12, strassenMultiply(S12_p, A22));
                
                auto S_p_inv_opt = fastInvert(S_p);
                if (!S_p_inv_opt) return std::nullopt; // A is singular
                BoolMatrix S_p_inv = *S_p_inv_opt;

                // Calculate blocks of (A')^-1, let's call them B11, B12, B21, B22
                BoolMatrix B12 = strassenMultiply(S21_p, S_p_inv);
                BoolMatrix B21 = strassenMultiply(S_p_inv, S12_p);
                BoolMatrix B11 = add(A21_inv, strassenMultiply(B12, S12_p));
                BoolMatrix B22 = S_p_inv;
                
                // We calculated B = (A')^-1 = (P*A)^-1 = A^-1 * P^-1
                // We need A^-1 = B * P, where P is the block row swap matrix.
                // Multiplying by P on the right swaps the block columns of B.
                // A^-1 = [[B12, B11], [B22, B21]]
                return combine(B12, B11, B22, B21);
            }

            // If all three paths fail, the original matrix A is truly singular.
            return std::nullopt;
        }

    } // end anonymous namespace

    // --- Public Function Implementations ---
    inline BoolMatrix strassenMultiply(const BoolMatrix& A, const BoolMatrix& B) {
        size_t n = A.size();
        if (n == 0) return {};
        if (A.empty() || B.empty() || A[0].size() != B.size()) throw std::runtime_error("Matrix dimensions are incompatible for multiplication.");
        
        if (n != A[0].size() || n != B.size() || n != B[0].size()){
             return standardMultiply(A, B);
        }

        if (n % 2 != 0 || n <= STRASSEN_CUTOFF) return standardMultiply(A, B);
        size_t k = n / 2;
        BoolMatrix A11(k,std::vector<bool>(k)), A12(k,std::vector<bool>(k)), A21(k,std::vector<bool>(k)), A22(k,std::vector<bool>(k));
        BoolMatrix B11(k,std::vector<bool>(k)), B12(k,std::vector<bool>(k)), B21(k,std::vector<bool>(k)), B22(k,std::vector<bool>(k));
        split(A, A11, A12, A21, A22); split(B, B11, B12, B21, B22);
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
        return combine(C11, C12, C21, C22);
    }

    inline std::optional<BoolMatrix> invertMatrixRobust(const BoolMatrix& A) {
        if (A.empty()) return BoolMatrix{};
        if (A.size() != A[0].size()) throw std::runtime_error("Matrix must be square.");
        return fastInvert(A);
    }

    inline void printMatrix(const BoolMatrix& A, std::ostream& out, const std::string& title) {
        if (!title.empty()) { out << title << " (" << A.size() << "x" << (A.empty() ? 0 : A[0].size()) << "):\n"; }
        if (A.empty()) { out << "[Empty Matrix]\n"; return; }
        for (const auto& row : A) {
            for (bool val : row) { out << val << " "; }
            out << "\n";
        }
        out << std::endl;
    }

    inline bool isIdentity(const BoolMatrix& A) {
        if (A.empty()) return true;
        size_t n = A.size();
        if (n == 0) return true;
        for (size_t i = 0; i < n; ++i) {
            if (A[i].size() != n) return false;
            for (size_t j = 0; j < n; ++j) {
                if (A[i][j] != (i == j)) return false;
            }
        }
        return true;
    }

} // namespace MatrixOps

#endif // ROBUST_INVERSION_HPP