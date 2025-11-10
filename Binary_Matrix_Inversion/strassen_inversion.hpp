#ifndef STRASSEN_INVERSION_HPP
#define STRASSEN_INVERSION_HPP

#include <vector>
#include <optional>
#include <string>
#include <iostream>
#include <numeric>
#include <stdexcept>

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

        struct LUPResult {
            BoolMatrix P, L, U;
        };

        std::optional<LUPResult> fastLUP(const BoolMatrix& A);
        std::optional<BoolMatrix> invertLowerTriangular(const BoolMatrix& L);
        std::optional<BoolMatrix> invertUpperTriangular(const BoolMatrix& U);

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

        // CORRECTED base case: Gaussian elimination with pivoting to find PA=LU
        std::optional<LUPResult> gaussJordanLUP(const BoolMatrix& A_in) {
            size_t n = A_in.size();
            if (n == 0) return LUPResult{{}, {}, {}};

            BoolMatrix U = A_in;
            BoolMatrix L = identity(n);
            std::vector<size_t> p_vec(n);
            std::iota(p_vec.begin(), p_vec.end(), 0);

            for (size_t j = 0; j < n; ++j) {
                size_t pivot_row = j;
                while (pivot_row < n && !U[pivot_row][j]) {
                    pivot_row++;
                }
                if (pivot_row == n) return std::nullopt; // Singular

                std::swap(U[j], U[pivot_row]);
                std::swap(p_vec[j], p_vec[pivot_row]);
                // Swap the computed part of L matrix as well
                for(size_t col = 0; col < j; ++col) {
                    bool temp = L[j][col];
                    L[j][col] = L[pivot_row][col];
                    L[pivot_row][col] = temp;
                }

                for (size_t i = j + 1; i < n; ++i) {
                    if (U[i][j]) {
                        L[i][j] = true;
                        for (size_t k = j; k < n; ++k) {
                            U[i][k] = U[i][k] ^ U[j][k];
                        }
                    }
                }
            }
            
            BoolMatrix P(n, std::vector<bool>(n, false));
            for(size_t i=0; i<n; ++i) P[i][p_vec[i]] = true;

            return LUPResult{P, L, U};
        }


        std::optional<BoolMatrix> invertLowerTriangular(const BoolMatrix& L) {
            size_t n = L.size(); if (n == 0) return BoolMatrix{};
            for (size_t i = 0; i < n; ++i) if (!L[i][i]) return std::nullopt;
            if (n % 2 != 0 || n <= STRASSEN_CUTOFF) {
                BoolMatrix L_inv = identity(n);
                for (size_t j = 0; j < n; ++j) for (size_t i = j + 1; i < n; ++i) {
                    bool sum = false; for (size_t k = j; k < i; ++k) sum ^= (L[i][k] & L_inv[k][j]);
                    L_inv[i][j] = sum;
                }
                return L_inv;
            }
            size_t k = n/2;
            BoolMatrix L11(k,std::vector<bool>(k)), L12(k,std::vector<bool>(k)), L21(k,std::vector<bool>(k)), L22(k,std::vector<bool>(k));
            split(L, L11, L12, L21, L22);
            auto L11_inv_opt = invertLowerTriangular(L11); auto L22_inv_opt = invertLowerTriangular(L22);
            BoolMatrix C21 = strassenMultiply(strassenMultiply(*L22_inv_opt, L21), *L11_inv_opt);
            return combine(*L11_inv_opt, BoolMatrix(k, std::vector<bool>(k, false)), C21, *L22_inv_opt);
        }

        std::optional<BoolMatrix> invertUpperTriangular(const BoolMatrix& U) {
            size_t n = U.size(); if (n == 0) return BoolMatrix{};
            for (size_t i = 0; i < n; ++i) if (!U[i][i]) return std::nullopt;
            if (n % 2 != 0 || n <= STRASSEN_CUTOFF) {
                BoolMatrix U_inv = identity(n);
                for (size_t j = n-1; j != (size_t)-1; --j) for (size_t i = j-1; i != (size_t)-1; --i) {
                    bool sum = false; for (size_t k = i+1; k <= j; ++k) sum ^= (U[i][k] & U_inv[k][j]);
                    U_inv[i][j] = sum;
                }
                return U_inv;
            }
            size_t k = n/2;
            BoolMatrix U11(k,std::vector<bool>(k)), U12(k,std::vector<bool>(k)), U21(k,std::vector<bool>(k)), U22(k,std::vector<bool>(k));
            split(U, U11, U12, U21, U22);
            auto U11_inv_opt = invertUpperTriangular(U11); auto U22_inv_opt = invertUpperTriangular(U22);
            BoolMatrix C12 = strassenMultiply(strassenMultiply(*U11_inv_opt, U12), *U22_inv_opt);
            return combine(*U11_inv_opt, C12, BoolMatrix(k, std::vector<bool>(k, false)), *U22_inv_opt);
        }

        std::optional<LUPResult> fastLUP(const BoolMatrix& A) {
            size_t n = A.size();
            if (n % 2 != 0 || n <= STRASSEN_CUTOFF) return gaussJordanLUP(A);
            size_t k = n/2;
            BoolMatrix A11(k,std::vector<bool>(k)), A12(k,std::vector<bool>(k)), A21(k,std::vector<bool>(k)), A22(k,std::vector<bool>(k));
            split(A, A11, A12, A21, A22);

            auto lup11 = fastLUP(A11); if (!lup11) return std::nullopt;
            auto& [P1, L1, U1] = *lup11;
            auto L1_inv_opt = invertLowerTriangular(L1); auto U1_inv_opt = invertUpperTriangular(U1);
            
            auto P1_A12 = strassenMultiply(P1, A12);
            auto U12 = strassenMultiply(*L1_inv_opt, P1_A12);
            auto L21 = strassenMultiply(A21, *U1_inv_opt);
            auto S = add(A22, strassenMultiply(L21, U12));
            
            auto lup_s = fastLUP(S); if (!lup_s) return std::nullopt;
            auto& [P2, L2, U2] = *lup_s;

            auto P1_full = combine(P1, BoolMatrix(k,std::vector<bool>(k, false)), BoolMatrix(k,std::vector<bool>(k, false)), identity(k));
            auto P2_full = combine(identity(k), BoolMatrix(k,std::vector<bool>(k, false)), BoolMatrix(k,std::vector<bool>(k, false)), P2);

            auto P_final = strassenMultiply(P2_full, P1_full);
            
            auto P2_L21 = strassenMultiply(P2, L21);
            auto L_final = combine(L1, BoolMatrix(k, std::vector<bool>(k, false)), P2_L21, L2);
            
            auto U_final = combine(U1, U12, BoolMatrix(k, std::vector<bool>(k, false)), U2);
            
            // P needs to be transposed to get the final correct row permutation matrix from the p_vec representation
            BoolMatrix P_final_transposed(n, std::vector<bool>(n));
            for(size_t i=0; i<n; ++i) for(size_t j=0; j<n; ++j) P_final_transposed[j][i] = P_final[i][j];

            return LUPResult{P_final_transposed, L_final, U_final};
        }
    } // end anonymous namespace

    // --- Public Function Implementations ---
    inline BoolMatrix strassenMultiply(const BoolMatrix& A, const BoolMatrix& B) {
        size_t n = A.size();
        if (n == 0) return {};
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
        auto lup_opt = fastLUP(A); if (!lup_opt) return std::nullopt;
        auto& [P, L, U] = *lup_opt;
        auto L_inv_opt = invertLowerTriangular(L);
        auto U_inv_opt = invertUpperTriangular(U);
        if (!L_inv_opt || !U_inv_opt) return std::nullopt;
        auto temp = strassenMultiply(*U_inv_opt, *L_inv_opt);
        return strassenMultiply(temp, P);
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
        for (size_t i = 0; i < n; ++i) {
            if (A[i].size() != n) return false;
            for (size_t j = 0; j < n; ++j) {
                if (A[i][j] != (i == j)) return false;
            }
        }
        return true;
    }

} // namespace MatrixOps

#endif // STRASSEN_INVERSION_HPP