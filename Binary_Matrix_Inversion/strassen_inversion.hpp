#ifndef STRASSEN_INVERSION_HPP
#define STRASSEN_INVERSION_HPP

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

    // --- Forward Declarations ---
    inline BoolMatrix strassenMultiply(const BoolMatrix& A, const BoolMatrix& B);
    inline std::optional<BoolMatrix> invertMatrixRobust(const BoolMatrix& A);

    // --- Helper Function Implementations ---
    inline std::optional<BoolMatrix> solveLower(const BoolMatrix& L, const BoolMatrix& B) {
        size_t n = L.size();
        if (n == 0) return BoolMatrix{};
        if (B.empty()) return BoolMatrix(n, std::vector<bool>(0));
        
        BoolMatrix X = B;
        for (size_t j = 0; j < n; ++j) {
             for (size_t i = j + 1; i < n; ++i) {
                if (L[i][j]) {
                    for(size_t k=0; k < X[0].size(); ++k) X[i][k] = X[i][k] ^ X[j][k];
                }
            }
        }
        return X;
    }
    
    inline std::optional<BoolMatrix> solveUpper(const BoolMatrix& U, const BoolMatrix& B) {
        size_t n = U.size();
        if (n == 0) return BoolMatrix{};
        if (B.empty()) return BoolMatrix(n, std::vector<bool>(0));

        BoolMatrix X = B;
        for (int j = n - 1; j >= 0; --j) {
            for (int i = 0; i < j; ++i) {
                if (U[i][(size_t)j]) {
                     for(size_t k=0; k < X[0].size(); ++k) X[i][k] = X[i][k] ^ X[(size_t)j][k];
                }
            }
        }
        return X;
    }

    // --- Internal Implementation Details ---
    namespace {

        const size_t LUP_CUTOFF = 32;

        struct LUPResult {
            BoolMatrix P, L, U;
        };
        
        BoolMatrix identity(size_t n) {
            BoolMatrix I(n, std::vector<bool>(n, false));
            for (size_t i = 0; i < n; ++i) I[i][i] = true;
            return I;
        }

        BoolMatrix add(const BoolMatrix& A, const BoolMatrix& B) {
            size_t n = A.size();
            size_t m = n > 0 ? A[0].size() : 0;
            BoolMatrix C(n, std::vector<bool>(m));
            for (size_t i = 0; i < n; ++i) for (size_t j = 0; j < m; ++j) C[i][j] = A[i][j] ^ B[i][j];
            return C;
        }
        
        BoolMatrix transpose(const BoolMatrix& A) {
            if (A.empty()) return {};
            size_t r = A.size();
            size_t c = A[0].size();
            BoolMatrix A_T(c, std::vector<bool>(r));
            for(size_t i = 0; i < r; ++i) for(size_t j = 0; j < c; ++j) A_T[j][i] = A[i][j];
            return A_T;
        }

        BoolMatrix combine(const BoolMatrix& C11, const BoolMatrix& C12, const BoolMatrix& C21, const BoolMatrix& C22) {
            size_t r1 = C11.empty() ? (C12.empty() ? 0 : C12.size()) : C11.size();
            size_t c1 = C11.empty() ? (C21.empty() ? 0 : C21[0].size()) : C11[0].size();
            size_t r2 = C21.empty() ? (C22.empty() ? 0 : C22.size()) : C21.size();
            size_t c2 = C12.empty() ? (C22.empty() ? 0 : C22[0].size()) : C12[0].size();

            size_t n_rows = r1 + r2;
            size_t n_cols = c1 + c2;
            BoolMatrix C(n_rows, std::vector<bool>(n_cols));

            if(r1 > 0 && c1 > 0) for(size_t i = 0; i < r1; ++i) for(size_t j = 0; j < c1; ++j) C[i][j] = C11[i][j];
            if(r1 > 0 && c2 > 0) for(size_t i = 0; i < r1; ++i) for(size_t j = 0; j < c2; ++j) C[i][j + c1] = C12[i][j];
            if(r2 > 0 && c1 > 0) for(size_t i = 0; i < r2; ++i) for(size_t j = 0; j < c1; ++j) C[i + r1][j] = C21[i][j];
            if(r2 > 0 && c2 > 0) for(size_t i = 0; i < r2; ++i) for(size_t j = 0; j < c2; ++j) C[i + r1][j + c1] = C22[i][j];
            
            return C;
        }

        std::optional<LUPResult> gaussJordanLUP(const BoolMatrix& A_in) {
            size_t n = A_in.size();
            if (n == 0) return LUPResult{identity(0), identity(0), identity(0)};
            BoolMatrix U = A_in;
            BoolMatrix L = identity(n);
            std::vector<size_t> p_vec(n);
            std::iota(p_vec.begin(), p_vec.end(), 0);
            for (size_t j = 0; j < n; ++j) {
                size_t pivot_row = j;
                while (pivot_row < n && !U[pivot_row][j]) pivot_row++;
                if (pivot_row == n) return std::nullopt;
                if (pivot_row != j) {
                    std::swap(U[j], U[pivot_row]);
                    std::swap(p_vec[j], p_vec[pivot_row]);
                    for(size_t col = 0; col < j; ++col) {
                        bool temp = L[j][col]; L[j][col] = L[pivot_row][col]; L[pivot_row][col] = temp;
                    }
                }
                for (size_t i = j + 1; i < n; ++i) {
                    if (U[i][j]) {
                        L[i][j] = true;
                        for (size_t k = j; k < n; ++k) U[i][k] = U[i][k] ^ U[j][k];
                    }
                }
            }
            BoolMatrix P(n, std::vector<bool>(n, false));
            for(size_t i = 0; i < n; ++i) P[i][p_vec[i]] = true;
            return LUPResult{P, L, U};
        }
        
        std::optional<LUPResult> fastLUP(const BoolMatrix& A) {
            size_t n = A.size();
            if (n <= LUP_CUTOFF) return gaussJordanLUP(A);

            size_t k = n / 2;
            size_t m = n - k;

            BoolMatrix A11(k, std::vector<bool>(k)), A12(k, std::vector<bool>(m));
            BoolMatrix A21(m, std::vector<bool>(k)), A22(m, std::vector<bool>(m));
            for(size_t r=0; r<k; ++r) for(size_t c=0; c<k; ++c) A11[r][c] = A[r][c];
            for(size_t r=0; r<k; ++r) for(size_t c=0; c<m; ++c) A12[r][c] = A[r][c+k];
            for(size_t r=0; r<m; ++r) for(size_t c=0; c<k; ++c) A21[r][c] = A[r+k][c];
            for(size_t r=0; r<m; ++r) for(size_t c=0; c<m; ++c) A22[r][c] = A[r+k][c+k];

            BoolMatrix A_left = combine(A11, {}, A21, {});
            auto lup_left = gaussJordanLUP(A_left);
            if (!lup_left) return std::nullopt;
            
            auto& [P_left, L_left_full, U11] = *lup_left;
            
            BoolMatrix L11(k, std::vector<bool>(k));
            BoolMatrix L21(m, std::vector<bool>(k));
            for(size_t r=0; r<k; ++r) for(size_t c=0; c<k; ++c) L11[r][c] = L_left_full[r][c];
            for(size_t r=0; r<m; ++r) for(size_t c=0; c<k; ++c) L21[r][c] = L_left_full[r+k][c];
            
            BoolMatrix A_right_panel = combine(A12, {}, A22, {});
            BoolMatrix A_right_permuted = strassenMultiply(P_left, A_right_panel);
            
            BoolMatrix A12_p(k, std::vector<bool>(m));
            BoolMatrix A22_p(m, std::vector<bool>(m));
            for(size_t r=0; r<k; ++r) for(size_t c=0; c<m; ++c) A12_p[r][c] = A_right_permuted[r][c];
            for(size_t r=0; r<m; ++r) for(size_t c=0; c<m; ++c) A22_p[r][c] = A_right_permuted[r+k][c];

            auto U12_opt = solveLower(L11, A12_p);
            if (!U12_opt) return std::nullopt;
            auto& U12 = *U12_opt;

            auto S = add(A22_p, strassenMultiply(L21, U12));
            auto lupS = fastLUP(S);
            if (!lupS) return std::nullopt;
            auto& [P_S, L22, U22] = *lupS;

            BoolMatrix P_bottom_ext = identity(n);
            for(size_t r=0; r<m; ++r) for(size_t c=0; c<m; ++c) P_bottom_ext[r+k][c+k] = P_S[r][c];
            BoolMatrix P_final = strassenMultiply(P_bottom_ext, P_left);
            
            // **** THE FINAL BUG FIX IS HERE ****
            // This is the correct way to construct L_final.
            BoolMatrix L_final = identity(n);
            BoolMatrix L21_permuted = strassenMultiply(P_S, L21);
            for(size_t r=0; r<k; ++r) for(size_t c=0; c<k; ++c) L_final[r][c] = L11[r][c];
            for(size_t r=0; r<m; ++r) for(size_t c=0; c<k; ++c) L_final[r+k][c] = L21_permuted[r][c];
            for(size_t r=0; r<m; ++r) for(size_t c=0; c<m; ++c) L_final[r+k][c+k] = L22[r][c];
            
            BoolMatrix U_final = combine(U11, U12, BoolMatrix(m, std::vector<bool>(k, false)), U22);

            return LUPResult{P_final, L_final, U_final};
        }
    } // End anonymous namespace

    inline std::optional<BoolMatrix> invertMatrixRobust(const BoolMatrix& A) {
        if (A.empty()) return {};
        if (A.size() != A[0].size()) throw std::runtime_error("Matrix must be square.");
        
        auto lup_opt = fastLUP(A);
        if (!lup_opt) return std::nullopt;
        
        auto& [P, L, U] = *lup_opt;
        
        auto U_inv_opt = solveUpper(U, identity(A.size()));
        if (!U_inv_opt) return std::nullopt;
        
        auto L_inv_U_inv_opt = solveLower(L, *U_inv_opt);
        if (!L_inv_U_inv_opt) return std::nullopt;

        return strassenMultiply(*L_inv_U_inv_opt, P);
    }
    
    inline void printMatrix(const BoolMatrix& A, std::ostream& out, const std::string& title = "") {
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

    inline BoolMatrix strassenMultiply(const BoolMatrix& A, const BoolMatrix& B) {
        size_t n = A.size();
        if (n == 0) return {};
        if (A.size() != A[0].size() || B.size() != B[0].size() || A.size() != B.size() || n % 2 != 0 || n <= LUP_CUTOFF) {
            BoolMatrix C(n, std::vector<bool>(n, false));
            for (size_t i = 0; i < n; ++i) {
                for (size_t j = 0; j < n; ++j) {
                    for (size_t l = 0; l < n; ++l) {
                        C[i][j] = C[i][j] ^ (A[i][l] & B[l][j]);
                    }
                }
            }
            return C;
        }
        
        size_t k = n / 2;
        BoolMatrix A11(k,std::vector<bool>(k)), A12(k,std::vector<bool>(k)), A21(k,std::vector<bool>(k)), A22(k,std::vector<bool>(k));
        BoolMatrix B11(k,std::vector<bool>(k)), B12(k,std::vector<bool>(k)), B21(k,std::vector<bool>(k)), B22(k,std::vector<bool>(k));
        for(size_t i=0; i<k; ++i) for(size_t j=0; j<k; ++j) { A11[i][j]=A[i][j]; A12[i][j]=A[i][j+k]; A21[i][j]=A[i+k][j]; A22[i][j]=A[i+k][j+k]; }
        for(size_t i=0; i<k; ++i) for(size_t j=0; j<k; ++j) { B11[i][j]=B[i][j]; B12[i][j]=B[i][j+k]; B21[i][j]=B[i+k][j]; B22[i][j]=B[i+k][j+k]; }
        
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
        
        BoolMatrix C(n, std::vector<bool>(n));
        for(size_t i=0; i<k; ++i) for(size_t j=0; j<k; ++j) { C[i][j]=C11[i][j]; C[i][j+k]=C12[i][j]; C[i+k][j]=C21[i][j]; C[i+k][j+k]=C22[i][j]; }
        return C;
    }

} // namespace MatrixOps

#endif // STRASSEN_INVERSION_HPP