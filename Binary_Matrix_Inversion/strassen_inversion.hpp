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

    using BoolMatrix = std::vector<std::vector<bool>>;

    inline std::optional<BoolMatrix> invertMatrixRobust(const BoolMatrix& A);
    inline BoolMatrix strassenMultiply(const BoolMatrix& A, const BoolMatrix& B);
    inline void printMatrix(const BoolMatrix& A, std::ostream& out, const std::string& title = "");
    inline bool isIdentity(const BoolMatrix& A);

    namespace {

        const size_t RECURSION_CUTOFF = 5;
        
        std::optional<BoolMatrix> fastInvertRecursive(const BoolMatrix& A);
        
        BoolMatrix identity(size_t n) {
            BoolMatrix I(n, std::vector<bool>(n, false));
            for (size_t i = 0; i < n; ++i) I[i][i] = true;
            return I;
        }

        BoolMatrix add(const BoolMatrix& A, const BoolMatrix& B) {
            size_t n = A.size();
            if (n == 0) return {};
            size_t m = A[0].size();
            BoolMatrix C(n, std::vector<bool>(m));
            for (size_t i = 0; i < n; ++i) {
                for (size_t j = 0; j < m; ++j) {
                    C[i][j] = A[i][j] ^ B[i][j];
                }
            }
            return C;
        }

        void split(const BoolMatrix& A, BoolMatrix& A11, BoolMatrix& A12, BoolMatrix& A21, BoolMatrix& A22) {
            size_t k = A.size() / 2;
            for (size_t i = 0; i < k; ++i) {
                for (size_t j = 0; j < k; ++j) {
                    A11[i][j] = A[i][j];
                    A12[i][j] = A[i][j + k];
                    A21[i][j] = A[i + k][j];
                    A22[i][j] = A[i + k][j + k];
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
                    C[i][j] = C11[i][j];
                    C[i][j + k] = C12[i][j];
                    C[i + k][j] = C21[i][j];
                    C[i + k][j + k] = C22[i][j];
                }
            }
            return C;
        }

        BoolMatrix standardMultiply(const BoolMatrix& A, const BoolMatrix& B) {
            size_t n = A.size();
            if (n == 0) return {};
            if (B.empty()) return BoolMatrix(n, std::vector<bool>(0));
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
                if (pivot_row == n) return std::nullopt;
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

        std::optional<BoolMatrix> fastInvertRecursive(const BoolMatrix& A) {
            size_t n = A.size();
            if (n == 0) return BoolMatrix{};

            if (n <= RECURSION_CUTOFF) {
                return gaussJordanInvert(A);
            }

            size_t piv_i = n, piv_j = n;
            for (size_t i = 0; i < n && piv_i == n; ++i) {
                for (size_t j = 0; j < n; ++j) {
                    if (A[i][j]) {
                        piv_i = i;
                        piv_j = j;
                        break;
                    }
                }
            }

            if (piv_i == n) {
                return std::nullopt;
            }

            BoolMatrix permutedA = A;
            if (piv_i != 0) {
                std::swap(permutedA[0], permutedA[piv_i]);
            }
            if (piv_j != 0) {
                for (size_t i = 0; i < n; ++i) {
                    bool temp = permutedA[i][0];
                    permutedA[i][0] = permutedA[i][piv_j];
                    permutedA[i][piv_j] = temp;
                }
            }

            size_t m = n - 1;
            BoolMatrix b(1, std::vector<bool>(m));
            BoolMatrix c(m, std::vector<bool>(1));
            BoolMatrix D(m, std::vector<bool>(m));

            for (size_t j = 0; j < m; ++j) b[0][j] = permutedA[0][j + 1];
            for (size_t i = 0; i < m; ++i) c[i][0] = permutedA[i + 1][0];
            for (size_t i = 0; i < m; ++i) {
                for (size_t j = 0; j < m; ++j) {
                    D[i][j] = permutedA[i + 1][j + 1];
                }
            }

            BoolMatrix cbT = strassenMultiply(c, b);
            BoolMatrix S = add(D, cbT);

            auto Sinv_opt = fastInvertRecursive(S);
            if (!Sinv_opt) return std::nullopt;
            const BoolMatrix& Sinv = *Sinv_opt;

            BoolMatrix t = strassenMultiply(Sinv, c);
            BoolMatrix uT = strassenMultiply(b, Sinv);
            BoolMatrix bT_t = strassenMultiply(b, t);
            bool alpha = true ^ bT_t[0][0];

            BoolMatrix permutedA_inv(n, std::vector<bool>(n));
            permutedA_inv[0][0] = alpha;
            for (size_t j = 0; j < m; ++j) permutedA_inv[0][j + 1] = uT[0][j];
            for (size_t i = 0; i < m; ++i) permutedA_inv[i + 1][0] = t[i][0];
            for (size_t i = 0; i < m; ++i) {
                for (size_t j = 0; j < m; ++j) {
                    permutedA_inv[i + 1][j + 1] = Sinv[i][j];
                }
            }

            BoolMatrix A_inv = permutedA_inv;
            if (piv_j != 0) {
                std::swap(A_inv[0], A_inv[piv_j]);
            }
            if (piv_i != 0) {
                for (size_t i = 0; i < n; ++i) {
                    bool temp = A_inv[i][0];
                    A_inv[i][0] = A_inv[i][piv_i];
                    A_inv[i][piv_i] = temp;
                }
            }
            
            return A_inv;
        }
    }

    inline std::optional<BoolMatrix> invertMatrixRobust(const BoolMatrix& A) {
        if (A.empty()) return BoolMatrix{};
        if (A.size() != A[0].size()) throw std::runtime_error("Matrix must be square.");
        return fastInvertRecursive(A);
    }

    inline BoolMatrix strassenMultiply(const BoolMatrix& A, const BoolMatrix& B) {
        size_t n1 = A.size();
        if (n1 == 0) return {};
        size_t p1 = A[0].size();
        size_t p2 = B.size();
        size_t n2 = B.empty() ? 0 : B[0].size();
        if (p1 != p2) throw std::runtime_error("Matrix dimensions are incompatible for multiplication.");

        if (n1 != p1 || n1 != n2 || n1 % 2 != 0 || n1 <= RECURSION_CUTOFF) {
            return standardMultiply(A, B);
        }

        size_t k = n1 / 2;
        BoolMatrix A11(k,std::vector<bool>(k)), A12(k,std::vector<bool>(k)), A21(k,std::vector<bool>(k)), A22(k,std::vector<bool>(k));
        BoolMatrix B11(k,std::vector<bool>(k)), B12(k,std::vector<bool>(k)), B21(k,std::vector<bool>(k)), B22(k,std::vector<bool>(k));
        
        split(A, A11, A12, A21, A22);
        split(B, B11, B12, B21, B22);
        
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
    
    inline void printMatrix(const BoolMatrix& A, std::ostream& out, const std::string& title = "") {
        if (!title.empty()) {
            out << title << " (" << A.size() << "x" << (A.empty() ? 0 : A[0].size()) << "):\n";
        }
        if (A.empty()) {
            out << "[Empty Matrix]\n";
            return;
        }
        for (const auto& row : A) {
            for (bool val : row) {
                out << val << " ";
            }
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

}

#endif 