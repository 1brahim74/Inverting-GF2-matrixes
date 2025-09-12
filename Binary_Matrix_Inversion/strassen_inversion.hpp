#ifndef STRASSEN_INVERSION_HPP
#define STRASSEN_INVERSION_HPP

//
// Copyright (c) 2023, Ibrahim Mammadov
// This code is a combination and refactoring of two original files.
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

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <optional>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iomanip>
#include <numeric>

namespace MatrixOps {

// --- Type Definitions ---
using IntMatrix = std::vector<std::vector<int>>;
using BoolMatrix = std::vector<std::vector<bool>>;

struct InversionResult {
    std::optional<BoolMatrix> inverse;
    size_t original_size;
    size_t padded_size;
};


// --- Internal Helper Functions (Anonymous Namespace) ---
namespace {
    // --- Helpers for Invertible Matrix Generation ---
    IntMatrix createMatrix(int n, int m) {
        return IntMatrix(n, std::vector<int>(m, 0));
    }

    IntMatrix multiplyIntMatrix(const IntMatrix& A, const IntMatrix& B) {
        int n = A.size();
        int m = B[0].size();
        int p = A[0].size();
        IntMatrix C = createMatrix(n, m);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                int sum = 0;
                for (int k = 0; k < p; k++) {
                    sum = (sum + (A[i][k] * B[k][j]));
                }
                C[i][j] = sum % 2;
            }
        }
        return C;
    }

    IntMatrix generateLowerTriangular(int n) {
        IntMatrix L = createMatrix(n, n);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j <= i; j++) {
                L[i][j] = (i == j) ? 1 : (std::rand() % 2);
            }
        }
        return L;
    }

    IntMatrix generateUpperTriangular(int n) {
        IntMatrix U = createMatrix(n, n);
        for (int i = 0; i < n; i++) {
            for (int j = i; j < n; j++) {
                U[i][j] = (i == j) ? 1 : (std::rand() % 2);
            }
        }
        return U;
    }

    // --- Helpers for Strassen Matrix Inversion (Reformatted for Readability) ---
    const size_t STRASSEN_CUTOFF = 4;

    // Forward declarations for recursive functions
    std::optional<BoolMatrix> strassenInverse(const BoolMatrix& A);
    BoolMatrix strassenMultiply(const BoolMatrix& A, const BoolMatrix& B);

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

    BoolMatrix subtract(const BoolMatrix& A, const BoolMatrix& B) {
        return add(A, B);
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

    BoolMatrix combine(const BoolMatrix& A11, const BoolMatrix& A12, const BoolMatrix& A21, const BoolMatrix& A22) {
        size_t k = A11.size();
        size_t n = 2 * k;
        BoolMatrix A(n, std::vector<bool>(n));
        for (size_t i = 0; i < k; ++i) {
            for (size_t j = 0; j < k; ++j) {
                A[i][j] = A11[i][j];
                A[i][j + k] = A12[i][j];
                A[i + k][j] = A21[i][j];
                A[i + k][j + k] = A22[i][j];
            }
        }
        return A;
    }

    std::optional<BoolMatrix> gaussJordanInverse(BoolMatrix A) {
        size_t n = A.size();
        if (n == 0) return BoolMatrix{};
        
        BoolMatrix I(n, std::vector<bool>(n, false));
        for (size_t i = 0; i < n; ++i) I[i][i] = true;

        for (size_t j = 0; j < n; ++j) {
            size_t pivot_row = j;
            while (pivot_row < n && !A[pivot_row][j]) {
                pivot_row++;
            }
            if (pivot_row == n) return std::nullopt; // Singular
            std::swap(A[j], A[pivot_row]);
            std::swap(I[j], I[pivot_row]);
            for (size_t i = 0; i < n; ++i) {
                if (i != j && A[i][j]) {
                    for (size_t k = 0; k < n; ++k) {
                        A[i][k] = A[i][k] ^ A[j][k];
                        I[i][k] = I[i][k] ^ I[j][k];
                    }
                }
            }
        }
        return I;
    }

    BoolMatrix strassenMultiply(const BoolMatrix& A, const BoolMatrix& B) {
        size_t n = A.size();
        if (n <= STRASSEN_CUTOFF) {
            BoolMatrix C(n, std::vector<bool>(n, false));
            for (size_t i = 0; i < n; ++i) {
                for (size_t j = 0; j < n; ++j) {
                    bool sum = false;
                    for (size_t k = 0; k < n; ++k) {
                        sum = sum ^ (A[i][k] & B[k][j]);
                    }
                    C[i][j] = sum;
                }
            }
            return C;
        }

        size_t k = n / 2;
        BoolMatrix A11(k, std::vector<bool>(k)), A12(k, std::vector<bool>(k)), A21(k, std::vector<bool>(k)), A22(k, std::vector<bool>(k));
        BoolMatrix B11(k, std::vector<bool>(k)), B12(k, std::vector<bool>(k)), B21(k, std::vector<bool>(k)), B22(k, std::vector<bool>(k));

        split(A, A11, A12, A21, A22);
        split(B, B11, B12, B21, B22);

        BoolMatrix P1 = strassenMultiply(add(A11, A22), add(B11, B22));
        BoolMatrix P2 = strassenMultiply(add(A21, A22), B11);
        BoolMatrix P3 = strassenMultiply(A11, subtract(B12, B22));
        BoolMatrix P4 = strassenMultiply(A22, subtract(B21, B11));
        BoolMatrix P5 = strassenMultiply(add(A11, A12), B22);
        BoolMatrix P6 = strassenMultiply(subtract(A21, A11), add(B11, B12));
        BoolMatrix P7 = strassenMultiply(subtract(A12, A22), add(B21, B22));

        BoolMatrix C11 = add(subtract(add(P1, P4), P5), P7);
        BoolMatrix C12 = add(P3, P5);
        BoolMatrix C21 = add(P2, P4);
        BoolMatrix C22 = add(subtract(add(P1, P3), P2), P6);

        return combine(C11, C12, C21, C22);
    }

    std::optional<BoolMatrix> strassenInverse(const BoolMatrix& A) {
        size_t n = A.size();
        if (n == 0) return BoolMatrix{};
        if (n <= STRASSEN_CUTOFF) return gaussJordanInverse(A);

        size_t k = n / 2;
        BoolMatrix A11(k, std::vector<bool>(k)), A12(k, std::vector<bool>(k)), A21(k, std::vector<bool>(k)), A22(k, std::vector<bool>(k));
        split(A, A11, A12, A21, A22);

        auto A11_inv_opt = strassenInverse(A11);
        if (A11_inv_opt) {
            BoolMatrix A11_inv = *A11_inv_opt;
            BoolMatrix A21_A11_inv = strassenMultiply(A21, A11_inv);
            BoolMatrix S = subtract(A22, strassenMultiply(A21_A11_inv, A12));

            auto S_inv_opt = strassenInverse(S);
            if (!S_inv_opt) return std::nullopt;
            BoolMatrix S_inv = *S_inv_opt;

            BoolMatrix A11_inv_A12 = strassenMultiply(A11_inv, A12);
            BoolMatrix C12 = strassenMultiply(A11_inv_A12, S_inv);
            BoolMatrix C21 = strassenMultiply(S_inv, A21_A11_inv);
            BoolMatrix C11 = add(A11_inv, strassenMultiply(C12, A21_A11_inv));
            BoolMatrix C22 = S_inv;

            return combine(C11, C12, C21, C22);
        }

        auto A22_inv_opt = strassenInverse(A22);
        if (A22_inv_opt) {
            BoolMatrix A22_inv = *A22_inv_opt;
            BoolMatrix A12_A22_inv = strassenMultiply(A12, A22_inv);
            BoolMatrix S_prime = subtract(A11, strassenMultiply(A12_A22_inv, A21));

            auto S_prime_inv_opt = strassenInverse(S_prime);
            if (!S_prime_inv_opt) return std::nullopt;
            BoolMatrix S_prime_inv = *S_prime_inv_opt;

            BoolMatrix A22_inv_A21 = strassenMultiply(A22_inv, A21);
            BoolMatrix C11 = S_prime_inv;
            BoolMatrix C12 = strassenMultiply(S_prime_inv, A12_A22_inv);
            BoolMatrix C21 = strassenMultiply(A22_inv_A21, S_prime_inv);
            BoolMatrix C22 = add(A22_inv, strassenMultiply(A22_inv_A21, C12));

            return combine(C11, C12, C21, C22);
        }

        return std::nullopt;
    }

    size_t nextPowerOfTwo(size_t n) {
        if (n == 0) return 1;
        if ((n & (n - 1)) == 0) return n; // Already a power of two
        size_t p = 1;
        while (p < n) {
            p <<= 1;
        }
        return p;
    }
} // end anonymous namespace

// --- Public API Functions ---
void printBoolMatrix(const BoolMatrix& A, std::ostream& out) {
    if (A.empty()) { out << "[Empty Matrix]\n"; return; }
    for (const auto& row : A) {
        for (bool val : row) { out << val << " "; }
        out << "\n";
    }
}

std::optional<IntMatrix> generateInvertibleMatrix(int n) {
    if (n <= 0) { return std::nullopt; }
    std::srand(static_cast<unsigned int>(std::time(0)));
    IntMatrix L = generateLowerTriangular(n);
    IntMatrix U = generateUpperTriangular(n);
    return multiplyIntMatrix(L, U);
}

InversionResult invertMatrix(const BoolMatrix& A) {
    size_t original_size = A.size();
    if (original_size == 0) {
        return { BoolMatrix{}, 0, 0 };
    }
    size_t padded_size = nextPowerOfTwo(original_size);
    std::optional<BoolMatrix> inverse_opt;
    if (original_size == padded_size) {
        inverse_opt = strassenInverse(A);
    } else {
        BoolMatrix padded_A(padded_size, std::vector<bool>(padded_size, false));
        for (size_t i = 0; i < original_size; ++i) {
            for (size_t j = 0; j < original_size; ++j) {
                padded_A[i][j] = A[i][j];
            }
        }
        for (size_t i = original_size; i < padded_size; ++i) {
            padded_A[i][i] = true;
        }
        auto padded_inv_opt = strassenInverse(padded_A);
        if (padded_inv_opt) {
            BoolMatrix result(original_size, std::vector<bool>(original_size));
            for (size_t i = 0; i < original_size; ++i) {
                for (size_t j = 0; j < original_size; ++j) {
                    result[i][j] = (*padded_inv_opt)[i][j];
                }
            }
            inverse_opt = result;
        }
    }
    return { inverse_opt, original_size, padded_size };
}

std::optional<BoolMatrix> readMatrixFromFile(const std::string& filename) {
    std::ifstream inFile(filename);
    if (!inFile.is_open()) { return std::nullopt; }
    BoolMatrix A;
    std::string line;
    while (std::getline(inFile, line)) {
        std::stringstream ss(line);
        int val;
        std::vector<bool> row;
        while (ss >> val) { row.push_back(val != 0); }
        if (!row.empty()) { A.push_back(row); }
    }
    if (A.empty() || A.size() != A[0].size()) { return std::nullopt; }
    return A;
}

template<typename MatrixType>
bool saveMatrixToFile(const MatrixType& matrix, const std::string& filename) {
    std::ofstream outFile(filename);
    if (!outFile.is_open()) { return false; }
    for (const auto& row : matrix) {
        for (size_t j = 0; j < row.size(); ++j) {
            outFile << row[j] << (j == row.size() - 1 ? "" : " ");
        }
        outFile << "\n";
    }
    return true;
}

} // namespace MatrixOps

#endif // STRASSEN_INVERSION_HPP