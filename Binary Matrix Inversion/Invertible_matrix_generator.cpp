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

/**
 * @file Invertible_matrix_generator.cpp
 * @author Ibrahim Mammadov
 * @contact ibrahim.22@intl.zju.edu.cn
 * @brief Generates a random, invertible square matrix over the finite field GF(2).
 * @version 1.0.0
 * This program creates an invertible matrix by generating a random lower-triangular
 * matrix (L) and a random upper-triangular matrix (U), both with ones on the
 * diagonal. The final matrix A = L * U is guaranteed to be invertible. The user
 * is prompted for a size, which must be a power of two to be compatible with
 * the accompanying Strassen inversion program. The resulting matrix is saved
 * to "matrix.txt".
 */

#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iomanip>

/**
 * @brief Defines a matrix as a 2D vector of integers.
 */
using Matrix = std::vector<std::vector<int>>;

/**
 * @brief Creates an n x m matrix initialized with zeros.
 * @param n The number of rows.
 * @param m The number of columns.
 * @return An n x m matrix filled with zeros.
 */
Matrix createMatrix(int n, int m) {
    return Matrix(n, std::vector<int>(m, 0));
}

/**
 * @brief Multiplies two matrices over the finite field GF(2).
 * @details Performs standard O(n^3) matrix multiplication, with all
 *          arithmetic performed modulo 2.
 * @param A The left-hand side matrix.
 * @param B The right-hand side matrix.
 * @return The resulting matrix C = (A * B) mod 2.
 */
Matrix multiplyMatrix(const Matrix& A, const Matrix& B) {
    int n = A.size();
    int m = B[0].size();
    int p = A[0].size();
    Matrix C = createMatrix(n, m);
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

/**
 * @brief Writes a matrix to a specified output stream.
 * @param A The matrix to write.
 * @param out The output stream (e.g., std::cout or an std::ofstream).
 */
void writeMatrixToStream(const Matrix& A, std::ostream& out) {
    for (size_t i = 0; i < A.size(); ++i) {
        for (size_t j = 0; j < A[i].size(); ++j) {
            out << A[i][j] << (j == A[i].size() - 1 ? "" : " ");
        }
        out << "\n";
    }
}

/**
 * @brief Generates a random lower-triangular matrix over GF(2).
 * @details The matrix has ones on the diagonal, which guarantees it is invertible.
 * @param n The size of the square matrix.
 * @return A random n x n lower-triangular matrix.
 */
Matrix generateLowerTriangular(int n) {
    Matrix L = createMatrix(n, n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j <= i; j++) {
            L[i][j] = (i == j) ? 1 : (std::rand() % 2);
        }
    }
    return L;
}

/**
 * @brief Generates a random upper-triangular matrix over GF(2).
 * @details The matrix has ones on the diagonal, which guarantees it is invertible.
 * @param n The size of the square matrix.
 * @return A random n x n upper-triangular matrix.
 */
Matrix generateUpperTriangular(int n) {
    Matrix U = createMatrix(n, n);
    for (int i = 0; i < n; i++) {
        for (int j = i; j < n; j++) {
            U[i][j] = (i == j) ? 1 : (std::rand() % 2);
        }
    }
    return U;
}

/**
 * @brief Main function to drive the matrix generation process.
 */
int main() {
    // Seed the random number generator.
    std::srand(static_cast<unsigned int>(std::time(0)));

    int n;
    std::cout << "Enter the size of the matrix (must be a power of 2): ";
    std::cin >> n;

    // Input validation: ensure n is a positive power of 2.
    if (!std::cin || n <= 0 || (n & (n - 1)) != 0) {
        std::cerr << "Error: Invalid size. Please enter a positive power of 2 (e.g., 2, 4, 8, 16...)." << std::endl;
        return 1;
    }

    // Generate lower and upper triangular matrices. Their product is always invertible.
    Matrix L = generateLowerTriangular(n);
    Matrix U = generateUpperTriangular(n);
    Matrix A = multiplyMatrix(L, U);

    // Write the final invertible matrix A to the output file.
    std::ofstream outFile("matrix.txt");
    if (!outFile.is_open()) {
        std::cerr << "Error opening output file: matrix.txt" << std::endl;
        return 1;
    }

    writeMatrixToStream(A, outFile);
    outFile.close();

    std::cout << "Successfully generated a " << n << "x" << n << " invertible matrix and saved it to matrix.txt." << std::endl;

    return 0;
}