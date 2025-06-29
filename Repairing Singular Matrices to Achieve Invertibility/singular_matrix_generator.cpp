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
 * @file singular_matrix_generator.cpp
 * @author Ibrahim Mammadov
 * @contact ibrahim.22@intl.zju.edu.cn
 * @brief Generates a random, singular (non-invertible) square matrix over GF(2).
 *
 * This program is designed to create test cases for algorithms that handle
 * singular matrices, such as matrix repair or inversion algorithms that need to
 * detect and report failure. It repeatedly generates random binary matrices
 * until it finds one that is non-invertible (i.e., its rank is less than its
 * dimension). The resulting singular matrix is saved to "matrix.txt" for use
 * by other programs.
 */

#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <cstdlib>
#include <ctime>
#include <utility> // For std::swap

/**
 * @brief Defines a matrix as a 2D vector of integers.
 */
using Matrix = std::vector<std::vector<int>>;

/**
 * @brief Generates a random n x n binary matrix.
 * @param n The dimension of the square matrix.
 * @return An n x n matrix with random 0 or 1 entries.
 */
Matrix generateRandomMatrix(int n) {
    Matrix A(n, std::vector<int>(n, 0));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            A[i][j] = std::rand() % 2;
        }
    }
    return A;
}

/**
 * @brief Writes a matrix to a file.
 * @param A The matrix to write.
 * @param filename The name of the file to save the matrix to.
 * @return True on success, false on failure to open the file.
 */
bool writeMatrixToFile(const Matrix &A, const std::string &filename) {
    std::ofstream fout(filename);
    if (!fout) {
        std::cerr << "Error: Cannot open " << filename << " for writing." << std::endl;
        return false;
    }
    for (size_t i = 0; i < A.size(); ++i) {
        for (size_t j = 0; j < A[i].size(); ++j) {
            fout << A[i][j] << (j == A[i].size() - 1 ? "" : " ");
        }
        fout << "\n";
    }
    fout.close();
    return true;
}

/**
 * @brief Computes the rank of a binary matrix over GF(2).
 * @details Uses Gaussian elimination with modulo-2 arithmetic. The number of
 *          pivots found equals the rank of the matrix.
 * @param A The input matrix.
 * @return The rank of the matrix.
 */
int computeRank(const Matrix &A) {
    int n = A.size();
    if (n == 0) return 0;
    Matrix M = A; // Work on a copy
    int rank = 0;
    for (int col = 0; col < n && rank < n; ++col) {
        // Find a pivot row (a row with a '1' in the current column)
        int pivotRow = rank;
        while (pivotRow < n && M[pivotRow][col] == 0) {
            pivotRow++;
        }

        if (pivotRow < n) { // A pivot was found
            std::swap(M[pivotRow], M[rank]);
            // Eliminate all other '1's in this column by XORing rows
            for (int row = 0; row < n; ++row) {
                if (row != rank && M[row][col] == 1) {
                    for (int k = col; k < n; k++) {
                        M[row][k] ^= M[rank][k];
                    }
                }
            }
            rank++;
        }
    }
    return rank;
}

/**
 * @brief Checks if a square matrix is invertible over GF(2).
 * @details A square matrix is invertible if and only if it has full rank.
 * @param A The matrix to check.
 * @return True if the matrix is invertible, false otherwise.
 */
bool isInvertible(const Matrix &A) {
    if (A.empty()) return false;
    int n = A.size();
    return computeRank(A) == n;
}

/**
 * @brief Main function to drive the singular matrix generation.
 */
int main() {
    int n;
    std::cout << "Enter matrix size n: ";
    std::cin >> n;
    if (!std::cin || n <= 0) {
        std::cerr << "Error: Matrix size must be a positive integer." << std::endl;
        return 1;
    }

    // Seed the random number generator.
    std::srand(static_cast<unsigned int>(std::time(0)));

    Matrix A;
    std::cout << "Generating random matrices until a singular one is found..." << std::endl;

    // Generate random matrices until a non-invertible (singular) one is found.
    do {
        A = generateRandomMatrix(n);
    } while (isInvertible(A));

    std::cout << "A singular (non-invertible) " << n << "x" << n << " matrix has been generated." << std::endl;

    // Write the singular matrix to "matrix.txt" for compatibility with other tools.
    if (writeMatrixToFile(A, "matrix.txt")) {
        std::cout << "Matrix saved to matrix.txt" << std::endl;
    } else {
        return 1; // Exit with an error if file writing failed.
    }

    return 0;
}