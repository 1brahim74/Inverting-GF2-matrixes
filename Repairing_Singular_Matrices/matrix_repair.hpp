#ifndef MATRIX_REPAIR_HPP
#define MATRIX_REPAIR_HPP

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
 * @file matrix_repair.hpp
 * @author Ibrahim Mammadov
 * @contact ibrahim.22@intl.zju.edu.cn
 * @brief A header-only library for generating, repairing, and analyzing square binary matrices.
 * @version 1.1.0
 */

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <bitset>
#include <numeric>   // For std::iota
#include <utility>   // For std::swap, std::pair
#include <algorithm> // For std::find
#include <stdexcept> // For std::runtime_error
#include <ctime>     // For std::time
#include <cstdlib>   // For std::srand, std::rand

namespace MatrixRepair {

/// Define a maximum matrix size for the bitset representation.
const int MAXN = 2048;
/// A type alias for a binary matrix, using std::bitset for efficiency.
using BinaryMatrix = std::vector<std::bitset<MAXN>>;
/// A type alias for storing the coordinates of flipped bits.
using FlipCoordinates = std::vector<std::pair<int, int>>;


/**
 * @brief Computes the rank of a square matrix over the finite field GF(2).
 * @param M The matrix to analyze.
 * @param N The dimension of the square matrix.
 * @return The rank of the matrix.
 */
int compute_rank(const BinaryMatrix& M, int N) {
    if (N == 0) return 0;
    auto A = M; // Work on a copy
    int rank = 0;
    for (int j = 0; j < N && rank < N; ++j) {
        int pivot_row = rank;
        while (pivot_row < N && !A[pivot_row].test(j)) {
            pivot_row++;
        }
        if (pivot_row < N) {
            std::swap(A[pivot_row], A[rank]);
            for (int i = 0; i < N; ++i) {
                if (i != rank && A[i].test(j)) {
                    A[i] ^= A[rank];
                }
            }
            rank++;
        }
    }
    return rank;
}

/**
 * @brief Checks if a square matrix is invertible (i.e., has full rank) over GF(2).
 * @param M The matrix to check.
 * @param N The dimension of the square matrix.
 * @return True if the matrix is invertible, false otherwise.
 */
bool is_invertible(const BinaryMatrix& M, int N) {
    if (N == 0) return false;
    return compute_rank(M, N) == N;
}

/**
 * @brief Repairs a singular binary matrix to make it invertible.
 *
 * This function takes a singular (non-invertible) square binary matrix and
 * performs an optimal number of bit-flips to make it invertible. The algorithm
 * is guaranteed to succeed by flipping exactly d bits, where d is the rank
 * deficiency (d = N - rank).
 *
 * @param M_in The input matrix to repair.
 * @param N The dimension of the square matrix.
 * @param flips A reference to a vector that will be populated with the (row, col)
 *              coordinates of the bits that were flipped.
 * @return The repaired, now-invertible matrix. If the matrix is already
 *         invertible, it is returned unchanged and the flips vector is empty.
 */
BinaryMatrix repair_matrix(const BinaryMatrix& M_in, int N, FlipCoordinates& flips) {
    flips.clear();
    BinaryMatrix mat_repaired = M_in;
    BinaryMatrix mat_gauss = M_in;

    std::vector<int> original_row_indices(N);
    std::iota(original_row_indices.begin(), original_row_indices.end(), 0);

    std::vector<int> pivot_cols;
    int rank = 0;
    for (int j = 0; j < N && rank < N; ++j) {
        int p = rank;
        while (p < N && !mat_gauss[p].test(j)) ++p;

        if (p < N) { // Found a pivot
            std::swap(mat_gauss[p], mat_gauss[rank]);
            std::swap(original_row_indices[p], original_row_indices[rank]);
            pivot_cols.push_back(j);
            for (int i = 0; i < N; ++i) {
                if (i != rank && mat_gauss[i].test(j)) mat_gauss[i] ^= mat_gauss[rank];
            }
            rank++;
        }
    }

    if (rank == N) {
        return mat_repaired; // Already invertible, no flips needed.
    }

    // Identify Dependent Rows (from original indices) and Free Columns
    std::vector<int> dependent_rows;
    for (int i = rank; i < N; ++i) {
        dependent_rows.push_back(original_row_indices[i]);
    }

    std::vector<int> free_cols;
    for (int j = 0; j < N; ++j) {
        if (std::find(pivot_cols.begin(), pivot_cols.end(), j) == pivot_cols.end()) {
            free_cols.push_back(j);
        }
    }
    
    // Pair the k-th dependent row with the k-th free column and flip the bit.
    int deficiency = N - rank;
    for (int k = 0; k < deficiency; ++k) {
        int row_to_fix = dependent_rows[k];
        int col_to_flip = free_cols[k];
        
        mat_repaired[row_to_fix].flip(col_to_flip);
        flips.emplace_back(row_to_fix, col_to_flip);
    }
    
    return mat_repaired;
}


/**
 * @brief Generates a random n x n singular binary matrix.
 *
 * This function repeatedly generates random binary matrices until it finds one
 * that is non-invertible.
 *
 * @param N The dimension of the square matrix.
 * @return An n x n singular matrix.
 */
BinaryMatrix generate_singular_matrix(int N) {
    if (N <= 0 || N > MAXN) {
        throw std::invalid_argument("Matrix size N must be between 1 and MAXN.");
    }
    
    // Seed the random number generator if it hasn't been seeded.
    static bool seeded = false;
    if (!seeded) {
        std::srand(static_cast<unsigned int>(std::time(0)));
        seeded = true;
    }

    BinaryMatrix A(N);
    do {
        for (int i = 0; i < N; i++) {
            A[i].reset(); // Clear previous bits
            for (int j = 0; j < N; j++) {
                if (std::rand() % 2) {
                    A[i].set(j);
                }
            }
        }
    } while (is_invertible(A, N));

    return A;
}


/**
 * @brief Reads a square binary matrix from a file.
 * @param filename The path to the input file.
 * @param N_out A reference to an integer that will store the detected matrix size.
 * @return The matrix read from the file.
 * @throws std::runtime_error if the file cannot be opened or is malformed.
 */
BinaryMatrix read_matrix_from_file(const std::string& filename, int& N_out) {
    std::ifstream fin(filename);
    if (!fin) {
        throw std::runtime_error("Error: Could not open input file '" + filename + "'.");
    }
    
    std::vector<std::string> lines;
    std::string line;
    while (std::getline(fin, line)) {
        if (!line.empty()) lines.push_back(line);
    }
    fin.close();

    N_out = lines.size();
    if (N_out == 0) {
        throw std::runtime_error("Error: Input file '" + filename + "' is empty.");
    }
    if (N_out > MAXN) {
        throw std::runtime_error("Error: Matrix size exceeds MAXN.");
    }

    BinaryMatrix matrix(N_out);
    for (int i = 0; i < N_out; ++i) {
        std::stringstream ss(lines[i]);
        int n_cols = 0;
        for (int j = 0; j < N_out; ++j) {
            int bit;
            if (!(ss >> bit)) {
                throw std::runtime_error("Error: Malformed matrix file. Row " + std::to_string(i) + " is too short.");
            }
            if (bit) matrix[i].set(j);
            n_cols++;
        }
        if (n_cols != N_out) {
             throw std::runtime_error("Error: Matrix is not square. Row " + std::to_string(i) + " has " + std::to_string(n_cols) + " columns, expected " + std::to_string(N_out));
        }
    }
    return matrix;
}

/**
 * @brief Writes a binary matrix to a file.
 * @param M The matrix to write.
 * @param N The dimension of the square matrix.
 * @param filename The name of the file to save the matrix to.
 * @return True on success, false on failure.
 */
bool write_matrix_to_file(const BinaryMatrix& M, int N, const std::string& filename) {
    std::ofstream fout(filename);
    if (!fout) {
        std::cerr << "Error: Cannot open " << filename << " for writing." << std::endl;
        return false;
    }
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            fout << M[i].test(j) << (j + 1 < N ? " " : "");
        }
        fout << "\n";
    }
    fout.close();
    return true;
}

} // namespace MatrixRepair

#endif // MATRIX_REPAIR_HPP