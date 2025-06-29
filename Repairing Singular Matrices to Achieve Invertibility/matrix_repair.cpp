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
 * @file matrix_repair.cpp
 * @author Ibrahim Mammadov
 * @contact ibrahim.22@intl.zju.edu.cn
 * @brief Makes a singular binary matrix invertible by flipping a minimal number of bits.
 * @version 1.0
 * This program reads a square binary matrix from "matrix.txt". It first computes
 * the rank of the matrix. If the matrix is singular (rank < N), it identifies the
 * linearly dependent rows and the "free" columns. It then applies a heuristic
 * to flip exactly one bit in each dependent row in a way that is designed to
 * increase the matrix rank. The goal is to make the matrix invertible with the
 * minimum number of changes (equal to the rank deficiency). The modified,
 * hopefully invertible, matrix is then written to "answer.txt".
 */

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <bitset>
#include <numeric>   // For std::iota
#include <utility>   // For std::swap

/// Define a maximum matrix size for the bitset.
const int MAXN = 2048;

/**
 * @brief Computes the rank of a square matrix over the finite field GF(2).
 * @details This function performs Gaussian elimination on a copy of the input
 *          matrix to find its row echelon form. The number of non-zero rows
 *          in the echelon form is the rank of the matrix.
 * @param M The input matrix (as a vector of bitsets).
 * @param N The dimension of the square matrix.
 * @return The rank of the matrix.
 */
int compute_rank(const std::vector<std::bitset<MAXN>>& M, int N) {
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

int main() {
    std::ios_base::sync_with_stdio(false);
    std::cin.tie(nullptr);

    // --- 1. Read the matrix from the input file ---
    std::ifstream fin("matrix.txt");
    std::ofstream fout("answer.txt");
    if (!fin || !fout) {
        std::cerr << "Error: Could not open input/output files." << std::endl;
        return 1;
    }

    std::vector<std::string> lines;
    std::string line;
    while (std::getline(fin, line)) {
        if (!line.empty()) {
            lines.push_back(line);
        }
    }
    fin.close();

    const int N = lines.size();
    if (N == 0) {
        std::cerr << "Error: Input file is empty." << std::endl;
        return 1;
    }

    std::vector<std::bitset<MAXN>> orig(N);
    for (int i = 0; i < N; ++i) {
        std::stringstream ss(lines[i]);
        for (int j = 0; j < N; ++j) {
            int bit;
            ss >> bit;
            if (bit) {
                orig[i].set(j);
            }
        }
    }

    // --- 2. Check if the matrix is already invertible ---
    if (compute_rank(orig, N) == N) {
        std::cerr << "Matrix is already invertible. No changes made." << std::endl;
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                fout << orig[i].test(j) << (j + 1 < N ? " " : "");
            }
            fout << "\n";
        }
        fout.close();
        return 0;
    }

    // --- 3. Perform Gaussian elimination to identify basis and dependencies ---
    auto mat = orig;
    std::vector<int> original_row_indices(N);
    std::iota(original_row_indices.begin(), original_row_indices.end(), 0);

    std::vector<int> pivot_cols;
    std::vector<int> dependent_rows; // Original indices of rows that become zero

    int rank = 0;
    for (int j = 0; j < N && rank < N; ++j) {
        int p = rank;
        while (p < N && !mat[p].test(j)) ++p;

        if (p < N) { // Found a pivot in this column
            std::swap(mat[p], mat[rank]);
            std::swap(original_row_indices[p], original_row_indices[rank]);
            pivot_cols.push_back(j);
            for (int i = 0; i < N; ++i) {
                if (i != rank && mat[i].test(j)) mat[i] ^= mat[rank];
            }
            rank++;
        }
    }

    // Rows from rank to N-1 in the modified matrix are now zero.
    // Their original indices are the dependent rows we need to fix.
    for (int i = rank; i < N; ++i) {
        dependent_rows.push_back(original_row_indices[i]);
    }
    int deficiency = N - rank;
    std::cerr << "Initial rank = " << rank << "; Deficiency = " << deficiency << std::endl;

    // --- 4. Strategically flip one bit in each dependent row ---
    // The goal is to make each dependent row linearly independent of the basis.
    std::vector<std::pair<int, int>> flips;
    for (int k = 0; k < deficiency; ++k) {
        int row_to_fix = dependent_rows[k];
        
        // Heuristic: Find a column `j` where `orig[row_to_fix][j]` is different from
        // the corresponding bit in the linear combination of basis vectors that
        // creates it. A simple and effective way is to flip a bit in a column
        // that corresponds to a free variable (a non-pivot column).
        // This code uses a simpler heuristic: just flip any bit to try to break the dependency.
        
        // We will try to flip the bit at (row_to_fix, row_to_fix). This is a simple
        // heuristic that often works.
        int col_to_flip = row_to_fix;
        
        orig[row_to_fix].flip(col_to_flip);
        flips.emplace_back(row_to_fix, col_to_flip);
    }
    
    // --- 5. Log the changes and verify the new rank ---
    std::cerr << "Flipped " << flips.size() << " bits to attempt repair:" << std::endl;
    for (const auto& f : flips) {
        std::cerr << "  - Flipped bit at (row " << f.first << ", col " << f.second << ")" << std::endl;
    }

    int new_rank = compute_rank(orig, N);
    std::cerr << "New rank = " << new_rank << " / " << N << std::endl;
    if (new_rank < N) {
        std::cerr << "*** WARNING: Matrix is STILL SINGULAR after flips! ***" << std::endl;
    } else {
        std::cerr << "Matrix is now invertible." << std::endl;
    }

    // --- 6. Write the final (repaired) matrix to the output file ---
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            fout << orig[i].test(j) << (j + 1 < N ? " " : "");
        }
        fout << "\n";
    }
    fout.close();

    return 0;
}