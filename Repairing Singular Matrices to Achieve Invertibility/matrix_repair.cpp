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
 * @file matrix_repair_guaranteed.cpp
 * @author Ibrahim Mammadov
 * @contact ibrahim.22@intl.zju.edu.cn
 * @version 1.0.0
 * @brief Makes a singular binary matrix invertible using a guaranteed, optimal algorithm.
 *
 * This program reads a square binary matrix from "matrix.txt". It first computes
 * the rank (r) and rank deficiency (d = n - r). It then makes the matrix
 * invertible by performing exactly d bit-flips.
 *
 * The algorithm is guaranteed to be correct and optimal. It identifies the d
 * dependent rows and d free columns from Gaussian elimination. It then creates a
 * one-to-one pairing between them and flips the bit at each intersection. This
 * systematically repairs the rank deficiency and is proven to result in an
 * invertible matrix.
 */

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <bitset>
#include <numeric>   // For std::iota
#include <utility>   // For std::swap
#include <algorithm> // For std::find

/// Define a maximum matrix size for the bitset.
const int MAXN = 2048;
constexpr const char* APP_VERSION = "1.0.0";

/**
 * @brief Computes the rank of a square matrix over the finite field GF(2).
 */
int compute_rank(const std::vector<std::bitset<MAXN>>& M, int N) {
    auto A = M;
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
    std::cout << "Guaranteed Matrix Repair | Version " << APP_VERSION << std::endl;
    std::cout << "------------------------------------------" << std::endl;
    
    std::ifstream fin("matrix.txt");
    if (!fin) {
        std::cerr << "Error: Could not open input file 'matrix.txt'." << std::endl;
        return 1;
    }
    
    std::vector<std::string> lines;
    std::string line;
    while (std::getline(fin, line)) {
        if (!line.empty()) lines.push_back(line);
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
            int bit; ss >> bit;
            if (bit) orig[i].set(j);
        }
    }

    // --- 1. Perform Gaussian elimination to identify pivots, dependent rows, and free columns ---
    auto mat = orig;
    std::vector<int> original_row_indices(N);
    std::iota(original_row_indices.begin(), original_row_indices.end(), 0);

    std::vector<int> pivot_cols;
    int rank = 0;
    for (int j = 0; j < N && rank < N; ++j) {
        int p = rank;
        while (p < N && !mat[p].test(j)) ++p;

        if (p < N) { // Found a pivot
            std::swap(mat[p], mat[rank]);
            std::swap(original_row_indices[p], original_row_indices[rank]);
            pivot_cols.push_back(j);
            for (int i = 0; i < N; ++i) {
                if (i != rank && mat[i].test(j)) mat[i] ^= mat[rank];
            }
            rank++;
        }
    }

    // --- 2. Check if the matrix is already invertible ---
    if (rank == N) {
        std::cerr << "Matrix is already invertible. No changes made." << std::endl;
        // Write original to answer.txt
        std::ofstream fout("answer.txt");
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) fout << orig[i].test(j) << (j + 1 < N ? " " : "");
            fout << "\n";
        }
        fout.close();
        return 0;
    }

    // --- 3. Identify Dependent Rows and Free Columns ---
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
    
    int deficiency = N - rank;
    std::cerr << "Initial rank = " << rank << "; Deficiency = " << deficiency << std::endl;
    
    // --- 4. Apply the Guaranteed Bit-Flip Algorithm ---
    // Pair the k-th dependent row with the k-th free column and flip the bit.
    std::vector<std::pair<int, int>> flips;
    for (int k = 0; k < deficiency; ++k) {
        int row_to_fix = dependent_rows[k];
        int col_to_flip = free_cols[k];
        
        orig[row_to_fix].flip(col_to_flip);
        flips.emplace_back(row_to_fix, col_to_flip);
    }
    
    // --- 5. Log changes and verify the new rank ---
    std::cerr << "Flipped " << flips.size() << " bits to repair the matrix:" << std::endl;
    for (const auto& f : flips) {
        std::cerr << "  - Flipped bit at (row " << f.first << ", col " << f.second << ")" << std::endl;
    }

    int new_rank = compute_rank(orig, N);
    std::cerr << "New rank = " << new_rank << " / " << N << std::endl;
    if (new_rank < N) {
        std::cerr << "*** FATAL ERROR: Algorithm failed to repair the matrix! ***" << std::endl;
    } else {
        std::cerr << "Matrix is now invertible." << std::endl;
    }

    // --- 6. Write the final (repaired) matrix ---
    std::ofstream fout("answer.txt");
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) fout << orig[i].test(j) << (j + 1 < N ? " " : "");
        fout << "\n";
    }
    fout.close();

    return 0;
}