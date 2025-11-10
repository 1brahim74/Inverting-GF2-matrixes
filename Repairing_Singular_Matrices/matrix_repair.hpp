#ifndef MATRIX_REPAIR_HPP
#define MATRIX_REPAIR_HPP

// --- INSTRUCTIONS ---
// For matrices larger than MAXN, this library requires the Boost library.
// Please install it (e.g., `sudo apt-get install libboost-all-dev`)
// You will also need to link against boost_system if required by your compiler setup.
// ---

#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <bitset>
#include <numeric>
#include <utility>
#include <algorithm>
#include <stdexcept>
#include <ctime>
#include <cstdlib>
#include <omp.h>

// For matrices larger than MAXN, we use boost::dynamic_bitset
#include <boost/dynamic_bitset.hpp>

namespace MatrixRepair {

/// Define a maximum matrix size for the high-performance std::bitset representation.
const int MAXN = 2048;

/// Type alias for the high-performance, fixed-size matrix
using StaticBinaryMatrix = std::vector<std::bitset<MAXN>>;

/// Type alias for the flexible, dynamic-size matrix
using DynamicBinaryMatrix = std::vector<boost::dynamic_bitset<>>;

/// A type alias for storing the coordinates of flipped bits.
using FlipCoordinates = std::vector<std::pair<int, int>>;


// --- Internal Template-Based Implementation ---
namespace internal {
    /**
     * @brief A struct to hold the detailed results of Gaussian elimination.
     */
    struct GaussianEliminationResult {
        int rank;
        std::vector<int> pivot_cols;
        std::vector<int> original_row_indices; // Tracks rows after swaps
    };

    /**
     * @brief Performs Gaussian elimination in parallel using OpenMP on a generic MatrixType.
     * @tparam MatrixType The type of the matrix (e.g., StaticBinaryMatrix or DynamicBinaryMatrix).
     */
    template <typename MatrixType>
    GaussianEliminationResult parallel_gaussian_elimination(const MatrixType& M_in, int N) {
        if (N == 0) return {0, {}, {}};
        
        MatrixType mat_gauss = M_in;
        std::vector<int> original_row_indices(N);
        std::iota(original_row_indices.begin(), original_row_indices.end(), 0);

        int rank = 0;
        std::vector<int> pivot_cols;
        
        for (int j = 0; j < N && rank < N; ++j) {
            int p = rank;
            while (p < N && !mat_gauss[p].test(j)) {
                p++;
            }

            if (p < N) {
                std::swap(mat_gauss[p], mat_gauss[rank]);
                std::swap(original_row_indices[p], original_row_indices[rank]);
                pivot_cols.push_back(j);
                
                #pragma omp parallel for
                for (int i = 0; i < N; ++i) {
                    if (i != rank && mat_gauss[i].test(j)) {
                        mat_gauss[i] ^= mat_gauss[rank];
                    }
                }
                rank++;
            }
        }
        return {rank, pivot_cols, original_row_indices};
    }
    
    /**
     * @brief Template implementation for repairing a matrix.
     */
    template <typename MatrixType>
    MatrixType repair_matrix_impl(const MatrixType& M_in, int N, FlipCoordinates& flips) {
        flips.clear();
        
        auto result = parallel_gaussian_elimination(M_in, N);
        int rank = result.rank;

        if (rank == N) {
            return M_in;
        }

        MatrixType mat_repaired = M_in;
        
        std::vector<int> dependent_rows;
        for (int i = rank; i < N; ++i) {
            dependent_rows.push_back(result.original_row_indices[i]);
        }

        std::vector<int> free_cols;
        for (int j = 0; j < N; ++j) {
            if (std::find(result.pivot_cols.begin(), result.pivot_cols.end(), j) == result.pivot_cols.end()) {
                free_cols.push_back(j);
            }
        }
        
        int deficiency = N - rank;
        for (int k = 0; k < deficiency; ++k) {
            int row_to_fix = dependent_rows[k];
            int col_to_flip = free_cols[k];
            
            mat_repaired[row_to_fix].flip(col_to_flip);
            flips.emplace_back(row_to_fix, col_to_flip);
        }
        
        return mat_repaired;
    }

} // end namespace internal


// --- Public API ---

/**
 * @brief Computes rank. Dispatches to the correct templated implementation.
 */
template <typename MatrixType>
int compute_rank(const MatrixType& M, int N) {
    return internal::parallel_gaussian_elimination(M, N).rank;
}

/**
 * @brief Repairs a matrix. Dispatches to the correct templated implementation.
 */
template <typename MatrixType>
MatrixType repair_matrix(const MatrixType& M_in, int N, FlipCoordinates& flips) {
    return internal::repair_matrix_impl(M_in, N, flips);
}

/**
 * @brief Generates a singular matrix. Chooses static or dynamic based on N.
 */
template <typename MatrixType>
MatrixType generate_singular_matrix_impl(int N) {
    static bool seeded = false;
    if (!seeded) {
        std::srand(static_cast<unsigned int>(std::time(0)));
        seeded = true;
    }

    MatrixType A(N);
    // For dynamic_bitset, rows must be resized
    if constexpr (std::is_same_v<MatrixType, DynamicBinaryMatrix>) {
        for(int i = 0; i < N; ++i) A[i].resize(N);
    }

    do {
        for (int i = 0; i < N; i++) {
            A[i].reset();
            for (int j = 0; j < N; j++) {
                if (std::rand() % 2) {
                    A[i][j] = 1;
                }
            }
        }
    } while (compute_rank(A, N) == N);

    return A;
}


// --- File I/O (Remains mostly the same, but must handle both matrix types) ---

/**
 * @brief Reads a matrix from a file, choosing static or dynamic based on size.
 * This is a factory function that returns a std::variant or similar would be cleaner,
 * but for this example we will handle the logic in main. Here we just read into strings.
 */
std::vector<std::string> read_lines_from_file(const std::string& filename, int& N_out) {
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
    // No MAXN check here anymore
    return lines;
}

// Function to parse lines into the chosen matrix type
void parse_lines_to_static_matrix(const std::vector<std::string>& lines, StaticBinaryMatrix& matrix, int N) {
    matrix.assign(N, std::bitset<MAXN>());
    for (int i = 0; i < N; ++i) {
        std::stringstream ss(lines[i]);
        for (int j = 0; j < N; ++j) {
            int bit;
            if (!(ss >> bit)) throw std::runtime_error("Malformed matrix file.");
            if (bit) matrix[i].set(j);
        }
    }
}

void parse_lines_to_dynamic_matrix(const std::vector<std::string>& lines, DynamicBinaryMatrix& matrix, int N) {
    matrix.assign(N, boost::dynamic_bitset<>(N));
    for (int i = 0; i < N; ++i) {
        std::stringstream ss(lines[i]);
        for (int j = 0; j < N; ++j) {
            int bit;
            if (!(ss >> bit)) throw std::runtime_error("Malformed matrix file.");
            if (bit) matrix[i][j] = 1;
        }
    }
}

/**
 * @brief Writes a matrix to a file. Templated to handle both types.
 */
template <typename MatrixType>
bool write_matrix_to_file(const MatrixType& M, int N, const std::string& filename) {
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