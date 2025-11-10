#include "matrix_repair.hpp" // Use the new flexible header
#include <iostream>
#include <string>
#include <limits>
#include <chrono> // For timing

// Use the library's namespace
using namespace MatrixRepair;

constexpr const char* APP_VERSION = "2.0.0"; // Version updated
constexpr const char* INPUT_FILENAME = "matrix.txt";
constexpr const char* OUTPUT_FILENAME = "answer.txt";

void clear_cin() {
    std::cin.clear();
    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
}

// A generic function to run the repair process on any matrix type
template <typename MatrixType>
void process_matrix(MatrixType& matrix, int N) {
    std::cout << "Matrix size: " << N << "x" << N << std::endl;
    if constexpr (std::is_same_v<MatrixType, StaticBinaryMatrix>) {
        std::cout << "Using high-performance static matrix (N <= " << MAXN << ")." << std::endl;
    } else {
        std::cout << "Using flexible dynamic matrix (N > " << MAXN << ")." << std::endl;
    }

    auto start_time = std::chrono::high_resolution_clock::now();
    int initial_rank = compute_rank(matrix, N);
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> rank_duration = end_time - start_time;
    std::cout << "Initial rank computation took " << rank_duration.count() << " seconds." << std::endl;

    if (initial_rank == N) {
        std::cout << "\nMatrix is already invertible (Rank = " << N << "). No changes needed." << std::endl;
        if (write_matrix_to_file(matrix, N, OUTPUT_FILENAME)) {
             std::cout << "Original matrix copied to '" << OUTPUT_FILENAME << "'." << std::endl;
        }
        return;
    }

    int deficiency = N - initial_rank;
    std::cout << "\nInitial rank = " << initial_rank << "; Deficiency = " << deficiency << std::endl;
    std::cout << "Repairing the matrix..." << std::endl;

    FlipCoordinates flips;
    start_time = std::chrono::high_resolution_clock::now();
    MatrixType repaired_matrix = repair_matrix(matrix, N, flips);
    end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> repair_duration = end_time - start_time;
    std::cout << "Repair process took " << repair_duration.count() << " seconds." << std::endl;


    std::cout << "Flipped " << flips.size() << " bits to repair the matrix:" << std::endl;
    for (const auto& f : flips) {
        std::cout << "  - Flipped bit at (row " << f.first << ", col " << f.second << ")" << std::endl;
    }

    int new_rank = compute_rank(repaired_matrix, N);
    std::cout << "\nNew rank = " << new_rank << " / " << N << std::endl;
    if (new_rank < N) {
        std::cerr << "*** FATAL ERROR: Algorithm failed to repair the matrix! ***" << std::endl;
    } else {
        std::cout << "Success! The matrix is now invertible." << std::endl;
    }

    if (write_matrix_to_file(repaired_matrix, N, OUTPUT_FILENAME)) {
        std::cout << "Repaired matrix saved to '" << OUTPUT_FILENAME << "'." << std::endl;
    }
}


int main() {
    std::cout << "Guaranteed Matrix Repair | Version " << APP_VERSION << std::endl;
    std::cout << "------------------------------------------" << std::endl;
    
    int choice = 0;
    while (choice != 1 && choice != 2) {
        std::cout << "Select an option:\n";
        std::cout << "  1. Repair a matrix from '" << INPUT_FILENAME << "'\n";
        std::cout << "  2. Generate a new singular matrix and repair it\n";
        std::cout << "Enter choice (1 or 2): ";
        std::cin >> choice;
        if (!std::cin || (choice != 1 && choice != 2)) {
            std::cerr << "Invalid input. Please enter 1 or 2." << std::endl;
            clear_cin();
            choice = 0;
        }
    }
    clear_cin();

    int N = 0;

    try {
        if (choice == 1) {
            std::cout << "\nReading matrix from '" << INPUT_FILENAME << "'..." << std::endl;
            auto lines = read_lines_from_file(INPUT_FILENAME, N);
            
            if (N <= MAXN) {
                StaticBinaryMatrix matrix;
                parse_lines_to_static_matrix(lines, matrix, N);
                process_matrix(matrix, N);
            } else {
                DynamicBinaryMatrix matrix;
                parse_lines_to_dynamic_matrix(lines, matrix, N);
                process_matrix(matrix, N);
            }

        } else { // choice == 2
            std::cout << "\nEnter matrix size N to generate: ";
            std::cin >> N;
            if (!std::cin || N <= 0) {
                std::cerr << "Error: Matrix size must be a positive integer." << std::endl;
                return 1;
            }
            clear_cin();
            
            std::cout << "Generating a singular " << N << "x" << N << " matrix..." << std::endl;
            
            if (N <= MAXN) {
                auto matrix = generate_singular_matrix_impl<StaticBinaryMatrix>(N);
                if (!write_matrix_to_file(matrix, N, INPUT_FILENAME)) {
                    std::cerr << "Warning: Could not save the generated matrix." << std::endl;
                }
                process_matrix(matrix, N);
            } else {
                auto matrix = generate_singular_matrix_impl<DynamicBinaryMatrix>(N);
                if (!write_matrix_to_file(matrix, N, INPUT_FILENAME)) {
                    std::cerr << "Warning: Could not save the generated matrix." << std::endl;
                }
                process_matrix(matrix, N);
            }
        }
    } catch (const std::runtime_error& e) {
        std::cerr << e.what() << std::endl;
        return 1;
    }

    return 0;
}