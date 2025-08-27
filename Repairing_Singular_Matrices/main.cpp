#include "matrix_repair.hpp"
#include <iostream>
#include <string>
#include <limits> // For numeric_limits

// Use the library's namespace
using namespace MatrixRepair;

constexpr const char* APP_VERSION = "1.1.0";
constexpr const char* INPUT_FILENAME = "matrix.txt";
constexpr const char* OUTPUT_FILENAME = "answer.txt";


/**
 * @brief Clears the input buffer to handle invalid input gracefully.
 */
void clear_cin() {
    std::cin.clear();
    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
}

/**
 * @brief Main function to run the matrix repair application.
 */
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
    clear_cin(); // Clear buffer after reading choice

    BinaryMatrix matrix;
    int N = 0;

    try {
        if (choice == 1) {
            std::cout << "\nReading matrix from '" << INPUT_FILENAME << "'..." << std::endl;
            matrix = read_matrix_from_file(INPUT_FILENAME, N);
        } else { // choice == 2
            std::cout << "\nEnter matrix size N to generate: ";
            std::cin >> N;
            if (!std::cin || N <= 0 || N > MAXN) {
                std::cerr << "Error: Matrix size must be a positive integer up to " << MAXN << "." << std::endl;
                return 1;
            }
            clear_cin();
            
            std::cout << "Generating a singular " << N << "x" << N << " matrix..." << std::endl;
            matrix = generate_singular_matrix(N);
            
            std::cout << "Saving generated matrix to '" << INPUT_FILENAME << "' for inspection." << std::endl;
            if (!write_matrix_to_file(matrix, N, INPUT_FILENAME)) {
                // Non-fatal error, we can still proceed
                std::cerr << "Warning: Could not save the generated matrix." << std::endl;
            }
        }
    } catch (const std::runtime_error& e) {
        std::cerr << e.what() << std::endl;
        return 1;
    }

    std::cout << "Matrix size: " << N << "x" << N << std::endl;

    // --- 1. Analyze and Repair the Matrix ---
    int initial_rank = compute_rank(matrix, N);
    
    if (initial_rank == N) {
        std::cout << "\nMatrix is already invertible (Rank = " << N << "). No changes needed." << std::endl;
        if (write_matrix_to_file(matrix, N, OUTPUT_FILENAME)) {
             std::cout << "Original matrix copied to '" << OUTPUT_FILENAME << "'." << std::endl;
        }
        return 0;
    }

    int deficiency = N - initial_rank;
    std::cout << "\nInitial rank = " << initial_rank << "; Deficiency = " << deficiency << std::endl;
    std::cout << "Repairing the matrix..." << std::endl;

    FlipCoordinates flips;
    BinaryMatrix repaired_matrix = repair_matrix(matrix, N, flips);

    // --- 2. Log changes and verify the new rank ---
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

    // --- 3. Write the final answer to the output file ---
    if (write_matrix_to_file(repaired_matrix, N, OUTPUT_FILENAME)) {
        std::cout << "Repaired matrix saved to '" << OUTPUT_FILENAME << "'." << std::endl;
    }

    return 0;
}