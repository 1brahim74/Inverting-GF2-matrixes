// you better use this command for running the code: g++ -std=c++17 -O3 -o generator_test main.cpp -fopenmp
#include <iostream>
#include <fstream>
#include <string>
#include <chrono>
#include <atomic> // For the thread-safe counter

#include "generate_all_invertible_matrices.hpp"

// Function to calculate the theoretical order of GL(n, 2)
unsigned long long get_theoretical_count(int n) {
    if (n < 0) return 0;
    unsigned long long count = 1;
    for (int i = 0; i < n; ++i) {
        unsigned long long power_of_2_n = (1ULL << n);
        unsigned long long power_of_2_i = (1ULL << i);
        if (power_of_2_n < power_of_2_i) { // Overflow check
             return 0;
        }
        count *= (power_of_2_n - power_of_2_i);
    }
    return count;
}

// Helper function to write a single matrix to an open file stream
void write_matrix_to_stream(const InvertibleMatrixGenerator::Matrix& matrix, std::ofstream& outFile) {
    for (const auto& row : matrix) {
        for (size_t j = 0; j < row.size(); ++j) {
            outFile << row[j] << (j == row.size() - 1 ? "" : " ");
        }
        outFile << "\n";
    }
    outFile << "\n"; // Blank line separator
}


int main() {
    int n;
    std::cout << "--- Invertible Matrix Toolkit (Iterator & Parallel) ---\n";
    std::cout << "Enter matrix size n (e.g., 3, 4, or 5): ";
    std::cin >> n;

    if (!std::cin || n <= 0) {
        std::cerr << "Error: Size must be a positive integer." << std::endl;
        return 1;
    }
    
    char save_choice;
    std::cout << "Do you want to save the matrices to a file? (y/n): ";
    std::cin >> save_choice;

    std::string filename;
    if (save_choice == 'y' || save_choice == 'Y') {
        filename = "unique_GF2_matrices_n" + std::to_string(n) + ".txt";
        std::cout << "Matrices will be saved to '" << filename << "'" << std::endl;
        if (n >= 5) {
            long double file_size_gb = (long double)get_theoretical_count(n) * (n * n + n) / (1024.0 * 1024.0 * 1024.0);
            std::cout << "Warning: The output file for n=" << n << " will be extremely large (approx. " 
                      << file_size_gb << " GB). Proceed with caution." << std::endl;
        }
    }
    
    unsigned long long theoretical_count = get_theoretical_count(n);
    std::cout << "Theoretical number of matrices for n=" << n << " is " << theoretical_count << std::endl;
    std::cout << "--------------------------------------------------------" << std::endl;
    std::cout << "Starting generation using the parallel method with " << omp_get_max_threads() << " threads...\n";

    std::atomic<long long> parallel_count = 0;
    auto start_par = std::chrono::high_resolution_clock::now();

    // Use a conditional to decide whether to open a file
    if (!filename.empty()) {
        std::ofstream outFile(filename);
        if (!outFile) {
            std::cerr << "Error: Could not open file for writing." << std::endl;
            return 1;
        }

        InvertibleMatrixGenerator::generateAllParallel(n, 
            [&](const InvertibleMatrixGenerator::Matrix& matrix) {
                parallel_count++;
                // This 'critical' section ensures only one thread writes to the file at a time.
                // It is necessary to prevent a corrupted file, but it will slow down
                // the process compared to just counting.
                #pragma omp critical
                {
                    write_matrix_to_stream(matrix, outFile);
                }
            });

    } else {
        // If not saving, just count as fast as possible.
        InvertibleMatrixGenerator::generateAllParallel(n, 
            [&](const InvertibleMatrixGenerator::Matrix& matrix) {
                parallel_count++;
            });
    }


    auto end_par = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_par = end_par - start_par;

    std::cout << "\nParallel generator finished." << std::endl;
    std::cout << "  - Matrices found: " << parallel_count << std::endl;
    std::cout << "  - Time taken: " << elapsed_par.count() << " seconds." << std::endl;
    std::cout << "  - Correctness: " << (parallel_count == theoretical_count ? "PASS" : "FAIL") << std::endl;
    
    if (!filename.empty()) {
        std::cout << "  - Results saved to '" << filename << "'" << std::endl;
    }
    std::cout << "--------------------------------------------------------" << std::endl;
    
    if (n >= 6) {
        std::cout << "\nNote: For n>=6, the process is possible due to memory efficiency, "
                  << "but it will take a very long time to complete." << std::endl;
    }

    return 0;
}
