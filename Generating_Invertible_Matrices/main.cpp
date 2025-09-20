#include <iostream>
#include <vector>
#include <string>
#include <chrono> // For timing
#include "generate_all_invertible_matrices.hpp"

int main() {
    int n;
    std::cout << "--- GF(2) Invertible Matrix Generator (Bruhat Decomposition Method) ---\n";
    std::cout << "Enter matrix size n (recommended <= 5): ";
    std::cin >> n;

    if (!std::cin || n <= 0) {
        std::cerr << "Error: Size must be a positive integer." << std::endl;
        return 1;
    }

    if (n > 5) {
        std::cout << "Warning: n > 5 will take an extremely long time.\nProceeding...\n";
    }

    std::cout << "\nGenerating all unique invertible matrices for n = " << n << "...\n";
    
    auto start = std::chrono::high_resolution_clock::now();
    std::vector<InvertibleMatrixGenerator::Matrix> allMatrices = InvertibleMatrixGenerator::generateAll(n);
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;


    if (allMatrices.empty() && n > 0) {
        std::cerr << "Generation failed. The size " << n << " is likely too large for this method.\n";
        return 1;
    }
    
    // Save the results to a file
    std::string filename = "unique_GF2_matrices_n" + std::to_string(n) + ".txt";
    if (!InvertibleMatrixGenerator::saveMatricesToFile(allMatrices, filename)) {
        std::cerr << "Error: Could not open or write to output file '" << filename << "'." << std::endl;
        return 1;
    }

    // Calculate the theoretical number of matrices for the final report
    unsigned long long theoretical_count = 1;
    for (int i = 0; i < n; ++i) {
        unsigned long long power_of_2_n = (1ULL << n);
        unsigned long long power_of_2_i = (1ULL << i);
        if (power_of_2_n > power_of_2_i) { // Avoid underflow for n > 63
             theoretical_count *= (power_of_2_n - power_of_2_i);
        } else {
            theoretical_count = 0; // Indicates overflow
            break;
        }
    }

    std::cout << "\nDone!" << std::endl;
    std::cout << "Generation took " << elapsed.count() << " seconds." << std::endl;
    std::cout << "Found " << allMatrices.size() << " unique invertible " << n << "x" << n << " matrices over GF(2)." << std::endl;
    std::cout << "Theoretical count: " << theoretical_count << std::endl;
    std::cout << "Results saved to '" << filename << "'" << std::endl;

    return 0;
}