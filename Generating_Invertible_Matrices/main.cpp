#include <iostream>
#include <vector>
#include <string>
#include <limits>
#include "generate_all_invertible_matrices.hpp"

int main() {
    int n;
    std::cout << "--- Exhaustive GF(2) Invertible Matrix Generator ---\n";
    std::cout << "Enter matrix size n (recommended <= 4): ";
    std::cin >> n;

    if (!std::cin || n <= 0) {
        std::cerr << "Error: Size must be a positive integer." << std::endl;
        return 1;
    }

    if (n > 4) {
        std::cout << "Warning: n > 4 will take a very long time and consume significant memory.\nProceeding...\n";
    }

    std::cout << "\nGenerating component matrices and finding all unique combinations...\n";
    std::vector<ExhaustiveGenerator::Matrix> allMatrices = ExhaustiveGenerator::generateAll(n);

    if (allMatrices.empty() && n > 0) {
        std::cerr << "Generation failed. The size " << n << " is likely too large for this exhaustive method.\n";
        return 1;
    }
    
    // Save the results to a file
    std::string filename = "unique_GF2_matrices_n" + std::to_string(n) + ".txt";
    if (!ExhaustiveGenerator::saveMatricesToFile(allMatrices, filename)) {
        std::cerr << "Error: Could not open or write to output file '" << filename << "'." << std::endl;
        return 1;
    }

    // Calculate the theoretical number of matrices for the final report
    unsigned long long theoretical_count = 1;
    for (int i = 0; i < n; ++i) {
        // Using unsigned long long to prevent overflow
        unsigned long long power_of_2_n = (1ULL << n);
        unsigned long long power_of_2_i = (1ULL << i);
        theoretical_count *= (power_of_2_n - power_of_2_i);
    }

    std::cout << "\nDone!" << std::endl;
    std::cout << "Found " << allMatrices.size() << " unique invertible " << n << "x" << n << " matrices over GF(2)." << std::endl;
    std::cout << "Theoretical count: " << theoretical_count << std::endl;
    std::cout << "Results saved to '" << filename << "'" << std::endl;

    return 0;
}