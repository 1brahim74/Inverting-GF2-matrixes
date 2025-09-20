#include <iostream>
#include <vector>
#include <string>
#include <limits>
#include "generate_all_invertible_matrices.hpp"

int main(int argc, char** argv) {
    int n;
    std::cout << "--- Exhaustive GF(2) Invertible Matrix Generator ---\n";
    std::cout << "Enter matrix size n (recommended <= 5): ";
    std::cin >> n;

    if (!std::cin || n <= 0) {
        std::cerr << "Error: Size must be a positive integer." << std::endl;
        return 1;
    }

    // Choose generator: default = Bruhat-aware; pass "--brute" to use old version
    bool useBrute = false;
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--brute") useBrute = true;
    }

    if (n > 6) {
        std::cout << "Warning: n > 6 will be very heavy even with Bruhat constraints.\nProceeding...\n";
    }

    std::cout << "\nGenerating component matrices and finding all unique combinations...\n";
    std::vector<ExhaustiveGenerator::Matrix> allMatrices;

    if (useBrute) {
        std::cout << "[Mode] Brute force (P * L * U) without Bruhat zero constraints on L.\n";
        allMatrices = ExhaustiveGenerator::generateAll(n);
    } else {
        std::cout << "[Mode] Bruhat-aware (P * L * U) with inversion-forced zeros in L.\n";
        allMatrices = ExhaustiveGenerator::generateAllBruhat(n);
    }

    if (allMatrices.empty() && n > 0) {
        std::cerr << "Generation failed or was capped to avoid explosion. "
                     "Try smaller n or '--brute' only for very small n.\n";
        return 1;
    }
    
    // Save the results to a file
    std::string filename = std::string(useBrute ? "brute_" : "bruhat_")
                         + "unique_GF2_matrices_n" + std::to_string(n) + ".txt";
    if (!ExhaustiveGenerator::saveMatricesToFile(allMatrices, filename)) {
        std::cerr << "Error: Could not open or write to output file '" << filename << "'." << std::endl;
        return 1;
    }

    // Theoretical |GL(n,2)| = ∏_{k=0}^{n-1} (2^n - 2^k)
    unsigned long long theoretical_count = 1;
    for (int k = 0; k < n; ++k) {
        unsigned long long pow2n = (1ULL << n);
        unsigned long long pow2k = (1ULL << k);
        theoretical_count *= (pow2n - pow2k);
    }

    std::cout << "\nDone!" << std::endl;
    std::cout << "Found " << allMatrices.size() << " unique invertible "
              << n << "x" << n << " matrices over GF(2)." << std::endl;
    std::cout << "Theoretical count: " << theoretical_count << std::endl;
    std::cout << "Results saved to '" << filename << "'\n";

    return 0;
}
