#include <iostream>
#include <vector>
#include <string>
#include <limits>
#include "matrix_repair.hpp" // CHANGED

int main() {
    std::cout << "--- Singular Matrix Generator and Repair Tool ---\n";

    // --- 1. Get user input and validate ---
    int n;
    std::cout << "Enter the desired matrix size n: ";
    std::cin >> n;

    if (!std::cin || n <= 0) {
        std::cerr << "Error: Matrix size must be a positive integer." << std::endl;
        return 1;
    }
    if (n > MatrixRepair::MAXN) {
        std::cerr << "Error: Matrix size " << n << " exceeds the maximum supported size of "
                  << MatrixRepair::MAXN << "." << std::endl;
        return 1;
    }

    // --- 2. Generate a singular matrix ---
    std::cout << "\nStep 1: Generating a singular " << n << "x" << n << " matrix..." << std::endl;
    auto singularMatrixOpt = MatrixRepair::generateSingularMatrix(n);
    
    if (!singularMatrixOpt) {
        std::cerr << "Error: Failed to generate a singular matrix." << std::endl;
        return 1;
    }
    
    const std::string singularFilename = "matrix.txt";
    if (MatrixRepair::saveMatrixToFile(*singularMatrixOpt, singularFilename)) {
        std::cout << "Successfully generated and saved singular matrix to '" << singularFilename << "'." << std::endl;
    } else {
        std::cerr << "Error: Failed to save singular matrix to file." << std::endl;
        return 1;
    }

    // --- 3. Repair the generated matrix ---
    std::cout << "\nStep 2: Repairing the singular matrix..." << std::endl;
    
    // Convert the IntMatrix to the BitsetMatrix required by the repair function
    MatrixRepair::BitsetMatrix matrixToRepair = MatrixRepair::convertToBitsetMatrix(*singularMatrixOpt);
    if (matrixToRepair.empty()) {
        std::cerr << "Error: Failed to convert matrix for repair process." << std::endl;
        return 1;
    }

    // Call the repair function
    MatrixRepair::RepairResult result = MatrixRepair::repairMatrix(matrixToRepair);

    // --- 4. Report the results ---
    std::cout << "\n--- Repair Report ---" << std::endl;
    std::cout << "Initial Rank: " << result.initial_rank << " / " << n << std::endl;
    std::cout << "Rank Deficiency: " << result.deficiency << std::endl;

    if (result.was_already_invertible) {
        std::cout << "The generated matrix was already invertible. No changes were made." << std::endl;
    } else {
        std::cout << "Flipped " << result.flips.size() << " bits to repair the matrix:" << std::endl;
        for (const auto& flip : result.flips) {
            std::cout << "  - Flipped bit at (row " << flip.first << ", col " << flip.second << ")" << std::endl;
        }
        std::cout << "Matrix is now guaranteed to be invertible." << std::endl;
    }
    
    // --- 5. Save the final repaired matrix ---
    const std::string repairedFilename = "answer.txt";
    if (MatrixRepair::saveMatrixToFile(result.repaired_matrix, n, repairedFilename)) {
        std::cout << "Repaired matrix saved to '" << repairedFilename << "'." << std::endl;
    } else {
        std::cerr << "Error: Failed to save repaired matrix to file." << std::endl;
        return 1;
    }

    std::cout << "\nProcess complete." << std::endl;

    return 0;
}