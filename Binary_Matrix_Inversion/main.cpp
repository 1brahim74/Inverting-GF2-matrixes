#include <iostream>
#include <limits>
#include <string>
#include <optional>
#include <vector> // Include vector for explicit type usage
#include "strassen_inversion.hpp"

/**
 * @brief Clears the input buffer to handle invalid input gracefully.
 */
void clear_cin() {
    std::cin.clear();
    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
}

int main() {
    // Define the filenames that will be used.
    const std::string matrixFile = "matrix.txt";
    const std::string inverseFile = "answer.txt";

    // NOTE: We assume your library has these type aliases.
    // If not, replace MatrixOps::BoolMatrix with std::vector<std::vector<bool>>
    // and MatrixOps::IntMatrix with std::vector<std::vector<int>>.
    using MatrixOps::BoolMatrix;
    using MatrixOps::IntMatrix;

    std::cout << "--- Matrix Inverter using Strassen's Algorithm ---\n";
    std::cout << "--------------------------------------------------\n";

    // --- Step 1: Obtain the matrix ---
    int choice = 0;
    while (choice != 1 && choice != 2) {
        std::cout << "Select an option:\n";
        std::cout << "  1. Invert a matrix from '" << matrixFile << "'\n";
        std::cout << "  2. Generate a new invertible matrix and invert it\n";
        std::cout << "Enter choice (1 or 2): ";
        std::cin >> choice;
        if (!std::cin || (choice != 1 && choice != 2)) {
            std::cerr << "Invalid input. Please enter 1 or 2.\n";
            clear_cin();
            choice = 0;
        }
    }
    clear_cin(); // Clear buffer after reading choice

    // This will hold the final matrix to be inverted, which MUST be a BoolMatrix.
    std::optional<BoolMatrix> matrixToInvertOpt;
    int n = 0;

    if (choice == 1) {
        std::cout << "\n--- Loading Matrix ---\n";
        matrixToInvertOpt = MatrixOps::readMatrixFromFile(matrixFile);
        if (!matrixToInvertOpt) {
            std::cerr << "Failed to read matrix from '" << matrixFile 
                      << "'. The file might be missing, empty, or malformed." << std::endl;
            return 1;
        }
    } 
    else { // choice == 2
        std::cout << "\n--- Generating Matrix ---\n";
        std::cout << "Enter the desired size of the square matrix (e.g., 5, 8, 25): ";
        std::cin >> n;

        // Input validation for size
        while (!std::cin || n <= 0) {
            std::cerr << "Error: Invalid input. Please enter a positive integer: ";
            clear_cin();
            std::cin >> n;
        }
        clear_cin();

        std::cout << "Generating a " << n << "x" << n << " invertible matrix...\n";
        // generateInvertibleMatrix returns an optional<IntMatrix>
        auto generatedIntMatrixOpt = MatrixOps::generateInvertibleMatrix(n);
        if (!generatedIntMatrixOpt) {
            std::cerr << "Failed to generate matrix. The requested size might be invalid." << std::endl;
            return 1;
        }
        
        // **KEY CHANGE: Convert the IntMatrix to a BoolMatrix**
        IntMatrix intMatrix = *generatedIntMatrixOpt;
        BoolMatrix boolMatrix(n, std::vector<bool>(n));
        for(int i = 0; i < n; ++i) {
            for(int j = 0; j < n; ++j) {
                boolMatrix[i][j] = (intMatrix[i][j] != 0);
            }
        }
        matrixToInvertOpt.emplace(boolMatrix); // Place the converted matrix into our optional

        // Save the generated (and converted) matrix for inspection
        if (!MatrixOps::saveMatrixToFile(*matrixToInvertOpt, matrixFile)) {
            std::cerr << "Failed to save generated matrix to '" << matrixFile << "'." << std::endl;
            return 1;
        }
        std::cout << "Successfully generated and saved matrix to '" << matrixFile << "'.\n";
    }

    // From this point on, we are guaranteed to be working with a BoolMatrix
    if (!matrixToInvertOpt) {
        std::cerr << "Error: No matrix was loaded or generated." << std::endl;
        return 1;
    }
    
    n = matrixToInvertOpt->size(); // Get the size from the successfully loaded/converted matrix

    // --- Step 2: Inversion ---
    std::cout << "\n--- Inverting Matrix ---\n";
    
    // Print the initial matrix if it's small enough to be readable
    if (n < 21) {
        std::cout << "\nInitial Matrix to Invert:\n";
        MatrixOps::printBoolMatrix(*matrixToInvertOpt, std::cout);
    }
    
    // Invert the matrix
    MatrixOps::InversionResult result = MatrixOps::invertMatrix(*matrixToInvertOpt);

    // Report padding if it occurred
    if (result.original_size != result.padded_size) {
        std::cout << "Note: Matrix was padded from " << result.original_size << "x" << result.original_size
                  << " to " << result.padded_size << "x" << result.padded_size << " for the inversion process.\n";
    }

    // Check if the inversion was successful
    if (result.inverse) {
        auto inverseMatrix = *result.inverse;

        // Conditionally print the result to console
        if (n < 21) {
            std::cout << "\nComputed Inverse:\n";
            MatrixOps::printBoolMatrix(inverseMatrix, std::cout);
        } else {
            std::cout << "\nInverse for large matrix computed. Skipping console output.\n";
        }

        // Save the inverse matrix to a file
        if (MatrixOps::saveMatrixToFile(inverseMatrix, inverseFile)) {
            std::cout << "Successfully computed the inverse and saved it to '" << inverseFile << "'.\n";
        } else {
            std::cerr << "Failed to save the inverse matrix to '" << inverseFile << "'." << std::endl;
            return 1;
        }
    } else {
        std::cerr << "\nError: The matrix is singular and cannot be inverted." << std::endl;
        std::ofstream outFile(inverseFile);
        if(outFile.is_open()) {
            outFile << "Matrix is singular.\n";
        }
        return 1; // Indicate failure
    }

    std::cout << "\nProcess complete.\n";
    return 0;
}