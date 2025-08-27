#include <iostream>
#include <limits>
#include "strassen_inversion.hpp"

int main() {
    // Define the filenames that will be used.
    const std::string matrixFile = "matrix.txt";
    const std::string inverseFile = "answer.txt";

    // --- Get user input ---
    int n;
    std::cout << "--- Invertible Matrix Generator & Inverter ---\n";
    std::cout << "Enter the desired size of the square matrix (e.g., 5, 8, 25): ";
    std::cin >> n;

    // Input validation
    while (!std::cin) {
        std::cerr << "Error: Invalid input. Please enter a positive integer: ";
        std::cin.clear();
        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        std::cin >> n;
    }

    // --- Step 1: Generation ---
    std::cout << "\n--- Step 1: Generation ---\n";
    std::cout << "Generating a " << n << "x" << n << " invertible matrix...\n";
    
    auto generatedMatrixOpt = MatrixOps::generateInvertibleMatrix(n);
    if (!generatedMatrixOpt) {
        std::cerr << "Failed to generate matrix. Please provide a positive size." << std::endl;
        return 1;
    }

    // Save the generated matrix to a file
    if (!MatrixOps::saveMatrixToFile(*generatedMatrixOpt, matrixFile)) {
        std::cerr << "Failed to save generated matrix to '" << matrixFile << "'." << std::endl;
        return 1;
    }
    std::cout << "Successfully generated and saved matrix to '" << matrixFile << "'.\n";

    // --- Step 2: Inversion ---
    std::cout << "\n--- Step 2: Inversion ---\n";
    
    // Read the matrix back from the file
    auto matrixToInvertOpt = MatrixOps::readMatrixFromFile(matrixFile);
    if (!matrixToInvertOpt) {
        std::cerr << "Failed to read matrix from '" << matrixFile << "'. File might be missing or malformed." << std::endl;
        return 1;
    }
    std::cout << "Successfully read " << n << "x" << n << " matrix from '" << matrixFile << "'." << std::endl;

    // =========================================================
    // NEWLY ADDED BLOCK: Print the initial matrix if it's small
    // =========================================================
    if (n < 21) {
        std::cout << "\nInitial Matrix to Invert:\n";
        MatrixOps::printBoolMatrix(*matrixToInvertOpt, std::cout);
    }
    // =========================================================
    
    // Invert the matrix
    MatrixOps::InversionResult result = MatrixOps::invertMatrix(*matrixToInvertOpt);

    // Report padding if it occurred
    if (result.original_size != result.padded_size) {
        std::cout << "Note: Matrix was padded from " << result.original_size << "x" << result.original_size
                  << " to " << result.padded_size << "x" << result.padded_size << " for inversion process.\n";
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
        std::cerr << "\nThe matrix is singular and cannot be inverted." << std::endl;
        // Optionally write a message to the output file
        std::ofstream outFile(inverseFile);
        if(outFile.is_open()) {
            outFile << "Matrix is singular.\n";
        }
        return 1;
    }

    std::cout << "\nProcess complete.\n";
    return 0;
}