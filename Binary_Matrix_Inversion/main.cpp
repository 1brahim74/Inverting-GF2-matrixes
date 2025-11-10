#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <optional>
#include <limits>
#include "strassen_inversion.hpp" // Changed from "robust_inversion.hpp"

// Use the MatrixOps namespace
using namespace MatrixOps;

// Function to clear cin buffer on invalid input
void clear_cin() {
    std::cin.clear();
    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
}

// Function to read a matrix from a file
std::optional<BoolMatrix> readMatrixFromFile(const std::string& filename) {
    std::ifstream inFile(filename);
    if (!inFile.is_open()) {
        std::cerr << "Error: Could not open file '" << filename << "'.\n";
        return std::nullopt;
    }

    BoolMatrix matrix;
    std::string line;
    size_t first_row_size = 0;

    while (std::getline(inFile, line)) {
        if (line.empty()) continue;

        std::stringstream ss(line);
        std::vector<bool> row;
        int val;
        while (ss >> val) {
            row.push_back(val != 0);
        }

        if (matrix.empty()) {
            first_row_size = row.size();
        } else if (row.size() != first_row_size) {
            std::cerr << "Error: Matrix is not rectangular.\n";
            return std::nullopt;
        }
        matrix.push_back(row);
    }

    if (matrix.empty()) {
        std::cerr << "Error: File is empty or contains no valid matrix data.\n";
        return std::nullopt;
    }
    if (matrix.size() != matrix[0].size()) {
        std::cerr << "Error: Matrix must be square.\n";
        return std::nullopt;
    }

    return matrix;
}

// Function to read a matrix from the terminal
std::optional<BoolMatrix> readMatrixFromTerminal() {
    std::cout << "\nEnter your matrix row by row (e.g., '1 0 1 0').\n";
    std::cout << "Press ENTER on an empty line when you are finished.\n";

    BoolMatrix matrix;
    std::string line;
    size_t first_row_size = 0;

    while (true) {
        std::cout << "> ";
        std::getline(std::cin, line);
        if (line.empty()) {
            break;
        }

        std::stringstream ss(line);
        std::vector<bool> row;
        int val;
        while (ss >> val) {
            if (val != 0 && val != 1) {
                 std::cerr << "Error: Please enter only 0s and 1s.\n";
                 return std::nullopt;
            }
            row.push_back(val != 0);
        }
        
        if (row.empty()) { // Could be a line with only whitespace
             if(matrix.empty()) {
                std::cout << "Input started with an empty line. Finishing input.\n";
                break;
             }
             continue;
        }

        if (matrix.empty()) {
            first_row_size = row.size();
        } else if (row.size() != first_row_size) {
            std::cerr << "Error: All rows must have the same number of columns.\n";
            return std::nullopt;
        }
        matrix.push_back(row);
    }

    if (matrix.empty()) {
        std::cerr << "Error: No matrix was entered.\n";
        return std::nullopt;
    }
    if (matrix.size() != matrix[0].size()) {
        std::cerr << "Error: Matrix must be square (rows must equal columns).\n";
        return std::nullopt;
    }

    return matrix;
}

int main() {
    const std::string MATRIX_FILE = "matrix.txt";
    const std::string ANSWER_FILE = "answer.txt";

    std::cout << "--- Robust Matrix Inverter using Strassen & LUP Decomposition ---\n";

    // --- Get Input Choice ---
    int input_choice = 0;
    while (input_choice != 1 && input_choice != 2) {
        std::cout << "\nSelect input method:\n";
        std::cout << "  1. Enter matrix from terminal\n";
        std::cout << "  2. Read matrix from '" << MATRIX_FILE << "'\n";
        std::cout << "Enter choice (1 or 2): ";
        std::cin >> input_choice;
        if (!std::cin || (input_choice != 1 && input_choice != 2)) {
            std::cerr << "Invalid input. Please enter 1 or 2.\n";
            clear_cin();
            input_choice = 0;
        }
    }
    clear_cin(); // Clear the rest of the line after reading the number

    // --- Load Matrix ---
    std::optional<BoolMatrix> matrix_opt;
    if (input_choice == 1) {
        matrix_opt = readMatrixFromTerminal();
    } else {
        matrix_opt = readMatrixFromFile(MATRIX_FILE);
    }

    if (!matrix_opt) {
        std::cerr << "\nMatrix loading failed. Exiting.\n";
        return 1;
    }

    std::cout << "\nMatrix loaded successfully.\n";

    // --- Get Output Choice ---
    int output_choice = 0;
    while (output_choice != 1 && output_choice != 2) {
        std::cout << "\nSelect output method for the result:\n";
        std::cout << "  1. Print to terminal\n";
        std::cout << "  2. Save to '" << ANSWER_FILE << "'\n";
        std::cout << "Enter choice (1 or 2): ";
        std::cin >> output_choice;
        if (!std::cin || (output_choice != 1 && output_choice != 2)) {
            std::cerr << "Invalid input. Please enter 1 or 2.\n";
            clear_cin();
            output_choice = 0;
        }
    }

    // --- Invert Matrix ---
    std::cout << "\nInverting matrix...\n";
    auto inverse_opt = invertMatrixRobust(*matrix_opt);

    // --- Handle and Display Output ---
    if (inverse_opt) {
        std::cout << "Inversion successful!\n";
        if (output_choice == 1) {
            printMatrix(*inverse_opt, std::cout, "\nComputed Inverse:");
        } else {
            std::ofstream outFile(ANSWER_FILE);
            if (!outFile) {
                std::cerr << "Error: Could not open '" << ANSWER_FILE << "' for writing.\n";
                return 1;
            }
            printMatrix(*inverse_opt, outFile);
            std::cout << "Result saved to '" << ANSWER_FILE << "'.\n";
        }
    } else {
        std::cerr << "Error: The matrix is singular and cannot be inverted.\n";
        if (output_choice == 2) {
            std::ofstream outFile(ANSWER_FILE);
            if (!outFile) {
                 std::cerr << "Error: Could not open '" << ANSWER_FILE << "' for writing.\n";
                return 1;
            }
            outFile << "Matrix is singular.\n";
            std::cout << "Status saved to '" << ANSWER_FILE << "'.\n";
        }
    }

    std::cout << "\nProcess complete.\n";
    return 0;
}