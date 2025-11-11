#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

using namespace std;

// Read a GF(2) matrix from a text file: rows with space-separated 0/1.
vector<vector<int>> readMatrix(const string& filename) {
    ifstream fin(filename);
    if (!fin.is_open()) {
        cerr << "Error: cannot open file " << filename << endl;
        exit(1);
    }

    vector<vector<int>> mat;
    string line;
    while (getline(fin, line)) {
        // skip completely empty lines
        bool onlySpaces = true;
        for (char c : line) {
            if (!isspace(static_cast<unsigned char>(c))) {
                onlySpaces = false;
                break;
            }
        }
        if (onlySpaces) continue;

        stringstream ss(line);
        vector<int> row;
        int x;
        while (ss >> x) {
            x &= 1; // ensure in {0,1} mod 2
            row.push_back(x);
        }
        if (!row.empty())
            mat.push_back(row);
    }

    if (mat.empty()) {
        cerr << "Error: file " << filename << " is empty or invalid." << endl;
        exit(1);
    }

    // Check rectangular
    size_t cols = mat[0].size();
    for (size_t i = 1; i < mat.size(); i++) {
        if (mat[i].size() != cols) {
            cerr << "Error: non-rectangular matrix in " << filename << endl;
            exit(1);
        }
    }

    return mat;
}

vector<vector<int>> multiplyGF2(const vector<vector<int>>& A,
                                const vector<vector<int>>& B) {
    size_t n = A.size();
    size_t m = A[0].size();
    size_t p = B[0].size();

    if (B.size() != m) {
        cerr << "Error: incompatible dimensions for multiplication: "
             << "A is " << n << "x" << m
             << ", B is " << B.size() << "x" << p << endl;
        exit(1);
    }

    vector<vector<int>> C(n, vector<int>(p, 0));

    // Standard triple loop in GF(2)
    for (size_t i = 0; i < n; i++) {
        for (size_t k = 0; k < m; k++) {
            if (A[i][k] == 0) continue;
            // When A[i][k] = 1, add row k of B (mod 2)
            for (size_t j = 0; j < p; j++) {
                C[i][j] ^= B[k][j]; // addition mod 2
            }
        }
    }

    return C;
}

bool isIdentityGF2(const vector<vector<int>>& M) {
    size_t n = M.size();
    if (M[0].size() != n) return false;

    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            int expected = (i == j) ? 1 : 0;
            if ((M[i][j] & 1) != expected) return false;
        }
    }
    return true;
}

void printMatrix(const vector<vector<int>>& M, const string& name) {
    cout << name << ":\n";
    for (const auto& row : M) {
        for (size_t j = 0; j < row.size(); j++) {
            cout << (row[j] & 1);
            if (j + 1 < row.size()) cout << ' ';
        }
        cout << '\n';
    }
    cout << '\n';
}

int main() {
    // Read A and B
    vector<vector<int>> A = readMatrix("matrix.txt");
    vector<vector<int>> B = readMatrix("answer.txt");

    // Check they are square and same size (for inverse test)
    if (A.size() != A[0].size()) {
        cerr << "Error: matrix.txt is not square.\n";
        return 1;
    }
    if (B.size() != B[0].size()) {
        cerr << "Error: answer.txt is not square.\n";
        return 1;
    }
    if (A.size() != B.size()) {
        cerr << "Error: A and B are different sizes.\n";
        return 1;
    }

    // Compute AB and BA in GF(2)
    auto AB = multiplyGF2(A, B);
    auto BA = multiplyGF2(B, A);

    printMatrix(AB, "A * B (mod 2)");
    printMatrix(BA, "B * A (mod 2)");

    cout << "A * B is " << (isIdentityGF2(AB) ? "" : "NOT ")
         << "the identity matrix over GF(2).\n";
    cout << "B * A is " << (isIdentityGF2(BA) ? "" : "NOT ")
         << "the identity matrix over GF(2).\n";

    return 0;
}
