#include "yigitulusoy_Matrix.h"
#include <fstream>
#include <sstream>
#include <cmath> // For std::sqrt, std::abs
#include <iomanip> // For std::fixed, std::setprecision

// Constructor
Matrix::Matrix(int r, int c) : rows(r), cols(c) {
    data = new double*[rows];
    for (int i = 0; i < rows; ++i) {
        data[i] = new double[cols];
        for (int j = 0; j < cols; ++j) {
            data[i][j] = 0.0; // Initialize with zeros
        }
    }
}

// Destructor
Matrix::~Matrix() {
    for (int i = 0; i < rows; ++i) {
        delete[] data[i];
    }
    delete[] data;
    data = nullptr; // Good practice to nullify after deletion
}

// Copy Constructor
Matrix::Matrix(const Matrix& other) : rows(other.rows), cols(other.cols) {
    data = new double*[rows];
    for (int i = 0; i < rows; ++i) {
        data[i] = new double[cols];
        for (int j = 0; j < cols; ++j) {
            data[i][j] = other.data[i][j];
        }
    }
}

// Assignment Operator
Matrix& Matrix::operator=(const Matrix& other) {
    if (this == &other) {
        return *this;
    }

    // Deallocate old memory
    for (int i = 0; i < rows; ++i) {
        delete[] data[i];
    }
    delete[] data;

    // Copy dimensions
    rows = other.rows;
    cols = other.cols;

    // Allocate new memory and copy data
    data = new double*[rows];
    for (int i = 0; i < rows; ++i) {
        data[i] = new double[cols];
        for (int j = 0; j < cols; ++j) {
            data[i][j] = other.data[i][j];
        }
    }
    return *this;
}

// Matrix-Vector Multiplication: result = M * v
double* Matrix::multiply(const double* vec, int vec_size) const {
    if (cols != vec_size) {
        std::cerr << "Error: Matrix-vector multiplication dimension mismatch." << std::endl;
        return nullptr; // Indicate error
    }

    double* result = new double[rows];
    for (int i = 0; i < rows; ++i) {
        result[i] = 0.0; // Initialize each element to 0
        for (int j = 0; j < cols; ++j) {
            result[i] += data[i][j] * vec[j];
        }
    }
    return result;
}

// Read matrix from file
void Matrix::readFromFile(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        // In a real application, you might throw an exception here.
        return;
    }

    std::string line;
    int currentRow = 0;
    while (std::getline(file, line) && currentRow < rows) {
        std::stringstream ss(line);
        double val;
        int currentCol = 0;
        while (ss >> val && currentCol < cols) {
            data[currentRow][currentCol] = val;
            currentCol++;
        }
        if (currentCol != cols) {
            std::cerr << "Warning: Incomplete row " << currentRow << " in " << filename << std::endl;
        }
        currentRow++;
    }
    if (currentRow != rows) {
        std::cerr << "Warning: Incomplete matrix in " << filename << std::endl;
    }
    file.close();
}

// Dot product of two vectors
double Matrix::dotProduct(const double* v1, const double* v2, int size) {
    if (!v1 || !v2) {
        std::cerr << "Error: Null vector passed to dotProduct." << std::endl;
        return 0.0;
    }
    double sum = 0.0;
    for (int i = 0; i < size; ++i) {
        sum += v1[i] * v2[i];
    }
    return sum;
}

// Normalize a vector
void Matrix::normalizeVector(double* vec, int size) {
    if (!vec) {
        std::cerr << "Error: Null vector passed to normalizeVector." << std::endl;
        return;
    }
    double norm = getVectorNorm(vec, size);
    if (norm == 0.0) return; // Avoid division by zero
    for (int i = 0; i < size; ++i) {
        vec[i] /= norm;
    }
}

// Get vector norm (Euclidean norm)
double Matrix::getVectorNorm(const double* vec, int size) {
    if (!vec) {
        std::cerr << "Error: Null vector passed to getVectorNorm." << std::endl;
        return 0.0;
    }
    double sum_sq = 0.0;
    for (int i = 0; i < size; ++i) {
        sum_sq += vec[i] * vec[i];
    }
    return std::sqrt(sum_sq);
}

// Subtract vectors
double* Matrix::subtractVectors(const double* v1, const double* v2, int size) {
    if (!v1 || !v2) {
        std::cerr << "Error: Null vector passed to subtractVectors." << std::endl;
        return nullptr;
    }
    double* result = new double[size];
    for (int i = 0; i < size; ++i) {
        result[i] = v1[i] - v2[i];
    }
    return result;
}

// Create a new dynamically allocated vector with an initial value
double* Matrix::createVector(int size, double initial_value) {
    double* vec = new double[size];
    for (int i = 0; i < size; ++i) {
        vec[i] = initial_value;
    }
    return vec;
}

// Copy contents of one vector to another
void Matrix::copyVector(double* dest, const double* src, int size) {
    if (!dest || !src) {
        std::cerr << "Error: Null vector passed to copyVector." << std::endl;
        return;
    }
    for (int i = 0; i < size; ++i) {
        dest[i] = src[i];
    }
}

// Delete a dynamically allocated vector
void Matrix::deleteVector(double* vec) {
    delete[] vec;
    vec = nullptr; // Good practice to nullify after deletion
}

// Print matrix (for debugging)
void Matrix::print() const {
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            std::cout << std::fixed << std::setprecision(4) << data[i][j] << "\t";
        }
        std::cout << std::endl;
    }
}