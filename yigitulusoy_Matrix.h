#ifndef YIGITULUSOY_MATRIX_H
#define YIGITULUSOY_MATRIX_H

#include <iostream> // For std::cerr, std::cout
#include <string>   // For std::string

class Matrix {
public:
    int rows;
    int cols;
    double** data; // Dynamically allocated 2D array for matrix data

    // Constructor
    Matrix(int r, int c);

    // Destructor
    ~Matrix();

    // Copy Constructor (for deep copy)
    Matrix(const Matrix& other);

    // Assignment Operator (for deep copy)
    Matrix& operator=(const Matrix& other);

    // Matrix-Vector Multiplication: result = M * v
    // M is this matrix, v is a dynamically allocated vector
    double* multiply(const double* vec, int vec_size) const;

    // Static helper functions for vector operations (since we can't use std::vector)
    static double dotProduct(const double* v1, const double* v2, int size);
    static void normalizeVector(double* vec, int size);
    static double getVectorNorm(const double* vec, int size);
    static double* subtractVectors(const double* v1, const double* v2, int size);
    static double* createVector(int size, double initial_value = 0.0);
    static void copyVector(double* dest, const double* src, int size);
    static void deleteVector(double* vec);


    // Read matrix from file
    void readFromFile(const std::string& filename);

    // Print matrix (for debugging)
    void print() const;
};

#endif // MATRIX_H