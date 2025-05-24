Readme.txt

**Project: Computing Eigenvalues and Eigenvectors using Normalized Power Iteration and Inverse Iteration**

This program implements the Normalized Power Iteration and Inverse Iteration algorithms to compute the dominant and smallest eigenvalues (in magnitude) and their corresponding eigenvectors for given matrices. It also verifies the relationship between the smallest eigenvalue of a matrix and the dominant eigenvalue of its inverse.

**How to Compile and Run:**

1.  **Source Files:** Ensure you have the following source files in the same directory:
    * `yigitulusoy_Matrix.h`
    * `yigitulusoy_Matrix.cpp`
    * `yigitulusoy_main.cpp`
    * `A.txt` (input matrix A)
    * `A_inv.txt` (input matrix A_inv)

2.  **Compilation:**
    Use a C++ compiler with C++11 standard support. Open your terminal or command prompt, navigate to the directory containing your source files, and execute the following command:
    ```bash
    g++ -std=c++11 -o eigenvalue_solver yigitulusoy_main.cpp yigitulusoy_Matrix.cpp -I.
    ```
    * `-o eigenvalue_solver`: This specifies the name of the executable file.
    * `yigitulusoy_main.cpp yigitulusoy_Matrix.cpp`: These are your source code files to be compiled.
    * `-I.`: This tells the compiler to look for header files (like `yigitulusoy_Matrix.h`) in the current directory.

3.  **Execution:**
    After successful compilation, run the program from your terminal or command prompt with the required four command-line arguments:
    ```bash
    ./eigenvalue_solver A.txt A_inv.txt 1e-6 results.txt
    ```
    **Required Arguments:**
    * `<A_file>`: Path to the input file containing matrix A.
    * `<A_inv_file>`: Path to the input file containing matrix A_inv.
    * `<tolerance>`: The tolerance value for the convergence criterion of the iterative algorithms.
    * `<output_file>`: The name of the output text file where the results will be written.

**Note on Program Completeness:**
This program is complete and implements all the required functionalities as specified in the project description.