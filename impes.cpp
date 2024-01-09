#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <fstream>
#include <string>

#include "shl_v2.cpp"
#include "shl_wh.cpp"
#include "dshl_v2.cpp"
#include "dshl2_v2.cpp"
#include "we_v2.cpp"
#include "pts_v2.cpp"

vector<vector<double>> init_Matrix(int n1, int n2){
    vector<vector<double>> M(n1, vector<double>(n2));
    for (int i = 0; i < n1; i++){
        for (int j = 0; j < n2; j++){
            M[i][j] = 0.0;
        }
    }
    return M;
}

vector<double> init_Vector(int n){
    vector<double> F(n);
    for (int i = 0; i < n; i++) F[i] = 0.0;
    return F;
}

void print_Matrix(vector<vector<double>> M, int dim){
    for (int i = 0; i < dim; i++){
        for (int j = 0; j < dim; j++){
            cout << fixed << setprecision(6) << M[i][j] << "\t";
        }
        cout << endl;
    }
}

void print_Vector(vector<double> F, int dim){
    for (int i = 0; i < dim; i++){
        cout << F[i] << "\n";
    }
}




double lambda_n(double Sn, double rho_n){
    // relative permeability of phase w
    // double k_rn = // TODO: Define function of Sw
    return 1.0;//k_rn / rho_n;
}

double lambda_w(double Sw, double rho_w){
    // relative permeability of phase w
    // double k_rw = // TODO: Define function of Sw
    return 1.0;//k_rw / rho_w;
}

double abs_perm(double xx){
    return 1.0;
}

double lambda_sum(double lambda_w, double lambda_n){
    return lambda_w + lambda_n;
}

double epsilon(double xx, double Sw, double rho_w, double rho_n){   // diffusive coefficient function
    double lamb_n = lambda_n(1.0-Sw,rho_n);
    double lamb_w = lambda_w(Sw,rho_w);
    return abs_perm(xx) * lambda_sum(lamb_w,lamb_n);
}





double f_pr(double xx, double rho_w, double rho_n){    // source term function
    double f_w = 0.0;
    double f_n = 0.0;

    return f_w/rho_w + f_n/rho_n;
}

// Function to perform Gaussian elimination and back substitution
void solveLinearSystem(vector<vector<double>>& A, vector<double>& b, vector<double>& x) {
    int n = A.size();

    // Forward elimination
    for (int i = 0; i < n; ++i) {
        // Make the diagonal element 1
        double divisor = A[i][i];
        for (int j = 0; j < n; ++j) {
            A[i][j] /= divisor;
        }
        b[i] /= divisor;

        // Eliminate other elements in the current column
        for (int k = 0; k < n; ++k) {
            if (k != i) {
                double factor = A[k][i];
                for (int j = 0; j < n; ++j) {
                    A[k][j] -= factor * A[i][j];
                }
                b[k] -= factor * b[i];
            }
        }
    }

    // Back substitution
    for (int i = 0; i < n; ++i) {
        x[i] = b[i];
    }
}

// Function to multiply a matrix and a vector
std::vector<double> matrixVectorMultiply(const std::vector<std::vector<double>>& matrix, const std::vector<double>& vector) {
    // Check if the matrix and vector dimensions are compatible for multiplication
    size_t numRows = matrix.size();
    size_t numCols = matrix[0].size();
    size_t vectorSize = vector.size();

    if (numCols != vectorSize) {
        // Dimensions are not compatible for multiplication
        throw std::invalid_argument("Matrix and vector dimensions are not compatible for multiplication");
    }

    // Resulting vector
    std::vector<double> result(numRows, 0.0);

    // Perform matrix-vector multiplication
    for (size_t i = 0; i < numRows; ++i) {
        for (size_t j = 0; j < numCols; ++j) {
            result[i] += matrix[i][j] * vector[j];
        }
    }

    return result;
}

// Function to add two matrices
std::vector<std::vector<double>> addMatrices(const std::vector<std::vector<double>>& matrix1, const std::vector<std::vector<double>>& matrix2) {
    // Check if the matrices have the same dimensions
    size_t numRows1 = matrix1.size();
    size_t numCols1 = matrix1[0].size();
    size_t numRows2 = matrix2.size();
    size_t numCols2 = matrix2[0].size();

    if (numRows1 != numRows2 || numCols1 != numCols2) {
        // Matrices have different dimensions, cannot perform addition
        throw std::invalid_argument("Matrices have different dimensions, cannot perform addition");
    }

    // Resulting matrix
    std::vector<std::vector<double>> result(numRows1, std::vector<double>(numCols1, 0.0));

    // Perform matrix addition
    for (size_t i = 0; i < numRows1; ++i) {
        for (size_t j = 0; j < numCols1; ++j) {
            result[i][j] = matrix1[i][j] + matrix2[i][j];
        }
    }

    return result;
}

// Function to add two vectors
std::vector<double> addVectors(const std::vector<double>& vector1, const std::vector<double>& vector2) {
    // Check if the vectors have the same size
    size_t size1 = vector1.size();
    size_t size2 = vector2.size();

    if (size1 != size2) {
        // Vectors have different sizes, cannot perform addition
        throw std::invalid_argument("Vectors have different sizes, cannot perform addition");
    }

    // Resulting vector
    std::vector<double> result(size1, 0.0);

    // Perform vector addition
    for (size_t i = 0; i < size1; ++i) {
        result[i] = vector1[i] + vector2[i];
    }

    return result;
}

int main(){

    // mesh
    int nel = 100;
    double a = 0.0;
    double b = 2.0;
    double h = (b-a)/nel;
    cout << "h = " << h << endl;

    // time
    double t = 0.0;
    double T = 1.0;
    double dt = h*h;
    double dt_save = 0.25;
    cout << "dt = " << dt << endl;

    // polynomials basis
    int k = 1;          // polynomial degree
    int np = k*nel+1;   // mesh total nodes
    int nen = k+1;      // number of element nodes
    int nint = k+1;     // number of integration points

    // Initializing the mesh and matrices
    vector<double> xl(np,0.0);
    xl[0] = a;
    for (int ii = 1; ii < xl.size(); ii++) xl[ii] = xl[ii-1] + h/(nen-1);

    // Global stiffines matrix
    vector<vector<double>> K(np, vector<double>(np));
    K = init_Matrix(np,np);
    // print_Matrix(K,np);

    // Global source vector
    vector<double> F(np);
    F = init_Vector(np);
    // print_Vector(F,np);



    // Solution vector - Sw
    vector<double> Sw_n(np);
    Sw_n = init_Vector(np);
    for (int ii = 1; ii < Sw_n.size(); ii++) Sw_n[ii] = 0.0;  
    // for (int ii = 0; ii < Sw_n.size(); ii++) cout << Sw_n[ii] << ",";
    // cout << endl;

    vector<double> Sw_np1(np);
    Sw_np1 = init_Vector(np);
    for (int ii = 1; ii < Sw_np1.size(); ii++) Sw_np1[ii] = 0.0;  
    // for (int ii = 0; ii < Sw_np1.size(); ii++) cout << Sw_np1[ii] << ",";
    // cout << endl;

    // Solution vector - pbar
    vector<double> pbar_n(np);
    pbar_n = init_Vector(np);
    for (int ii = 1; ii < pbar_n.size(); ii++) pbar_n[ii] = 0.0;  
    // for (int ii = 0; ii < pbar_n.size(); ii++) cout << pbar_n[ii] << ",";
    // cout << endl;

    vector<double> pbar_np1(np);
    pbar_np1 = init_Vector(np);
    for (int ii = 1; ii < pbar_np1.size(); ii++) pbar_np1[ii] = 0.0;  
    // for (int ii = 0; ii < pbar_np1.size(); ii++) cout << pbar_np1[ii] << ",";
    // cout << endl;




    // Local matrix and vector (element)
    vector<vector<double>> Ke(nen, vector<double>(nen));
    Ke = init_Matrix(nen,nen);

    vector<double> Fe(nen);
    Fe = init_Vector(nen);



    // Basis functions and numerical integration points
    vector<vector<double>> shg(nint, vector<double>(nint));
    shg = init_Matrix(nint,nint);

    vector<vector<double>> shg_wh(nint, vector<double>(nint));
    shg_wh = init_Matrix(nint,nint);

    vector<vector<double>> dshg(nint, vector<double>(nint));
    dshg = init_Matrix(nint,nint);

    vector<vector<double>> dshg2(nint, vector<double>(nint));
    dshg2 = init_Matrix(nint,nint);

    vector<double> w(nint);
    w = init_Vector(nint);

    vector<double> pt(nint);
    pt = init_Vector(nint);

    shg = shl(nen,nint);
    // print_Matrix(shg,nen);
    // shg_wh = shl_wh(nen,nint,beta);
    // print_Matrix(shg_wh,nen);
    dshg = dshl(nen,nint);
    // print_Matrix(dshg,nen);
    dshg2 = dshl2(nen,nint);
    // print_Matrix(dshg2,nen);
    w = we(nint);
    // print_Vector(w,nint);
    pt = pts(nint);
    // print_Vector(pt,nint);



    // Solver - IMPES
    while (t <= T) { 
        // cout << "t = " << t << endl;  

        // RESOLVING FOR PRESSURE - pbar_np1    

        // restart global stiffines matrix and source vector
        K = init_Matrix(np,np);
        F = init_Vector(np);

        // Global problem construction
        for (int n = 0; n < nel; n++){
            // restart local (element) stiffines matrix and source vector
            Ke = init_Matrix(nen,nen);
            Fe = init_Vector(nen);

            for (int l = 0; l < nint; l++){

                // Evaluate x as function of reference element parameter, i.e., x(t)
                double xx = 0.0;
                for (int i = 0; i < nen; i++){
                    xx += shg[i][l]*xl[n*(nen-1) + i];
                }

                // Local source vector and stiffines matrix construction 
                for (int j = 0; j < nen; j++){
                    Fe[j] = Fe[j] + f_pr(xx, rho_w, rho_n)*shg[j][l]*w[l]*h/2.0; 

                    for (int i = 0; i < nen; i++){
                        // TODO: Evaluate Sw_mean elementwise    
                        Ke[i][j] = Ke[i][j] + epsilon(xx, Sw_mean, rho_w, rho_n)*(dshg[i][l]*2.0/h)*(dshg[j][l]*2.0/h)*w[l]*h/2.0;

                        // TODO: coding ...

                    }
                }
            }

            // TODO: coding ...
        }

        // TODO: coding ...

        // updating results and time step
        for (int ii = 0; ii < Sw_n.size(); ii++) Sw_n[ii] = Sw_np1[ii];
        for (int ii = 0; ii < pbar_n.size(); ii++) pbar_n[ii] = pbar_np1[ii];
        t += dt;

    } 

    cout << "OK!" << endl;

    return 0;
}