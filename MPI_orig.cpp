// Zbiór Mandelbrota - MPI - C++
// Uzywane polecenie do kompilacji pliku: mpicxx kolokwium.cpp
// Uzywane polecenie do uruchamiania programu: mpiexec -np 6 --hostfile nodes.txt ./a.out
// Wejscie: mandelbrot.in
// Wyjscie: mandelbrot.out
 
#include <iostream>
#include <vector>
#include <mpi.h>
#include <stdio.h>
#include <cmath>
#include <fstream>
using namespace std;
 
struct Complex {
    double a, b;
   
    Complex(double a = 0, double b = 0) {
        this->a = a;
        this->b = b;
    }
   
    double mod() {
        return sqrt(a*a + b*b);
    }
   
    Complex operator + (Complex z) {
        return Complex(a+z.a, b+z.b);
    }
   
    Complex operator * (Complex z) {
        return Complex(a*z.a - b*z.b, a*z.b + b*z.a);
    }
};
 
int main() {
    MPI_Init(NULL, NULL);
 
    int n;
    MPI_Comm_size(MPI_COMM_WORLD, &n);
 
    int id;
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
   
    MPI_Status status;
   
    int W = 0, H = 0, N = 0;
   
    if (id == 0) {
        fstream file("mandelbrot.in", ios::in);
        file >> W >> H >> N;
        file.close();
       
        for (int i = 1; i < n; ++i) {
            MPI_Send(&W, 1, MPI_INT, i, 13, MPI_COMM_WORLD);
            MPI_Send(&H, 1, MPI_INT, i, 13, MPI_COMM_WORLD);
            MPI_Send(&N, 1, MPI_INT, i, 13, MPI_COMM_WORLD);
        }
    }
    else {
        MPI_Recv(&W, 1, MPI_INT, 0, 13, MPI_COMM_WORLD, &status);
        MPI_Recv(&H, 1, MPI_INT, 0, 13, MPI_COMM_WORLD, &status);
        MPI_Recv(&N, 1, MPI_INT, 0, 13, MPI_COMM_WORLD, &status);
    }
   
    int size = (W+1)*(H+1);
    char output[size];
   
    for (int i = id*size/n; i < (id+1)*size/n; ++i) {
        double x =  i % (W+1), y = i / (W+1);
        Complex z = Complex();
        Complex p = Complex(-2 + 3.0*x/W, 1 - 2.0*y/H);
       
        output[i] = 'X';
        for (int k = 0; k <= N; ++k) {
            if (z.mod() > 2) {
                output[i] = ' ';
                break;
            }
            z = z * z + p;
        }
    }
   
    if (id == 0) {
        for (int i = 1; i < n; ++i) {
            int start = i*size/n, meta = (i+1)*size/n;
            MPI_Recv(output+start, meta-start, MPI_CHAR, i, 13, MPI_COMM_WORLD, &status);
        }
       
        fstream file("mandelbrot.out", ios::out);
        for (int i = 0; i < size; ++i) {
            file << output[i];
            if (i % (W+1) == W) file << "\n";
        }
        file.close();
    }
    else {
        int start = id*size/n, meta = (id+1)*size/n;
        MPI_Send(output+start, meta-start, MPI_CHAR, 0, 13, MPI_COMM_WORLD);
    }
 
    MPI_Finalize();
}
