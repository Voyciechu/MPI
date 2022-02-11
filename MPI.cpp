// Zbi�r Mandelbrota - MPI - C++
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
 
struct Zespolona 
{
    double x, i;
   
    Zespolona(double x = 0, double i = 0) 
    {
        this->x = x;
        this->i = i;
    }
   
    Zespolona operator + (Zespolona input) 
    {
        return Zespolona (x + input.x, i + input.i);
    }
   
    Zespolona operator * (Zespolona input) 
    {
        return Zespolona (x * input.x - i * input.i, x * input.i + i * input.x);
    }

    double mod() 
    {
        return sqrt((x*x) + (i*i));
    }
};
 
int main() 
{
    MPI_Init(NULL, NULL);
 
    int n, id, width = 0, height = 0, N = 0;

    MPI_Comm_size(MPI_COMM_WORLD, &n);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
   
    MPI_Status status;
      
    if (id == 0) 
    {
        //Wczytanie danych z pliku przez główny wątek
        fstream input_file("mandelbrot.in", ios::in);
        input_file >> width >> height >> N;
        input_file.close();
       
        //Rozesłanie wczytanych danych do pozostałych wątków
        for (int i = 1; i < n; ++i) 
        {
            MPI_Send(&width, 1, MPI_INT, i, 286261, MPI_COMM_WORLD);
            MPI_Send(&height, 1, MPI_INT, i, 286261, MPI_COMM_WORLD);
            MPI_Send(&N, 1, MPI_INT, i, 286261, MPI_COMM_WORLD);
        }
    }
    else 
    {
        MPI_Recv(&width, 1, MPI_INT, 0, 286261, MPI_COMM_WORLD, &status);
        MPI_Recv(&height, 1, MPI_INT, 0, 286261, MPI_COMM_WORLD, &status);
        MPI_Recv(&N, 1, MPI_INT, 0, 286261, MPI_COMM_WORLD, &status);
    }
   
    int size = (width+1)*(height+1);
    char output[size];

    int part_start = id * size / n;
    int part_end = (id + 1) * size / n;

    for (int i = part_start; i < part_end; ++i) 
    {
        double x =  i % (width+1), y = i / (width+1);
        Zespolona z = Zespolona();
        Zespolona p = Zespolona(-2 + 3.0*x/width, 1 - 2.0*y/height);
       
        output[i] = 'X';
        for (int k = 0; k <= N; ++k) {
            if (z.mod() > 2) {
                output[i] = ' ';
                break;
            }
            z = z * z + p;
        }
    }
   
    if (id == 0) 
    {
        for (int i = 1; i < n; ++i) 
        {
            int start = i*size/n, meta = (i+1)*size/n;
            MPI_Recv(output+start, meta-start, MPI_CHAR, i, 286261, MPI_COMM_WORLD, &status);
        }
       
        fstream file("mandelbrot.out", ios::out);
        for (int i = 0; i < size; ++i) 
        {
            file << output[i];
            if (i % (width+1) == width) file << "\n";
        }
        file.close();
    }
    else 
    {
        int start = id*size/n, meta = (id+1)*size/n;
        MPI_Send(output+start, meta-start, MPI_CHAR, 0, 286261, MPI_COMM_WORLD);
    }
 
    MPI_Finalize();
    return 0;
}
