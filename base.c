#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <time.h>
#include "mpi.h"

int main(int argc, char** argv)
{
    //Startup
    MPI_Init(&argc, &argv);

    //Data
    int rank, size, N=10, *tablica, x, reszta, po_ile, *tablica_buforowa, *tablica_wynikowa, minimum;
    MPI_File fh;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank); // Thread number
    MPI_Comm_size(MPI_COMM_WORLD, &size); // Threads amount
    
    MPI_File_open(MPI_COMM_WORLD, "dane.txt", MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &fh); 

    printf("Reszta %d\n", reszta);

    if(reszta == 0)
    {
        po_ile = N / size;
    } else 
    {
        po_ile = (N + size -reszta) /size;
    }

    tablica_buforowa = malloc(sizeof(int) * po_ile);
    
    if(rank == 0)
    {
        tablica = (int*)malloc(N * sizeof(int));
        tablica_wynikowa = (int*)malloc(size * sizeof(int));
        //Wypełnianie tablicap[x] danymi
        //Reszta kodu wątku głównego
    }
    
    MPI_Scatter(tablica, po_ile, MPI_INT, tablica_buforowa, po_ile, MPI_INT, 0, MPI_COMM_WORLD);

    //Mielenie
    MPI_File_write(fh, buf, 1000, MPI_INT, MPI_STATUS_IGNORE);
    MPI_File_close(&fh); 

    MPI_Gather(&minimum, 1,MPI_INT, tablica_wynikowa, 1,MPI_INT, 0, MPI_COMM_WORLD);
    
    if(rank == 0)
    {
        minimum = tablica_wynikowa[0];
        for(x=0; x < size; x +=1)
        {
            if(minimum > tablica_wynikowa[x])
                minimum = tablica_wynikowa[x];
        }
        printf("Minimum to: %d\n", minimum);
    }

    // Exit
    MPI_Finalize();
    return 0;
}

