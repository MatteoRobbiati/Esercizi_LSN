#include "mpi.h"
#include <iostream>

using namespace std;

int main(int argc, char* argv[]){
  int size, rank;
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  cout << "ehi ciao sto soltanto eseguendo gli ordini ciao ciao ciao sto allungando la frase per vedere se ci sono problemi " << endl;

  MPI_Finalize();
  return 0;
}
