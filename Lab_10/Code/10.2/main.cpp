#include <iostream>
#include <vector>
#include <fstream>
#include <algorithm>
#include <string>
#include "random.h"
#include "TSP.h"
#include "mpi.h"

using namespace std;

int main(int argc, char *argv[]){

  // ~~~~~~~~~~~~~~~ Preparing the cities into the square ~~~~

  int N = 32;
  int M = 1000;        // 250 chromos for each walker (4x250=1000 as in Lab_09)
  int Nstep = 1000;
  int Nswap = 10;

  // ~~~~~~~~~~~~~~ some useful objects ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  vector<int> best_path_1(N);
  vector<int> best_path_2(N);

  int tag1, tag2;

  string square_circuit = "../../Results/square_circuit.dat";
  string guys[]   = {"Jim", "Al", "John", "Jack"};
  vector<int> idx = {0, 1, 2, 3};
  Random rnd;
  rnd.Init();
  rnd.generate_cities_in_square(1, N, square_circuit);


  // ~~~~~~~~~~~~~~~ introducing MPI ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  int size, rank;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Status stat;

  string id = to_string(rank);

  // ~~~~~~~~~~~~~~~ all the guys working on the same circuit ~~~~~~~~~~~~~~~~~~~~~~~
  // rank+1 is used for setting the random seed: i.e. srand(rank+1)
  // because i need that the four guys work with different pop of chromos

  Salesman *guy = new Salesman(N, M, square_circuit, rnd, rank+1);

  // ~~~~~~~~~~~~~~~ saving the four initial configuration ~~~~~~~~~~~~~~~~~~~~~~~~~

  guy->show_best_chromo();
  guy->save_chromo(0, "../../Results/first_found_by"+guys[rank]+".dat", true);

  // ~~~~~~~~~~~~~~ the parallel developement of the GA ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // The guys make a run of Nswaps moves, then they exchange the best chromo found
  // sharing with the others the progresses done

  for(int i=0; i<Nstep; i=i+Nswap){
    if(rank==0 && i%100==0) cout << "Generation " << i+1 << endl;
    //shuffle idx vector if rank=0, in order to change the perm of the exchanges each time
    if(rank==0) random_shuffle(idx.begin(), idx.end());
    MPI_Bcast(&idx.front(), 4, MPI_INT, 0, MPI_COMM_WORLD);      // sharing with the ranks
    guy->run(Nswap);                                             // the run of 10 steps
    //now i have done Nswap genetic step and i've sorted the population
    for(int j=0; j<2; j++){
      tag1 = j*2;
      tag2 = j*2+1;
      if(rank == idx[2*j]){
        best_path_1 = guy->get_best_path();
        MPI_Send(&best_path_1.front(), N, MPI_INT, idx[2*j+1], tag1, MPI_COMM_WORLD);
        MPI_Recv(&best_path_2.front(),  N, MPI_INT, idx[2*j+1], tag2, MPI_COMM_WORLD, &stat);
        guy->set_path(0, best_path_2);
      }
      if(rank == idx[2*j+1]){
        best_path_2 = guy->get_best_path();
        MPI_Recv(&best_path_1.front(), N, MPI_INT, idx[2*j], tag1, MPI_COMM_WORLD, &stat);
        MPI_Send(&best_path_2.front(), N, MPI_INT, idx[2*j], tag2, MPI_COMM_WORLD);
        guy->set_path(0, best_path_1);
      }
    }

    guy->sort_pop();
  }


  // ~~~~~~~~~~~~~~~ here i show the results ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  cout << "This is the result obtained by: " << guys[rank] << "(rank " << rank << ")" << endl;
  guy->show_best_chromo();
  guy->save_chromo(0, "../../Results/best_found_by" + guys[rank]+".dat", true);

  rnd.SaveSeed();

  MPI_Finalize();

  return 0;
}
