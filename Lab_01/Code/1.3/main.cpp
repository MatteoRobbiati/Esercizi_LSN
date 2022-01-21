#include <iostream>
#include "statistic.h"
#include "random.h"
#include <vector>

using namespace std;

int main(int argc, char* argv[]){

  if(argc < 2){
    cout << "You have to enter the number of throws!" << endl;
    return -1;
  }

  Statistic mystat;

  mystat.blocking(atoi(argv[1]), 100, "../../Results/buffon.dat");
  return 0;
}
