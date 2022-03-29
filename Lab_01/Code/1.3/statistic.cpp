#include "statistic.h"
#include "random.h"

using namespace std;

//~~~~~~~~~~~~~~~~ costruttore e distruttore ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Statistic::Statistic(){
  _rnd.Init();
};

Statistic::~Statistic(){};

//~~~~~~~~~~~~~~~~ metodi ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

void Statistic::blocking(int M, int N, const char* filename){
  ofstream out;
  out.open(filename, ios::out | ios::trunc);
  vector<double> sum;
  vector<double> sum2;

  for(int i=0; i<N; i++){
    vector<double> new_block = block_step(M/N);
    sum.resize(new_block.size());
    sum2.resize(new_block.size());
    for(unsigned int j=0; j<new_block.size(); j++){
      sum.at(j)+=new_block.at(j);
      sum2.at(j)+=pow(new_block.at(j),2);
      out << sum.at(j)/(i+1) << "   " << error(sum.at(j)/(i+1), sum2.at(j)/(i+1), i) << "   ";
    }
    out << endl;
  }
  out.close();
  return;
}


// ultimo esercizio in cui uso block_step
// la misura in questo caso è la stima del pigreco fatta L volte

vector<double> Statistic::block_step(int L){
  vector<double> measures;
  double d = 5;
  double l = 3;

  int counter = 0;
  for(int i=0; i<L; i++){
    double y_1 = _rnd.Rannyu(0,d);
    double sin = _rnd.generate_sin_angle();              // rapporto tra y e raggio della circ
                                                         // concentrica a quella unitaria
    double y_2 = y_1+sin*l;
    if(y_2>d || y_2<0.) counter++;
  }

  measures.push_back(2.*l*L/(double(counter)*d));
  return measures;
}


void Statistic::blocking_on_vector(vector<double> vec, const char* filename){
  ofstream out;
  out.open(filename, ios::out | ios::trunc);
  double sum=0;
  double sum2=0;

  for(unsigned int i=0; i<vec.size(); i++){
      sum+=vec.at(i);
      sum2+=pow(vec.at(i),2);
      out << sum/(i+1) << "   " << error(sum/(i+1), sum2/(i+1), i) << endl;
  }
  out.close();
  return;
}
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ funzioni generiche ~~~~~~~~~~~~~~~~


double Statistic::error(double val, double val2, unsigned int k){
  if(k==0) return 0;
  else return sqrt((val2-val*val)/k);
};


double Statistic::mean(vector<double> vec){
  int sum = 0;
  for (unsigned int i=0; i<vec.size(); i++) sum+=vec.at(i);
  return sum/vec.size();
}

double Statistic::chiquad(vector<double> vec, double EV, bool approx){
  ofstream out;
  out.open("../../Results/chiquad.dat");
  double chiq=0;
  if(approx==true){
    cout << "Approximated evaluation of chi_squared is starting. " << endl;
    cout << "We will use Expected Value as Variance in this approximation." << endl;
    for(unsigned int i=0; i<vec.size(); i++){
      chiq+=pow((vec.at(i)-EV),2)/EV;
      out << vec.at(i) << endl;
    }
  }
  else{
    cout << "Standard evaluation of chi_squared is starting." << endl;
    cout << "If you want to evaluate chi_sqared on this way you have to implement a method." << endl;
  }
  out.close();
  return chiq;
}




//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ TLC ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

void Statistic::TLC(int M, const char* distribution, vector<double> par, const char* filename){

  int N[] = {1,2,10,100};
  ofstream out;
  out.open(filename);
  for(int i=0; i<M; i++){                               // M realizzazioni
    double temp=0;
    for(int j=0; j<4; j++){
      for(int k=0; k<N[j]; k++){                        // ottenute sommando N[j] elementi
        temp+=_rnd.generate_by_inversion(distribution, par);
      }
      out << temp/N[j] << "    ";
    }
    out << endl;
  }
  out.close();
}
