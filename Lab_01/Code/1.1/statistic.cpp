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


vector<double> Statistic::block_step(int L){
  vector<double> measures;
  double meas = 0;
  double var = 0;
  for(int i=0; i<L; i++){
    double rand = _rnd.Rannyu();
    meas+=rand;
    var+= pow((rand-0.5),2);
  }
  measures.push_back(meas/L);
  measures.push_back(var/L);
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

double Statistic::uniform_sampling(int min, int max){
  return _rnd.Rannyu(min,max);
}

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
  double chiq=0;
  if(approx==true){
    for(unsigned int i=0; i<vec.size(); i++){
      chiq+=pow((vec.at(i)-EV),2)/EV;
    }
  }
  else cout << "Per applicare il chi-quadro senza approssimazione utilizzare il metodo che coinvolge le incertezze sulle misure " << endl;
  return chiq;
}
