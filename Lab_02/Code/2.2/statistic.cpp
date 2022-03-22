#include "statistic.h"
#include "measure.h"
#include "random.h"

using namespace std;

//~~~~~~~~~~~~~~~~ costruttore e distruttore ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Statistic::Statistic(){
  _rnd.Init();
};

Statistic::~Statistic(){};

//~~~~~~~~~~~~~~~~ metodi ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


void Statistic::blocking(int M, int N, Measure *measure, const char* filename){
  ofstream out;
  int L = M/N;
  out.open(filename, ios::out | ios::trunc);
  int dim = measure->get_dimension();                          // la dimensione della mia misura
  vector<double> sum(dim,0);
  vector<double> sum2(dim,0);

  for(int i=0; i<N; i++){
    vector<double> meas(dim,0);
    for(int k=0; k<L; k++){
      vector<double> onethrow = measure->get_measure();
      for(int j=0; j<dim; j++) meas.at(j) += onethrow.at(j);
    }
    for(int j=0; j<dim; j++){
      sum.at(j)+=meas.at(j)/L;
      sum2.at(j)+=(meas.at(j)/L)*(meas.at(j)/L);
      double errore = error(sum.at(j)/(i+1), sum2.at(j)/(i+1), i);
      if(i==N-1) out << sqrt(sum.at(j)/(i+1)) << "   " << errore / (2*sqrt(sum.at(j)/(i+1))) << endl;
    }
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
