/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
//parameters, observables
#include <string>
using namespace std;

const int m_props=4;
const double K_B = 1.380649e-23;
const double sigma = 0.34e-9;
const double mass = 39.948*1.66054e-27;
const double epsilon_on_K_B = 120;

int n_props;
int iv,ik,it,ie;
double stima_pot, stima_kin, stima_etot, stima_temp;

// averages
double acc,att;

//configuration
const int m_part=108;
double x[m_part],y[m_part],z[m_part],xold[m_part],yold[m_part],zold[m_part];
double vx[m_part],vy[m_part],vz[m_part];

// thermodynamical state
int npart;
double energy,temp,vol,rho,box,rcut;
string restart;    // useful for restarting with r(t-dt) and r(t)

// simulation
int nstep, iprint, seed;
double delta;

//blocking

void blocking_on_MD(int, int, string);
double error(double, double, unsigned int);

//functions
void Equilibrate_system();
void set_restart(string, string);
void rescale_velocities();
void Input(void);
void Move(void);
void ConfFinal(void);
void ConfXYZ(int);
void Measure(bool);
double Force(int, int);
double Pbc(double);
/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
