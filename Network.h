/* Network.h */
#include <algorithm>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <ctime>
 #include <sstream>
#include <iomanip>
#include <new>
using namespace std;

int N = 500; // Network size
int T = 300000; // number of time steps
double t = 0.0; // time (initialized to zero)
// object to define each neuron

struct Neuron {
double v; // voltage values
double* tsp_excite;
 int tsp_Ecount;
vector <double> tspN; // spike times of this Neuron, for raster plot
vector <int> connections; // list of Neuron* Network Neurons to which this Neurons is connected (i.e. receives  input)
int type; // excitatory (0), inhibitory (1)

double g_excite0; // excitatory conductance
double g_excite1; // excitatory conductance
double h_excite; // excitatory h value
double S_excite; //S/N value associated with excitatory --> type
double f; // strength of external poisson spikes
double fback; // strength of external background spikes
bool Espike; // keeps track of whether neuron spike at given  time step to be used other
};

Neuron* Network = new Neuron[N];
double oldv;
double gE0, gE1;
double a0,a1,b0,b1,k1,k2,vtilda,k1tilda,k2tilda,tspike;
double vT = 1.0; // threshold voltage
double vR = 0.0; // reset voltage
double Tfinreal = 3000.0; // should not be less than 1500, to see the entire odor on period
double Tfin = Tfinreal/20.0; // final time
double delt = Tfin/T; // step size
double vE = (14.0/3.0);
double gL = 1.0;
double C = 1.0;
double overC = 1.0/C;
double sigma_Ereal = 1.0;
double sigma_E = sigma_Ereal/20.0; // rise & decay raste of  conductance pulse (excite)
double oversigmaE = 1.0/sigma_E;
double alpha_sigmaE = exp(-delt/sigma_E);
// background input data
double nubackreal = 2700.0;
double nuback = nubackreal*0.02; // Poisson rate for background input
double overNnuback = 1.0/(N*nuback);
double feback = 0.0039; // strength of background input to PNs
int background;// 0: odor spikes, 1: background spikes

//coupling strengths between neurons of different types
double S_ee = 0.27225; // from excitatory to excitatory
// keeps track of S values; scaled by NE & NI
double SFrom[][2] = { {S_ee/(0.75*N), 0},{0, 0} };
/*SFrom[0][0] = S ee/N;
SFrom[1][1] = S ii/N; SFrom[0][1] = S ei/N; SFrom[1][0] = S ie/N);*/
double p_ee = 0.1; // probability of E-->E connection
double pfrom[][2]= { {p_ee, 0}, {0, 0} };

int E = 0;
int EorI;
int m_e =0;
int M_e = 0;
double lfp_e = 0;
//create all files for reading out data

ofstream M("m.csv");
ofstream V("voltages.csv");
ofstream L_E("lfp_e.csv");
ofstream GE("Excitatoryconductances.csv");

// function prototypes
void create_neurons (Neuron* Network, int N,int T);
void connect_network (Neuron* Network, int N);
void connect_network_sparse (Neuron* Network, int N);
void generalupdate_gh(Neuron &n, double alpha_sigmaE, double delt, int EorI, double oversigmaE); // general update for all g, h values for all neurons
void vspike(Neuron* Network, double tspike, int k); // put internal spiketimes in other neuron's spike lists
void spikeupdate_gh(Neuron &n, int EorI, double time, double oversigmaE); // update any neurons that  received spikes during this time step
void spikeupdate_poisson(Neuron* Network, int N, double t, double delt, double oversigmaE, bool background);
void raster_plot(Neuron* Network, int N); //create data output for  raster plot (to be sent to matlab)
// functions
void create_neurons (Neuron* Network, int N, int T) {
  double volt;
  for (int i = 0; i < N; i++) {
    Neuron n; // all voltages start at 0
    n.v = 0.0;
    // each neuron can start with a random initial voltage in [0,1] // volt = (double)rand()/RAND MAX;
    // n.v = volt;
           // exactly 25% inhibitory neurons (same 25%)

    n.type = 0; // for an all-excitatory network (comment out if-else above) // n.type = 0;
    // initialize conductance vectors g(t) = h(t) = 0 since  there were previously no spikes
    n.g_excite0 = 0.0;
    n.g_excite1 = 0.0;
    n.h_excite = 0.0;
           // determine S values based on neuron's type
    n.S_excite = SFrom[0][n.type];
    n.Espike = false;
    // added terms
    n.tsp_excite = new double[1000];
     n.tsp_Ecount = 0;
    // record strength of poisson spikes
    // set 1/3 of the LNs and PNs to have non-zero odor input strengths
    n.f = 0;
       // set all PNs to have the same background input strengths
    n.fback = feback;
    //Uniformly Forces All E-cells and I cells (comment out  if-else above)
    /*if (n.type == 1){
    n.f = fi; n.fback =fiback;
    }
    else{
    n.f = fe; n.fback =feback;
    }*/
    Network[i] = n; // add this neuron to the network
  }
  //cout << "There are " << INHIB << "/" << N << " (" <<  (double)100*INHIB/N << "%) inhibitory neurons in the  network." << endl;
}
void connect_network (Neuron* Network, int N) { // all-to-all connection
  // all to all network
  for (int i = 0; i < N; i++) { // go through all neurons
    for (int j = 0; j < i; j++) {// for all neurons w/ index  less than i, connect to neuron i
      Network[i].connections.push_back(j);
    }
    for (int j = i+1; j < N; j++) { // for all neurons w/ index  more than i, connect to neuron i
      Network[i].connections.push_back(j);
    }
  }
}
void connect_network_sparse (Neuron* Network, int N) { // sparse  network connection
  ofstream C("connections.csv");
  double r = rand();
  r = r/RAND_MAX;
  int rowtype, coltype;
  for (int i = 0; i < N; i++) {// go through all neurons
    rowtype = Network[i].type;
    //C << i;
    for (int j = 0; j < i; j++) { // for all neurons w/ index  less than i, connect to neuron i
      coltype = Network[j].type;
      if (r <= pfrom[rowtype][coltype]) {
        Network[i].connections.push_back(j);
        C << "," << 1;
      }
      else {
        C << "," << 0;
      }
      r = rand();
      r = r/RAND_MAX;
    }
    C << "," << 0;
    for (int j = i+1; j < N; j++) { // for all neurons w/ index  more than i, connect to neuron i
      coltype = Network[j].type;
      if (r <= pfrom[rowtype][coltype]) {
        Network[i].connections.push_back(j);
        C << "," << 1;
      }
      else {
        C << "," << 0;
      }
      r = rand();
      r = r/RAND_MAX;
    }
    C << "," << endl;
  }
  C.close();
}
void generalupdate_gh(Neuron &n, double alpha_sigmaE, double  delt, int EorI, double oversigmaE){
  double G, H;
  if (EorI == 0) {
    G = n.g_excite0*alpha_sigmaE + n.h_excite*alpha_sigmaE*delt*oversigmaE;
    H = n.h_excite*alpha_sigmaE;
    n.g_excite1 = G;
    n.h_excite = H;
  }
}
void vspike(Neuron* Network, double tspike, int k) {
  int foo, foo2;
  if (Network[k].type == 0) {//(excite)
    for (int i = 0; i < Network[k].connections.size(); i++) {//  for all of the connected neurons i connected to neuron k
    foo = Network[k].connections[i]; // neuron number of  connection i
    foo2 = Network[foo].tsp_Ecount;
    Network[foo].tsp_excite[foo2] = tspike; // mark this time down as the foo2-nd spiketime
    Network[foo].tsp_Ecount = Network[foo].tsp_Ecount+1; // increment tsp Ecount (i.e. neuron foo received an additional spike)
    if (Network[foo].tsp_Ecount >= 999) cout << "We are over the allocated length of tsp_Ecount" << endl; Network[foo].Espike = true;
    }
  }
}
void spikeupdate_gh(Neuron &n, int EorI, double time, double  oversigmaE){
  if (EorI == 0) {
    // we will simply add onto the *next* g and h values (at  tn+1) to account for the spikes that occured before tn+1
    double Esumg = 0; double Esumh = 0;
    for (int j = 0; j< n.tsp_Ecount; j++) {
      Esumh = Esumh+ exp(-(time - n.tsp_excite[j])*oversigmaE);
      Esumg = Esumg+ exp(-(time - n.tsp_excite[j])*oversigmaE)*(time-n.tsp_excite[j])*  oversigmaE;
    }
    if (n.type == 0) { // excitatory
      n.g_excite1 = n.g_excite1 + Esumg*n.S_excite*oversigmaE;
      n.h_excite = n.h_excite + Esumh*n.S_excite*oversigmaE;
    }
  }
}

void spikeupdate_poisson(Neuron* Network, int N, double t, double  delt, double oversigmaE, double overNnu, int background){
  double tsp = t;
  double r, u, f;
  int n;
  while (tsp <= t+delt) {
    r = rand();
    u = r/RAND_MAX;
    while (u == 0) {
      r = rand();
      u = r/RAND_MAX;
    }
      tsp=-log(u)*overNnu+tsp;
    if (tsp <= t+delt) {
      n = rand()%N; // which neuron will receive the spike
      f = Network[n].fback;
      Network[n].g_excite1 = Network[n].g_excite1 + f*exp(-(t  - tsp)*oversigmaE)*(t -tsp)*oversigmaE*oversigmaE;
      Network[n].h_excite = Network[n].h_excite + f*exp(-(t -  tsp)*oversigmaE)*oversigmaE;
    }
  }
}
void raster_plot(Neuron* Network, int N) {
  ofstream raster("raster_plot.csv");
  for (int i = 0; i < N; i++) {
    if (Network[i].type == 0) {
      raster << 0 << ",";
      if (Network[i].tspN.size() != 0) {
        for (int j = 0; j < Network[i].tspN.size()-1; j++) {
          raster << Network[i].tspN[j] << ",";
        }
        raster << Network[i].tspN[Network[i].tspN.size()-1]  << endl;
      }
      else {
        raster << " ," << endl;
      }
    }
  }
  raster.close();
}
