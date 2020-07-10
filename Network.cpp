//
//  Network.cpp
//  Alpha_Coupling_Conductance_Based
//
//  Created by Yiyin Hu on 2020/7/9.
//  Copyright © 2020 Yiyin Hu. All rights reserved.
//

/* Network.cpp */
#include "Network.h"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <ctime>
#include <new>
using namespace std;

int main(){
    srand((int)time(NULL)); // random seed
// compile list of select neurons for which we will keep data
    int jump = N/3; // to keep all neuron data, jump = 1;
    if (jump%4 == 0) jump++; // avoids choosing all LNs
    vector <int> dataneurons;
    for (int j = 0; j < N; j = j+jump) // add select neurons to vector
    {
        dataneurons.push_back(j);
    }
    create_neurons(Network, N, T); // creates new neurons  and initializes them
    connect_network(Network, N); // connects neurons, all-to-all
    //connect network sparse(Network, N); // connects neurons, sparse
    V << t << ",";
    GE << t << ",";
    for (int a = 0; a < dataneurons.size()-1; a++) {
        V << Network[dataneurons[a]].v << ",";
        GE << Network[dataneurons[a]].g_excite1 << ",";
    }
    V << Network[dataneurons.size()-1].v << endl;
    GE << Network[dataneurons.size()-1].g_excite1 << endl;
    // local field potential
    lfp_e = 0;
    for (int c = 0; c < N; c++) {
        if (Network[c].type == 0) // excite
        {
            lfp_e = lfp_e + Network[c].v;
        }
    }
    L_E << t << "," << lfp_e/N<< endl;
    
    // spike rate
    M << t << "," << m_e << endl;
    for (int i = 0; i < T; i++) // i increments with time step
    {
       // read out the step and time every 20000 time steps
        if (i%20000 == 0) {
            cout << "i = " << i << endl;
            cout << "time: " << t << endl;
        }
        // calculates t(n+1)'s g and h values, as if no spikes were  received during that time step
        for (int j = 0; j < N; j++) {
            generalupdate_gh(Network[j], alpha_sigmaE, delt, E,oversigmaE);
        }
        // calculate variables: a0,a1,b0,b1,k1,k2 for RK2 method
        for (int k = 0; k < N; k++) {
            oldv = Network[k].v;
            gE0 = Network[k].g_excite0;
            
            a0 = (gL+gE0)*overC;
            b0 = (gL*vR+gE0*vE)*overC;
            
            gE1 = Network[k].g_excite1;
            a1 = (gL+gE1)*overC;
            b1 = (gL*vR+gE1*vE)*overC;
            k1 = -a0*oldv+b0;
            k2 = -a1*(oldv+delt*k1)+b1;
            
            Network[k].v = oldv+(delt/2)*(k1+k2);
            if (Network[k].v > vT) // Neuron has spiked. Reset  potential and integrate to end of time step.
            {
                tspike = t+delt*((vT-oldv)/(Network[k].v-oldv));
                vtilda = (vR - 0.5*(tspike-t)*(b0+b1-delt*a1*b0))/ (1.0-0.5*(tspike-t)*(a0+a1-delt*a0*a1));
                k1tilda = -a0*vtilda+b0;
                k2tilda = -a1*(vtilda+k1tilda*delt)+b1;
                Network[k].v = vtilda+(delt/2)*(k1tilda+k2tilda);// replace with new spiked voltage value
                Network[k].tspN.push_back(tspike);
                vspike(Network, tspike, k); // labels tspike time as a "current spike" for neuron's connections
                // update count for firing rates
                if (Network[k].type == 0) // excite
                {
                    m_e++;
                    M_e++;
                }
                if (Network[k].v > vT) {
                    cout << "Two spikes occurred in one time step. Vtilda was above threshold. Time step may be too large."<<endl;
                    raster_plot(Network, N);
                    return 0;
                }
            }
        }
        spikeupdate_poisson(Network, N, t, delt, oversigmaE, overNnuback, 1); // background spikes
        for (int l = 0; l < N; l++) { // check which neurons have received spikes, if so,  update their g and h values accordingly
            if (Network[l].Espike) {
                spikeupdate_gh(Network[l], E, t+delt, oversigmaE);
                Network[l].tsp_Ecount = 0;
                Network[l].Espike = false;
            }
        }
        //cout<<"Here";
        // collect data from this time step
        t = t+delt;
        
        // voltage and conductances (these .v values are actually  v[i+1] values )
        V << t << ",";
        GE << t << ",";
        
        for (int a = 0; a < dataneurons.size()-1; a++) {
            V << Network[dataneurons[a]].v << ",";
            GE << Network[dataneurons[a]].g_excite1 << ",";
        }
        V << Network[dataneurons.size()-1].v << endl;
        GE << Network[dataneurons.size()-1].g_excite1 << endl;
        // shift all current conductance values to be old g  allowing new g to be calculated at the next time step
        for (int b = 0; b < N; b++) {
            Network[b].g_excite0 = Network[b].g_excite1;
        }
        // local field potential
        lfp_e = 0;
        for (int c = 0; c < N; c++) {
            if (Network[c].type == 0) // excite
            {
                lfp_e = lfp_e + Network[c].v;
            }
        }
        L_E << t << "," << lfp_e/N << endl;
        // spike rate
        M << t << "," << m_e<<endl;
        m_e = 0;
    }
    raster_plot(Network, N); //collect data for rastor plots
    
    V.close();
    GE.close();
    L_E.close();
    M.close();
    for (int n =1; n < N; n++)
    {
        delete[] Network[n].tsp_excite;
    }
    delete[] Network;
    cout << clock()/CLOCKS_PER_SEC << " seconds " << endl;
}

