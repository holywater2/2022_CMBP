#ifndef ____METROPOLIS_LR____
#define ____METROPOLIS_LR____ 

// File & IO System
#include <iostream>

// Data Structure
#include <vector>
#include <tuple>
#include <string>
#include <tuple>
// Mathmatics
#include <math.h>
#include <random>
// Etc.
#include <stdlib.h>
#include <ctime>
// Ewald sum
#include "../Ewald_sum.cpp"

using namespace std;

// static random_device rd;  // Will be used to obtain a seed for the random number engine
// static mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
// static uniform_real_distribution<> dis(0.0, 1.0);

static long unsigned int seed = static_cast<long unsigned int>(time(0)); 
static mt19937 gen(seed); // Standard mersenne_twister_engine seeded with time()
static uniform_real_distribution<> dis(0.0, 1.0);

//In 2D Ising model => 2/log(1+sqrt(2))
// const double T_CRIT = 2.269185;

typedef tuple<double,int> duo;

class Metropolis_LR_2D{

    public:
        const int L;
        const int N;
        const int Bin;
        const double B;
        const double J;
        const double alpha;
        bool isTinf;

        double cur_beta;

        const int XNN = 1;
        const int YNN;

        double HH;
        int sigma;
        clock_t __start__, __finish__;

        ewald_ND e2d;

        vector<double> MV;
        vector<double> CV;
        vector<double> TV;
        vector<double> BetaV;
        vector<double> res;

        short* sc; // Square lattice configuration of 2D Ising model
        // double prob[5];

        long Fliped_Step = 0;
        long Total_Step  = 0;
        long Calc_call = 0;

        static string Name(){return "Metropolis_LR";}
        Metropolis_LR_2D(int L, int bin, double B, double J, double alpha, double Tsrt, double Tfin, bool isTinf);
        Metropolis_LR_2D(vector<double> args);
        ~Metropolis_LR_2D(){
            // delete e2d;
            __finish__ = clock();
            cout << "------------------------------------------------------------------------------------------------------------------\n";
            cout << "Model calculation finished. Spent time: " << (double)(__finish__-__start__)/CLOCKS_PER_SEC << "\n";
            cout << "------------------------------------------------------------------------------------------------------------------\n";
        }
        long double Prob(double delta);
        void Initialize(double beta);
        // void Initialzie(int idx);
        int SweepHelical(int i);
        int BoundaryHelical(int i);
        duo Measure();
        duo Measure_fast();
        void Calculate(int _n = 0,bool Random = true);
        void IterateUntilEquilibrium(int equil_time,bool random = true);
};

Metropolis_LR_2D::Metropolis_LR_2D(int L, int bin, double B, double J, double alpha, double Tsrt, double Tfin, bool isTinf)
:L(L), N(L*L), Bin(bin), B(B), J(J), YNN(L), alpha(alpha) {
    this-> isTinf = isTinf;
    this-> sc = new short[N];

    this-> MV = vector<double>(Bin);
    this-> CV = vector<double>(Bin);
    this-> TV = vector<double>(Bin);
    this-> BetaV = vector<double>(Bin);

    __start__ = clock();

    for(int i = 0; i < bin; i++){ 
        if(!Tsrt){
            this->TV[i] = Tsrt + ((Tfin-Tsrt)/(double)(bin))*(i+1);
        } else if(bin == 1){
            this->TV[i] = Tsrt;
        } else {
            this->TV[i] = Tsrt + ((Tfin-Tsrt)/(double)(bin-1))*(i);
        }
        this->BetaV[i] = 1/TV[i];
    }
    // Ewald sum
    e2d = ewald_ND(2,vector<int>({L,L}),alpha);
}

Metropolis_LR_2D::Metropolis_LR_2D(vector<double> args):
Metropolis_LR_2D(args[0],args[1],args[2],args[3],args[4],args[5],args[6],args[7]){}

long double Metropolis_LR_2D::Prob(double delta){
    return expl(-this->cur_beta*delta);
}

void Metropolis_LR_2D::Initialize(double beta){
    this-> Calc_call = 0;
    this-> Fliped_Step = 0;
    this-> Total_Step = 0;

    for(int i = 0; i < N; i++){
        // T = 0 start
        sc[i] = 1;

        // T = \inf start
        if(this->isTinf) this->sc[i] -= int(dis(gen)*2)*2;
    }
    // this->ProbCalc(beta);
    cur_beta = beta;

    this->Measure();
}

int Metropolis_LR_2D::SweepHelical(int i){
    int nn, sum = 0;
    // int XNN = 1, YNN = L;

    if((nn = i - XNN) < 0) nn += this->N;
    sum += this->sc[nn];
    if((nn = i + XNN) >= this->N) nn -= this->N;
    sum += this->sc[nn];
    if((nn = i - YNN) < 0) nn += this->N;
    sum += this->sc[nn];
    if((nn = i + YNN) >= this->N) nn -= this->N;
    sum += this->sc[nn];

    return sum;
}

int Metropolis_LR_2D::BoundaryHelical(int i){
    int nn, sum = 0;
    // int XNN = 1, YNN = L;
    
    if((nn = i + XNN) == N) nn = 0;
    sum += this->sc[nn];
    if((nn = i + YNN) >= N) nn -= N;
    sum += this->sc[nn];

    return sum;
}

duo Metropolis_LR_2D::Measure(){ //O(N^2)
    int i, j;
    int sigma = 0;
    double HH, result = 0;
    // cout << e2d.pi_ij(0,56) << '\n';
    for(i = 0; i < N; i++){
        double temp = 0;
        for(j = i+1; j < N; j++){    
            result += J*e2d.pi_ij(i,j)*sc[i]*sc[j];
        }
        sigma += sc[i];
    }
    // result*= 0.5;   // 전체의 절반만 계산했음. 해밀토니안에 2를 안곱해도 계수 조정해서 찾을 수 있지 않을까?
                    // i=0, j=i
    HH = -result-B*sigma;
    this->HH = HH;
    this->sigma = sigma;

    return make_tuple(HH,sigma);
}

duo Metropolis_LR_2D::Measure_fast(){
    return make_tuple(this->HH,this->sigma);
}

void Metropolis_LR_2D::Calculate(int _n, bool Random){ //O(N^2)
    int i, k, n;
    n = !_n ? (this->N) : _n;
    for(i = 0; i < n; i++){
        // Sweep Randomly
        if(Random){
            k = (this->N)*dis(gen);
        // Sweep Sequential
        } else if(n%2 == 0){
            k = 2*i;
            if(k < N) k = (int(k/L))%2 == 0 ? k+1 : k;
            else k = (int(k/L))%2 == 0 ? k-N : k-N+1;
        } else{
            k = 2*i >= N ? 2*i-N: 2*i;
        }

        double delta = 0;
        for(int jj = 0; jj < N; jj++){ // O(N)
            delta += e2d.pi_ij(jj,k)*sc[jj];
            // cout << jj << " " << k << " " << e2d.pi_ij(jj,k) << " " << delta << ' ';
        }
        delta *= 2*sc[k];
        // cout << "res" << k << ' ' << delta << '\n';


        this->Total_Step++;

        // double pb = dis(gen);
        // double target = Prob(delta);

        if((delta <= 0) || (dis(gen) < Prob(delta))){
            this->Fliped_Step++;
            this->sc[k] *= -1;
            this->sigma += 2*(this->sc[k]);
            this->HH += delta;
        }
    }
}

void Metropolis_LR_2D::IterateUntilEquilibrium(int equil_time,bool random){
    for(int j = 0; j < equil_time; j++)
        Calculate(0,random);
}

#endif // ____METROPOLIS_LR____ 