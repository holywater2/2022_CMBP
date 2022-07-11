#ifndef ____AA_Metropolis____
#define ____AA_Metropolis____ 

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
#include "../../Ewald_sum.cpp"

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

class AA_Metropolis{

    public:
        const int Lx; // Short range (Time domain)
        const int Ly; // Long range  (Spatial domain)
        const int N;
        const int Bin;
        const double B;
        const double Jx;
        const double Jy;
        const double alpha; // dimension = 1
        bool isTinf;

        double cur_beta;

        const int XNN = 1;
        const int YNN;

        double HH;
        int sigma;
        int staggered;
        clock_t __start__, __finish__;

        ewald_ND e2d;

        vector<double> MV;
        vector<double> CV;
        vector<double> TV;
        vector<double> BetaV;
        vector<double> res;

        short* sc; // Square lattice configuration of 2D Ising model
        // double prob[5];
        bool* sign;

        unsigned long long Fliped_Step = 0;
        unsigned long long Total_Step  = 0;
        long Calc_call = 0;

        static string Name(){return "AA_Metropolis";}
        AA_Metropolis(int Lx, int Ly, int bin, double B, double Jx, double Jy, double alpha, double Tsrt, double Tfin, bool isTinf);
        AA_Metropolis(vector<double> args);
        ~AA_Metropolis(){
            // delete e2d;
            __finish__ = clock();
            cout << "------------------------------------------------------------------------------------------------------------------\n";
            cout << "Model calculation finished. Spent time: " << (double)(__finish__-__start__)/CLOCKS_PER_SEC << "\n";
            cout << "------------------------------------------------------------------------------------------------------------------\n";
        }
        long double Prob(double delta);
        void Initialize(double beta);
        // void Initialzie(int idx);
        void Measure();
        void Measure_fast();
        void Calculate(int _n = 0,bool Random = true);
        void IterateUntilEquilibrium(int equil_time,bool random = true);
};

AA_Metropolis::AA_Metropolis(int Lx, int Ly, int bin, double B, double Jx, double Jy, double alpha, double Tsrt, double Tfin, bool isTinf)
:Lx(Lx), Ly(Ly), N(Lx*Ly), Bin(bin), B(B), Jx(Jx), Jy(Jy), YNN(Lx), alpha(alpha) {
    this-> isTinf = isTinf;
    this-> sc = new short[N];
    this-> sign = new bool[N];
    // Checkboard style
    // if(Lx%2 == 1)
    //     for(int i = 0; i < N; i++)
    //         sign[i] = i%2;
    // else
    //     for(int i = 0; i < N; i++)
    //         sign[i] = (i+(i/Lx))%2;
    // for(int i = 0; i < N; i++){
    //     cout << -sign[i]*2+1;
    // }
    for(int i = 0; i < N; i++)
        sign[i] = (i/Lx)%2;
    
    

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
    e2d = ewald_ND(2,vector<int>({Lx,Ly}),alpha);
}

AA_Metropolis::AA_Metropolis(vector<double> args):
AA_Metropolis(args[0],args[1],args[2],args[3],args[4],args[5],args[6],args[7],args[8],args[9]){}

long double AA_Metropolis::Prob(double delta){
    return expl(-this->cur_beta*delta);
}

void AA_Metropolis::Initialize(double beta){
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

void AA_Metropolis::Measure(){ //O(N^2)
    int i, j;
    int sigma = 0;
    int staggered = 0;
    double HH, result = 0;
    // cout << e2d.pi_ij(0,56) << '\n';
    for(i = 0; i < N; i++){
        double temp = 0;
        for(j = i+Lx; j < N; j+=Lx) // y 방향 long range 계산
            temp += e2d.pi_ij_1D(i,j)*sc[j];
        temp *= Jy;
        if(i%Lx != Lx-1) // x 방향 short ragne 계산
            temp += Jx*sc[i+1];
        else
            temp += Jx*sc[i+1-Lx];
        sigma += sc[i];
        staggered += sc[i]*(-sign[i]*2+1);
        result += temp*sc[i];
    }
    // result*= 0.5;   // 전체의 절반만 계산했음. 해밀토니안에 2를 안곱해도 계수 조정해서 찾을 수 있지 않을까?
                    // i=0, j=i
    HH = -result-B*sigma;
    this->HH = HH;
    this->sigma = sigma;
    this->staggered = staggered;
}

void AA_Metropolis::Measure_fast(){
    // double res[] = {HH,sigma,staggered};
    // return vector<double>(res,res +sizeof(res)/sizeof(res[0]));
    return; // do nothing
}

void AA_Metropolis::Calculate(int _n, bool Random){ //O(N^2)
    int i, k, n;
    n = !_n ? (this->N) : _n;
    for(i = 0; i < n; i++){
        // Sweep Randomly
        if(Random){
            k = (this->N)*dis(gen);
        // Sweep Sequential
        } else if(n%2 == 0){
            k = 2*i;
            if(k < N) k = (int(k/Lx))%2 == 0 ? k+1 : k;
            else k = (int(k/Lx))%2 == 0 ? k-N : k-N+1;
        } else{
            k = 2*i >= N ? 2*i-N: 2*i;
        }

        double delta = 0;
        for(int jj = k%Lx; jj < N; jj+=Lx)  // Long range diff
            delta += Jy*e2d.pi_ij(jj,k)*sc[jj];

        int nn;                             // Short range diff
        if(((nn = k - XNN)+1)%Lx == 0) nn += this->Lx; 
        delta += Jx*sc[nn];
        if(((nn = k + XNN)-1)%Lx == Lx-1) nn -= this->Lx;
        delta += Jx*sc[nn];
        
        delta *= 2*sc[k];
        this->Total_Step++;

        if((delta <= 0) || (dis(gen) < Prob(delta))){
            this->Fliped_Step++;
            this->sc[k] *= -1;
            this->sigma += 2*(this->sc[k]);
            this->staggered +=2*(this->sc[k]*(-sign[k]*2+1));
            this->HH += delta;
        }
    }
}

void AA_Metropolis::IterateUntilEquilibrium(int equil_time,bool random){
    for(int j = 0; j < equil_time; j++)
        Calculate(0,random);
}

#endif // ____AA_Metropolis____ 