#ifndef ____Anti_Anisotropy____
#define ____Anti_Anisotropy____ 

// File & IO System
#include <iostream>

// Data Structure
#include <vector>
#include <tuple>
#include <string>
#include <tuple>
// Mathmatics
#include <math.h>
#include <cmath>
#include <random>
// Etc.
#include <stdlib.h>
#include <ctime>

using namespace std;

// static random_device rd;  // Will be used to obtain a seed for the random number engine
// static mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
// static uniform_real_distribution<> dis(0.0, 1.0);

static long unsigned int seed = static_cast<long unsigned int>(time(0)); 
static mt19937 gen(seed); // Standard mersenne_twister_engine seeded with time()
static uniform_real_distribution<> dis(0.0, 1.0);

//In 2D Ising model => 2/log(1+sqrt(2))
// const double T_CRIT = 2.269185;

typedef tuple<int,int> duo;

class Anti_Anisotropy_2D{
    public:
        const int L;
        const int N;
        const int Bin;
        const int B;
        const double Jx;
        const double Jy;
        const double Alpha;

        bool isTinf;

        const int XNN = 1;
        const int YNN;

        int HH;
        int sigma;
        clock_t __start__, __finish__;


        // #ifdef _WIN32
        // static string Filename = ".\\Result\\Anti_Anisotorpy_c_"+to_string(L)+"_int"+to_string(bin);
        // #endif
        // #ifdef linux
        // static string Filename = "./Result/Anti_Anisotorpy_c_"+to_string(L)+"_int"+to_string(bin);
        // #endif

        vector<double> MV;
        vector<double> CV;
        vector<double> TV;
        vector<double> BetaV;
        vector<double> res;
        vector<double> Jij;

        short* sc; // Square lattice configuration of 2D Ising model
        double prob[5];

        long Fliped_Step = 0;
        long Total_Step  = 0;
        long Calc_call = 0;

        static string Name(){return "Anti_Anisotorpy";}
        Anti_Anisotropy_2D(int L, int bin, int B, double Jx, double Jy, ,double alpha, double Tsrt, double Tfin, bool isTinf);
        Anti_Anisotropy_2D(vector<double> args);
        ~Anti_Anisotropy_2D(){
            __finish__ = clock();
            cout << "------------------------------------------------------------------------------------------------------------------\n";
            cout << "Model calculation finished. Spent time: " << (double)(__finish__-__start__)/CLOCKS_PER_SEC << "\n";
            cout << "------------------------------------------------------------------------------------------------------------------\n";
        }
        void ProbCalc(double beta);
        void Initialize(double beta);
        // void Initialzie(int idx);
        int SweepHelical(int i);
        int BoundaryHelical(int i);
        duo Measure();
        duo Measure_fast();
        void Calculate(int _n = 0,bool Random = true);
        void IterateUntilEquilibrium(int equil_time,bool random = true);
};

Anti_Anisotropy_2D::Anti_Anisotropy_2D(int L, int bin, int B, double Jx, double Jy, ,double alpha, double Tsrt, double Tfin, bool isTinf)
:L(L), N(L*L), Bin(bin), B(B), Jx(Jx), Jy(Jy), YNN(L), Alpha(alpha) {
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

    this-> Jij = vector<vector<double>>(N,vector<double>(N,0));
    vector<double> site_cache(L,-1);
    for(int i = 0; i < N; i++){
        for(int j = 0; j < N; j++){
            if(i == j) Jij[i][j] = 0;
            else if(j/L == i/L){
                if(site_cache[abs(j%L-i%L)] == -1)
                    site_cac   he[abs(j%L-i%L)] = pow((abs(j%L-i%L)),-2-alpah);
                Jij[i][j] = Jx/site_cache[abs(j%L-i%L)];
            }
            else Jij[i][j] = Jy;
        }
    }
}

Anti_Anisotropy_2D::Anti_Anisotropy_2D(vector<double> args): Anti_Anisotropy_2D(args[0],args[1],args[2],args[3],args[4],args[5],args[6],args[7],args[8],args[9]){}

void Anti_Anisotropy_2D::ProbCalc(double beta){
    for(int i = 2; i < 5; i += 2){
        this->prob[i] = exp(-2*beta*i);
    }
}

void Anti_Anisotropy_2D::Initialize(double beta){
    this-> Calc_call = 0;
    this-> Fliped_Step = 0;
    this-> Total_Step = 0;

    for(int i = 0; i < N; i++){
        // T = 0 start
        sc[i] = 1;

        // T = \inf start
        if(this->isTinf) this->sc[i] -= int(dis(gen)*2)*2;
    }
    this->ProbCalc(beta);

    this->Measure();
}

int Anti_Anisotropy_2D::SweepHelical(int i){
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

int Anti_Anisotropy_2D::BoundaryHelical(int i){
    int nn, sum = 0;
    // int XNN = 1, YNN = L;
    
    if((nn = i + XNN) == N) nn = 0;
    sum += this->sc[nn];
    if((nn = i + YNN) >= N) nn -= N;
    sum += this->sc[nn];

    return sum;
}

duo Anti_Anisotropy_2D::Measure(){
    int i, sum, HH;
    int res = 0, sigma = 0;
    for(i = 0; i < N; i++){
        sum = this->BoundaryHelical(i);
        res += J*sum*sc[i];
        sigma += sc[i];
    }
    
    HH = -res-B*sigma;
    this->HH = HH;
    this->sigma = sigma;

    return make_tuple(HH,sigma);
}

duo Anti_Anisotropy_2D::Measure_fast(){
    return make_tuple(this->HH,this->sigma);
}


void Anti_Anisotropy_2D::Calculate(int _n, bool Random){
    int i, k, delta, n;
    n = !_n ? (this->N) : _n;
    double a;
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
        
        delta = (this->SweepHelical(k))*(this->sc[k]);
        
        this->Total_Step++;

        if((delta <= 0) || (dis(gen) < prob[delta])){
            this->Fliped_Step++;
            this->sc[k] *= -1;
            this->sigma += 2*(this->sc[k]);
            this->HH += 2*delta;
        }
    }
}

void Anti_Anisotropy_2D::IterateUntilEquilibrium(int equil_time,bool random){
    for(int j = 0; j < equil_time; j++)
        Calculate(0,random);
}

#endif // ____Anti_Anisotropy____ 