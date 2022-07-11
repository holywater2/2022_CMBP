#ifndef ____SwendsenWang____
#define ____SwendsenWang____ 

// File & IO System
// #include <iostream>
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

using namespace std;

// static random_device rd;  // Will be used to obtain a seed for the random number engine
// static mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
// static uniform_real_distribution<> dis(0.0, 1.0);

static long unsigned int seed = static_cast<long unsigned int>(time(0)); 
static mt19937 gen(seed); // Standard mersenne_twister_engine seeded with time()
static uniform_real_distribution<> dis(0.0, 1.0);

//In 2D Ising model => 2/log(1+sqrt(2))
const double T_CRIT = 2.269185; 

typedef tuple<int,int> duo;

class SwendsenWang_2D{ 
    public:
        const int L;
        const int N;
        const int Bin;
        const int B;
        const int J;
        bool isTinf;

        const static int init = 1.01;
        const int XNN = 1;
        const int YNN; 

        clock_t __start__, __finish__;

        // #ifdef _WIN32
        // static string Filename = ".\\Result\\Metropolis_c_"+to_string(L)+"_int"+to_string(bin);
        // #endif
        // #ifdef linux
        // static string Filename = "./Result/Metropolis_c_"+to_string(L)+"_int"+to_string(bin);
        // #endif

        vector<double> MV;
        vector<double> CV;
        vector<double> TV;
        vector<double> BetaV;
        vector<double> res;

        short* sc; // Square lattice configuration of 2D Ising model
        double prob;

        unsigned long Fliped_Step = 0;
        unsigned long Total_Step  = 0;
        unsigned long Calc_call = 0;

        static string Name(){return "SwendsenWang";}
        SwendsenWang_2D(int L, int bin, double B, double J, double Tsrt, double Tfin, bool isTinf);
        SwendsenWang_2D(vector<double> args);
        ~SwendsenWang_2D(){
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
        int Calculate();
        int Calculate(int site);
        void IterateUntilEquilibrium(int equil_time);
};

SwendsenWang_2D::SwendsenWang_2D(int L, int bin, double B, double J, double Tsrt, double Tfin, bool isTinf) :L(L), N(L*L), Bin(bin), B(B), J(J), YNN(L){
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
}

SwendsenWang_2D::SwendsenWang_2D(vector<double> args): SwendsenWang_2D(args[0],args[1],args[2],args[3],args[4],args[5],args[6]){}

void SwendsenWang_2D::ProbCalc(double beta){
    this->prob = 1-exp(-2*beta*J);
}

void SwendsenWang_2D::Initialize(double beta){
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

}

int SwendsenWang_2D::SweepHelical(int i){
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

int SwendsenWang_2D::BoundaryHelical(int i){
    int nn, sum = 0;
    // int XNN = 1, YNN = L;
    
    if((nn = i + XNN) == N) nn = 0;
    sum += this->sc[nn];
    if((nn = i + YNN) >= N) nn -= N;
    sum += this->sc[nn];

    return sum;
}

duo SwendsenWang_2D::Measure(){
    int i, sum, HH;
    int res = 0, sigma = 0;
    for(i = 0; i < N; i++){
        sum = this->BoundaryHelical(i);
        res += J*sum*sc[i];
        sigma += sc[i];
    }
    
    HH = -res-B*sigma;

    return make_tuple(HH,sigma);
}

// void SwendsenWang_2D::Calculate(int _n){
//     int i, k, delta, n;
//     n = !_n ? (this->N) : _n;
//     double a;
//     for(i = 0; i < n; i++){
//         // Sweep Randomly
//         k = (this->N)*dis(gen);
//         // Sweep Sequential
//         // k = n

//         delta = (this->SweepHelical(k))*(this->sc[k]);
        
//         a = dis(gen);
//         this->Total_Step++;

//         if(delta <= 0){
//             this->Fliped_Step++;
//             this->sc[k] = -(this->sc[k]);
//         } else if(a < prob[delta]){
//             this->Fliped_Step++;
//             this->sc[k] = -(this->sc[k]);
//         }
//     }
// }

int SwendsenWang_2D::Calculate(){
    /*In SwendsenWang algorithm, each clustering and flipping process is one mcs*/
    int i, k, l, sp, cur, size = 0;
    short snew, sold;
    vector<int> stack(N);
    vector<int> clusterd(N); // -2, +2 Flip, -1, +1 Not Flip, 0 Not clusterd
    vector<int> clusterd_check(N); // -2, +2 Flip, -1, +1 Not Flip, 0 Not clusterd
    vector<short> sc_copy;
    sc_copy.insert(sc_copy.end(),&sc[0],&sc[N]);
    vector<int> nn(4);

    /*k: Choosing random initial site of SwendsenWang step */
    // k = N*dis(gen);

    int cnt = 0;

    for(int k = 0; k < N; k++){
        // k += 2;
        // if(k >= N) k -= N;
        // cout << k << endl;
        if(clusterd[k] != 0) continue;
        cnt = cnt+1;

        stack[0] = k;
        sp = 1; // stack pointer
        
        sold = sc_copy[k];
        snew = -sc_copy[k];
        
        sc_copy[k] = snew;

        // Determine whether flip or not
        int Flip_factor = 1; // Not Flip
        double temp = dis(gen);

        if(temp<0.5){
            Flip_factor = 2; // Flip
            Fliped_Step++;
            size++;
        }

        clusterd[k] = Flip_factor*snew;
        clusterd_check[k] = cnt;


        while(sp){
            cur = stack[--sp];
            nn[0] = (cur+XNN+N)%N;
            nn[1] = (cur+YNN+N)%N;
            nn[2] = (cur-XNN+N)%N;
            nn[3] = (cur-YNN+N)%N;

            for(l = 0; l < 4; l++){
                if(clusterd[nn[l]] == 0 && sc_copy[nn[l]] == sold){
                    if(dis(gen) < prob){
                        clusterd[nn[l]] = Flip_factor*snew;
                        clusterd_check[nn[l]] = cnt;

                        stack[sp++] = nn[l];
                        sc_copy[nn[l]] = snew;
                        
                        if(Flip_factor == 2){
                            Fliped_Step++;
                            size++;
                        }
                    }
                }
            }
        }
    }

    // string res = "";
    // for(int i = 0; i < L; i++){
    //     for(int j = 0; j < L; j++){
    //         if(clusterd[i*L+j] == 2) res += "��";
    //         else if(clusterd[i*L+j] == 1) res += "��";
    //         else if(clusterd[i*L+j] == -2) res += "��";
    //         else if(clusterd[i*L+j] == -1) res += "��";	
    //     }
    //     res += "\n";
    // }
    // cout << res << "\n";

    // res = "";
    // for(int i = 0; i < L; i++){
    //     for(int j = 0; j < L; j++){
    //         if(sc[i*L+j] == 1) res += "��";
    //         else if(sc[i*L+j] == -1) res += "��";
    //     }
    //     res += "\n";
    // }
    // cout << res << "\n";

    // res = "";
    // for(int i = 0; i < L; i++){
    //     for(int j = 0; j < L; j++){
    //         res += to_string(clusterd_check[i*L+j]) + " ";
    //     }
    //     res += "\n";
    // }
    // cout << res << "\n";

    for(int k = 0; k < N; k++){
        this->sc[k] += (clusterd[k]/2)*2;
    }

    // res = "";
    // for(int i = 0; i < L; i++){
    //     for(int j = 0; j < L; j++){
    //         if(sc[i*L+j] == 1) res += "��";
    //         else if(sc[i*L+j] == -1) res += "��";
    //     }
    //     res += "\n";
    // }
    // cout << res << "\n";

    Total_Step++;
    return size;
}

int SwendsenWang_2D::Calculate(int site){
    /*In SwendsenWang algorithm, each clustering and flipping process is one mcs*/
    int i, k, l, sp, cur, size = 0;
    short snew, sold;
    vector<int> stack(N);
    vector<int> nn(4);

    /*k: Choosing random initial site of SwendsenWang step */
    // k = N*dis(gen);
    k = site;
    
    stack[0] = k;
    sp = 1; // stack pointer
    
    sold = sc[k];
    snew = -sc[k];
    
    sc[k] = snew;

    while(sp){
        cur = stack[--sp];
        if((nn[0] = cur - XNN) <  0) nn[0] += N;
        if((nn[1] = cur + XNN) >= N) nn[1] -= N;
        if((nn[2] = cur - YNN) <  0) nn[2] += N;
        if((nn[3] = cur + YNN) >= N) nn[3] -= N;
        
        // if((nn[0] = cur - XNN) <  0) nn[0] += N;
        if((nn[0] = cur + XNN) >= N) nn[0] -= N;
        // if((nn[2] = cur - YNN) <  0) nn[2] += N;
        if((nn[1] = cur + YNN) >= N) nn[1] -= N;


        for(l = 0; l < 2; l++){
            if(sc[nn[l]] == sold){
                if(dis(gen) < prob){
                    stack[sp++] = nn[l];
                    sc[nn[l]] = snew;
                    Fliped_Step++;
                    size++;
                }
            }
        }
    }
    
    Total_Step++;
    return size;
}


void SwendsenWang_2D::IterateUntilEquilibrium(int equil_time){
    double k = 0;
    while(k < equil_time){
        k += Calculate()/(double)N;
    }
}

#endif // ____SwendsenWang____ 