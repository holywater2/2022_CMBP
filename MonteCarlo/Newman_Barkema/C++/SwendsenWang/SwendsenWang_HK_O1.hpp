// Unoptimized version : It works now!
#ifndef ____SwendsenWang_HK____
#define ____SwendsenWang_HK____ 

// Macro that will be discharged at the end of header
#define cor(x,y) ((y)*L+x+N)%N              // HPBC
// #define cor(x,y) ((y)*L+((x+L)%L)+N)%N   // PBC

#define rng() dis(gen)

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

using namespace std;

/* Random Number Generator */
// static random_device rd;  // Will be used to obtain a seed for the random number engine - device dependent
// static mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()

static long unsigned int seed = static_cast<long unsigned int>(time(0)); // Set the seed of rng by current time.
static mt19937 gen(seed); // Standard mersenne_twister_engine seeded with time()
static uniform_real_distribution<> dis(0.0, 1.0);

// typedef
typedef tuple<int,int> duo;

// In 2D Ising model => T_CRIT = 2/log(1+sqrt(2))
const double T_CRIT = 2.269185; 

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

        vector<double> MV;
        vector<double> CV;
        vector<double> TV;
        vector<double> BetaV;
        vector<double> res;

        // short* sc;
        vector<short> sc;                      // Square lattice configuration of 2D Ising model
        double prob;

        vector<short> iBondLeft  = vector<short>(N);
        vector<short> jBondAbove = vector<short>(N);

        vector<int> label;
        vector<int> labels;

        unsigned long Fliped_Step = 0;
        unsigned long Total_Step  = 0;
        unsigned long Calc_call = 0;

        // Return name of Model
        static string Name(){return "SwendsenWang_HK";}

        // Constructor of Model
        SwendsenWang_2D(int L, int bin, double B, double J, double Tsrt, double Tfin, bool isTinf);
        SwendsenWang_2D(vector<double> args);
        
        // Destructor of Model
        ~SwendsenWang_2D(){
            __finish__ = clock();
            cout << "------------------------------------------------------------------------------------------------------------------\n";
            cout << "Model calculation finished. Spent time: " << (double)(__finish__-__start__)/CLOCKS_PER_SEC << "\n";
            cout << "------------------------------------------------------------------------------------------------------------------\n";
        }

        // Initialize the Ising model according to beta = 1/kT
        void Initialize(double beta);
        // Initialize the Ising model according to beta = 1/kT
        void ProbCalc(double beta);

        // Boundary condition used for Measurement (calculate right and below interaction)
        int BoundaryHelical(int i);
        // Boundary condition on Sweep (calculate)
        int SweepHelical(int i);

        // Measure sum of magnetizaiton and Hamiltonian of Spin configuration on current situation
        duo  Measure();

        // Sub step of Monte Carlo step.
        int  Calculate();
        void IterateUntilEquilibrium(int equil_time);

        // HK algorithms parts for SwendsenWang
        void MakeBonds();
        void ClusterAndLabel();
        void Union(int x, int y);
        int  Find(int x);
        void FlipCluster();
        string ClusterBoxDrawing();

        vector<string> A = {"·","╴","╶","─",\
                    "╷","┐","┌","┬",\
                    "╵","┘","└","┴",\
                    "│","┤","├","┼"};
};

SwendsenWang_2D::SwendsenWang_2D(int L, int bin, double B, double J, double Tsrt, double Tfin, bool isTinf) :L(L), N(L*L), Bin(bin), B(B), J(J), YNN(L){
    this-> isTinf = isTinf;
    // this-> sc = new short[N];
    this-> sc = vector<short>(N);

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
        if(this->isTinf) this->sc[i] -= int(rng()*2)*2;
    }
    this->ProbCalc(beta);
}

int SwendsenWang_2D::SweepHelical(int i){
    int nn, sum = 0;

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

int SwendsenWang_2D::Calculate(){
    /*In SwendsenWang algorithm, each clustering and flipping process is one mcs*/
    int size = Fliped_Step;
    
    MakeBonds();
    ClusterAndLabel();
    FlipCluster();
    
    Total_Step++;
    return Fliped_Step-size;
}

void SwendsenWang_2D::IterateUntilEquilibrium(int equil_time){
    int k = 0;
    while(k++ < equil_time){
        Calculate();
    }
}

void SwendsenWang_2D::MakeBonds(){
    for(int i = 0; i < N; i++){
        this->iBondLeft[i]     = false;
        this->jBondAbove[i]    = false;
    }

    int nn;
    for(int i = 0; i < N; i++){
        nn = (i-XNN+N)%N;
        if((sc[i] == sc[nn]) && (rng() < prob))
            iBondLeft[i] = true;
        nn = (i-YNN+N)%N;
        if((sc[i] == sc[nn]) && (rng() < prob))
            jBondAbove[i] = true;
    }
}

/* Debugging function*/
string SwendsenWang_2D::ClusterBoxDrawing(){
    string res = "";
    for(int j = 0; j < L; j++){
        for(int i = 0; i < L; i++){
            res += A[iBondLeft[cor(i,j)] + 2*iBondLeft[cor(i+1,j)] + 4*jBondAbove[cor(i,j+1)]+ 8*jBondAbove[cor(i,j)]];
        }
        res += "\n";
    }
    return res;
}

void SwendsenWang_2D::ClusterAndLabel(){
    // Largest label is used checked to make unused labeling number.
    int largest_label = 0;
    // Array that contains the labeling number of each site.
    label = vector<int>(N);
    // Array that contains the matching number of label (ex) label[4] = 2 means #4 cluster is same label with #2 cluster.
    // labels need to be initialized to label[x] = x which imply each cluster group is initialiy set independent.
    labels = vector<int>(N+1);
    for(int cnt = 0; cnt < N+1; cnt++) labels[cnt] = cnt;

    // bond[] contains the site of neareast bond.
    int bond[4];
    // bonds_num is the number of bonding from current site.
    int bonds_num;

    string box = ClusterBoxDrawing(); // For debugging purpose.

    for(int j = 0; j < L; j++){
        for(int i = 0; i < L; i++){
            bonds_num = 0;
            
            // HPBC(Helical periodic boundary condition) is apllied -> so right most site is conneted the right left most site of the next line.
            if(!(i == 0 && j == 0)      && iBondLeft [cor(i,j)])    bond[bonds_num++] = cor(i-1,j); // Left should be checked everywhere except 0,0
            if((j > 0)                  && jBondAbove[cor(i,j)])    bond[bonds_num++] = cor(i,j-1); // Above should be checked except 0th row
            if((j == L-1)               && jBondAbove[cor(i,j+1)])  bond[bonds_num++] = cor(i,j+1); // Below should be checked at last row
            if((i == L-1 && j == L-1)   && iBondLeft [0])           bond[bonds_num++] = cor(i+1,j); // Right shoudl be checked at last site

            // // PBC(periodic boundary condition)
            // if((i > 0)                  && iBondLeft [cor(i,j)])    bond[bonds_num++] = cor(i-1,j); // Left should be checked except 0th column
            // if((j > 0)                  && jBondAbove[cor(i,j)])    bond[bonds_num++] = cor(i,j-1); // Above should be checked except 0th row
            // if((j == L-1)               && jBondAbove[cor(i,j+1)])  bond[bonds_num++] = cor(i,j+1); // Below should be checked at last row
            // if((i == L-1)               && iBondLeft [cor(i+1,j)])  bond[bonds_num++] = cor(i+1,j); // Right shoudl be checked at last column


            // If there is no bond at current site, set it as new cluster number.
            if(bonds_num == 0){
                label[cor(i,j)] = ++largest_label;
            } else {
                // If there are more than one bonds, then find the minimum number of label and union all label of near bonds to that minimum label.
                // And set current site as the label of the minimum bond.
                
                // Step of finding minimum (proper) label among bonds
                int min_label = largest_label;
                for(int ll = 0; ll < bonds_num; ll++){
                    int proper_label = Find(label[bond[ll]]);
                    min_label = min_label > proper_label ? proper_label : min_label;
                }
                
                // Step of union near label to minimum label
                for(int ll = 0; ll < bonds_num; ll++){
                    Union(label[bond[ll]],min_label);
                }
                // Common part: Set the current site label
                label[cor(i,j)] = Find(min_label);
            }
        }
    }
}

void SwendsenWang_2D::FlipCluster(){
    bool sNewChosen[N] = {};
    short sNew[N];
    for(int i = 0; i < N; i++){
        label[i] = Find(label[i]);
    }

    for(int i = 0; i < N; i++){
        if (!sNewChosen[label[i]]) {    
            sNew[label[i]] = 2*int(2*rng())-1;
            sNewChosen[label[i]] = true;
        }
        if(this->sc[i] != sNew[label[i]]){
            this->sc[i] = sNew[label[i]];
            Fliped_Step++;
        }
    }
}

void SwendsenWang_2D::Union(int x, int y){
    // Maybe y should be smaller I guess : labels[x] = y should satisfy x >= y.
    int xx = Find(x);
    int yy = Find(y);
    if (xx > yy)
      labels[xx] = yy;
    else
      labels[yy] = xx;

    // labels[Find(x)] = Find(y);
}

// This function has two roles. First is find the minimum matching label of the number "x".
// Second it make direct data set by substitue current matching label to number "x".
int SwendsenWang_2D::Find(int x){
    int y = x;
    // First : Find minimum matching label
    while(labels[y] != y)
        y = labels[y];
    // Second: Substitute the minmum matching label directly to labels[x] = y.
    while(labels[x] != x){
        int z = labels[x];
        labels[x] = y;
        x = z;
    }
    return y;
}

#undef cor
#undef rng
#endif // ____SwendsenWang_HK____ 