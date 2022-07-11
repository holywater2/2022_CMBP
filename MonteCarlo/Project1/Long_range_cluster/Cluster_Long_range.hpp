#ifndef ____CLUSTER_LR____
#define ____CLUSTER_LR____ 

/**
 * @file Cluster_Long_range.hpp
 * @author Seongsu, Kim (holywater@gist.ac.kr)
 * @brief 
 * @version 0.1
 * @date 2022-05-05
 * 
 * @copyright Copyright (c) 2022
 * 
 * Ref1: Luijten and Blöte method
 */

// File & IO System
#include <iostream>
// Data Structure
#include <vector>
#include <tuple>
#include <string>
#include <tuple>
#include <set>
#include <algorithm>
#include <queue>
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

class Cluster_LR_2D{

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
        double J_tot;
        vector<int> i_of_bond;
        vector<long double> Walker_Table_P;
        vector<int> Walker_Table_A;

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

        static string Name(){return "Cluster_LR";}
        Cluster_LR_2D(int L, int bin, double B, double J, double alpha, double Tsrt, double Tfin, bool isTinf);
        Cluster_LR_2D(vector<double> args);
        ~Cluster_LR_2D(){
            // delete e2d;
            __finish__ = clock();
            cout << "------------------------------------------------------------------------------------------------------------------\n";
            cout << "Model calculation finished. Spent time: " << (double)(__finish__-__start__)/CLOCKS_PER_SEC << "\n";
            cout << "------------------------------------------------------------------------------------------------------------------\n";
        }
        double Prob(double delta);
        void Initialize(double beta);
        // void Initialzie(int idx);
        int SweepHelical(int i);
        int BoundaryHelical(int i);
        duo Measure();
        duo Measure_fast();
        void Calculate(int _n = 0,bool Random = true);
        void IterateUntilEquilibrium(int equil_time,bool random = true);

        int PoissonNumberGenerator(double lambdaTot);
        double JTot();
        void WalkerTableGenerator(double eps = 1e-10);
};

Cluster_LR_2D::Cluster_LR_2D(int L, int bin, double B, double J, double alpha, double Tsrt, double Tfin, bool isTinf)
:L(L), N(L*L), Bin(bin), B(B), J(J), YNN(L), alpha(alpha) {
    this-> isTinf = isTinf;
    this-> sc = new short[N];

    this-> MV = vector<double>(Bin);
    this-> CV = vector<double>(Bin);
    this-> TV = vector<double>(Bin);
    this-> BetaV = vector<double>(Bin);
    this-> i_of_bond = vector<int>(N*(N-1)/2);

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
    // Ewald sum in 2D
    e2d = ewald_ND(2,vector<int>({L,L}),alpha);

    J_tot = JTot();
    WalkerTableGenerator();
}

Cluster_LR_2D::Cluster_LR_2D(vector<double> args):
Cluster_LR_2D(args[0],args[1],args[2],args[3],args[4],args[5],args[6],args[7]){}

double Cluster_LR_2D::Prob(double delta){
    return exp(-this->cur_beta*delta);
}

// https://kr.mathworks.com/matlabcentral/answers/28161-poisson-random-number-generator
// Mean이 lambdatot보다 1크게 나옴 -> -1 더해서 사용
// int Cluster_LR_2D::PoissonNumberGenerator(double lambdaTot){
//     int k = 1;
//     long double prod = dis(gen);
//     long double exp_lambda = expl(-lambdaTot);
//     while(prod > exp_lambda){
//         prod *= dis(gen);
//         k += 1;
//     }
//     return k;
// }
// https://en.wikipedia.org/wiki/Poisson_distribution#Random_drawing_from_the_Poisson_distribution
// 폐기처분 -> poisson_distribution<int> d(lambda 사용);
int Cluster_LR_2D::PoissonNumberGenerator(double lambdaTot){
    int k = 0;
    long double prod = dis(gen);
    long double s = expl(-lambdaTot);
    long double exp_lambda = -lambdaTot;
    while(prod > s){
        k += 1;
        exp_lambda += logl(lambdaTot/k);
        s += expl(exp_lambda);
    }
    return k;
}


double Cluster_LR_2D::JTot(){ // Sum of J_ij where j > i
    double J_tot = 0;
    for(int i = 0;   i < N; i++)
    for(int j = i+1; j < N; j++){       
        i_of_bond[N*i - i*(i+1)/2 + (j-i-1)] = i;
        // cout << N*i - i*(i+1)/2 + (j-i-1) << '\n';
        // J_tot += pow(e2d.dist_ij(i,j),-alpha);
        J_tot += e2d.pi_ij(i,j);
    }
    cout << J_tot << '\n';
    return J_tot*J;
}

void Cluster_LR_2D::Initialize(double beta){
    this-> Calc_call = 0;
    this-> Fliped_Step = 0;
    this-> Total_Step = 0;

    for(int i = 0; i < N; i++){
        // T = 0 start
        sc[i] = 1;

        // T = \inf start
        if(this->isTinf) this->sc[i] -= int(dis(gen)*2)*2;
    }
    cur_beta = beta;

    this->Measure();
}

int Cluster_LR_2D::SweepHelical(int i){
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

int Cluster_LR_2D::BoundaryHelical(int i){
    int nn, sum = 0;
    // int XNN = 1, YNN = L;
    
    if((nn = i + XNN) == N) nn = 0;
    sum += this->sc[nn];
    if((nn = i + YNN) >= N) nn -= N;
    sum += this->sc[nn];

    return sum;
}

duo Cluster_LR_2D::Measure(){ //O(N^2)
    int i, j;
    int sigma = 0;
    double HH, result = 0;
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

duo Cluster_LR_2D::Measure_fast(){
    // Measure_fast is not implemneted
    return Measure();
}
void Cluster_LR_2D::Calculate(int _n, bool Random){
    double lambda_tot = 2*cur_beta*J_tot;
    // cout << lambda_tot << '\n';
    poisson_distribution<int> PNG(lambda_tot);
    int iter_k = PNG(lambda_tot);
    // cout << "iter_k " << iter_k << endl;
    vector<short> bond_list = vector<short>(N*(N-1)/2,0); // hashset을 써야하나?
    vector<set<int>> adj = vector<set<int>>(N);
    for(int it = 0; it < iter_k; it++){ // O(lambda)
        // Walker's Alias Method
        int l = ((N*(N-1))/2)*dis(gen);    
        if(dis(gen) > Walker_Table_P[l])
            l = Walker_Table_A[l];

        int i = i_of_bond[l];
        int j = l- (N*i - i*(i+1)/2 + (-i-1));

        if(sc[i]*sc[j] == 1){
            adj[i].insert(j);
            adj[j].insert(i);
            // bond_list[l]++;
        }
        // cout << i << ' ' << j << '\n';
    }

    // cout << endl;

    vector<bool> visited = vector<bool>(N); // 0 = false
    queue<int> que;
    for(int i = 0; i < N; i++){ //O(N+lambda) = O(N+E)
        // cout << "i is " << i << "\n";
        if(visited[i] == true)
            continue;
        if(adj[i].empty()){
            // cout << "i is empty" << '\n';
            visited[i] = true;
            continue;
        }

        int mul = 2*(int)(dis(gen)*2) -1;

        que.push(i);
        visited[i] = true;
        // BFS
        while(!que.empty()){
            int s = que.front();
            que.pop();
            
            // cout << s << endl;

            sc[s] *= mul; // Flip
            this-> Total_Step++;
            if(mul == -1)
                this->Fliped_Step++;

            auto itr = adj[s].begin();
            while(itr != adj[s].end()){
                int n = *itr++;
                if(!visited[n]){
                    visited[n] = true;
                    que.push(n);
                }
            }
        }
    }
}

void Cluster_LR_2D::WalkerTableGenerator(double eps){ // 제대로 작동하고 있는 듯 (Test 완료)
    this->Walker_Table_P = vector<long double>(N*(N-1)/2);
    this->Walker_Table_A = vector<int>(N*(N-1)/2,-1);
    vector<int> overfull;
    vector<int> underfull;

    // double a = 0;
    for(int i = 0;   i < N; i++)
    for(int j = i+1; j < N; j++){
        // double prob = pow(e2d.dist_ij(i,j),-alpha)*N*(N-1)/2/J_tot;
        long double prob = J*e2d.pi_ij(i,j)*N*(N-1)/2/J_tot;

        int cur_point = N*i - i*(i+1)/2 + (j-i-1);
        // cout << cur_point << " " << prob << '\n';
        if(prob > 1) overfull.push_back(cur_point);
        else underfull.push_back(cur_point);
        Walker_Table_P[cur_point] = prob;
        // a += prob;
    }

    // cout << "a is " << a << "and N*(N-1)/2 is " << N*(N-1)/2 <<'\n';


    // vector<int>::iterator iter_over  = overfull.begin();
    // vector<int>::iterator iter_under = underfull.begin();
    int iter_over = 0;
    int iter_under = 0;

    
    while(iter_over < overfull.size() && iter_under < underfull.size()){
        long double& P_over = Walker_Table_P[overfull[iter_over]];
        long double& P_under = Walker_Table_P[underfull[iter_under]];
        // Case 1: P(iter_over) = 1 + (del + alpha) & P(iter_under) = 1 - (del)
        // Case 2: P(iter_over) = 1 + (del) & P(iter_under) = 1 - (del)
        // Case 3: P(iter_over) = 1 + (del) & P(iter_under) = 1 - (del + alpha)
        // cout << overfull[iter_over] << " " << P_over << '\n';
        // cout << underfull[iter_under] << " " << P_under << '\n';

        P_over = P_over - (1 - P_under);
        Walker_Table_A[underfull[iter_under]] = overfull[iter_over];
        // cout << overfull[iter_over] << " " << P_over << '\n';
        // cout << underfull[iter_under] << " " << P_under << '\n';
        // cout << '\n';
        
        iter_under++;
        if(P_over < 1+eps && P_over > 1-eps){ // == , equal case
            P_over = 1;
            iter_over++;
        } else if(P_over < 1){
            underfull.push_back(overfull[iter_over]);
            iter_over++;
        }
    }
    // Iterator 사용하면서 element 넣으면 오류가 날 수 있음 -> address reallocation
}

void Cluster_LR_2D::IterateUntilEquilibrium(int equil_time,bool random){
    for(int j = 0; j < equil_time; j++)
        Calculate(0,random);
}

#endif // ____CLUSTER_LR____ 