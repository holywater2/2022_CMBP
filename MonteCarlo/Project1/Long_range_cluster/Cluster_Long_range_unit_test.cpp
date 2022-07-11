#include "Cluster_Long_range.hpp"
#include "../../headers/Writer.hpp"
#include <iostream>
#include <iomanip>

/***************** Parameters 1 *****************/
const int kL          = 32;          /*Parameter: lattice size*/
const int kN          = kL*kL;
const int kBin        = 20;          /*Parameter: Change binning of temperature*/
const int kB          = 0;
const int kJ          = 1;
const double alpha    = 2+1;

// const double Tsrt = T_CRIT*(1-0.08);
// const double Tfin = T_CRIT*(1+0.08);

const double Tsrt = 5.0;
const double Tfin = 5.5;

const double isTinf = false;
const bool Random = true;

int equil_time = 1000;
const int equil_time_base = 200;
int mcs = 1e4;
/***************** Parameters 1 *****************/

typedef Cluster_LR_2D Model;

// Filename Base: '\Result\(Model Name)_c_(kL)_int[erval]_(kBin) + (blahblah)
#ifdef _WIN32
static string kFilename = ".\\Result\\"+Model::Name()+"_c_"+to_string(kL)+"_int"+to_string(kBin)+"_mcs"+to_string(mcs);
#endif
#ifdef linux
static string kFilename = "./Result/"+Model::Name()+"_c_"+to_string(kL)+"_int"+to_string(kBin)+"_mcs"+to_string(mcs);
#endif

// clock used to measure time
clock_t __start__, __finish__;

void Greetings(){
    string Tat = isTinf ? "inf" : "0";

    cout << Model::Name() + "Algorithm\n";
    cout << "Radnomness test(seed): " << seed << '\n';
    cout << "L = " << kL << ", " << "bin = " << kBin << ", Start T at " << Tat << "\n";
    cout << "------------------------------------------------------------------------------------------------------------------" << "\n";
    cout << "--index--||---Temp----||EQ:sig------HH----------||magnetization---specific heat||Fliped Step------Total Step------" << "\n";
    cout << "------------------------------------------------------------------------------------------------------------------" << endl;
    cout << fixed <<setprecision(6);
    
    // Show +/- sign
    cout << showpos;


    __start__ = clock();
}

void Farewell(int N = 0){
    __finish__ = clock();
    if(N)
        cout << "\nProgram Abonormally Exit. Spent time: " << (double)(__finish__-__start__)/CLOCKS_PER_SEC << "\n";
    else
        cout << "Program Exit. Spent time: " << (double)(__finish__-__start__)/CLOCKS_PER_SEC << "\n";
    cout << "-------------------------------------------------------------------------------------------\n";
}

// Event handler that notify the spent time when program abonormally stop
void handler(int A){
    cout << endl;
    Farewell(1);
    exit(A);
}

// arguments list that helps to pass the args to model
vector<double> args = {kL,kBin,kB,kJ,alpha,Tsrt,Tfin,isTinf};
int main(){
    signal(SIGSEGV, &handler);
    signal(SIGINT, &handler);
    Greetings();

    for(int gg = 0; gg < 1; gg++){
        Model model = Model(args);
        model.Initialize(model.BetaV[0]);
        // model.IterateUntilEquilibrium(equil_time);

        cout << model.J_tot << '\n';

        double MM, HH;
        double mcs_i = 1/double(mcs);
        double kNi = 1/double(kN);

        // cout << "blah" << '\n';
        
        int a = 0;
        double cnt = 10;
        poisson_distribution<int> d(40000);
        static mt19937 gen(123124); // Standard mersenne_twister_engine seeded with time()
        static uniform_real_distribution<> dis(0.0, 1.0);
        for(int i = 0; i < cnt; i++){
            int temp = d(gen);
            // int temp = model.PoissonNumberGenerator(40000);
            cout << temp << endl;
            a += temp;        
        }
        cout << a/cnt << '\n';
    }
    Farewell();
}

// int main(){
//     signal(SIGSEGV, &handler);
//     signal(SIGINT, &handler);
//     Greetings();

//     for(int gg = 0; gg < 1; gg++){
//         Model model = Model(args);
//         model.Initialize(model.BetaV[0]);
//         model.IterateUntilEquilibrium(equil_time);

//         cout << model.J_tot << '\n';

//         double MM, HH;
//         double mcs_i = 1/double(mcs);
//         double kNi = 1/double(kN);

//         cout << "blah" << '\n';
//         int N = kN;
//         int iter_k = 15;

//         vector<short> bond_list = vector<short>(kN*(N-1)/2,0); // hashset을 써야하나?
//         vector<set<int>> adj = vector<set<int>>(kN);
//         for(int it = 0; it < iter_k; it++){ // O(lambda)
//             // Walker's Alias Method
//             int l = ((kN*(kN-1))/2)*dis(gen);    
//             if(dis(gen) > model.Walker_Table_P[l])
//                 l = model.Walker_Table_A[l];

//             int i = model.i_of_bond[l];
//             int j = l- (N*i - i*(i+1)/2 + (-i-1));
//             cout << i << " " << j << " ";
//             cout << model.sc[i] << " " << model.sc[j] << '\n';


//             if(model.sc[i]*model.sc[j] == 1){
//                 adj[i].insert(j);
//                 adj[j].insert(i);
//                 // bond_list[l]++;
//             }
//             // cout << i << ' ' << j << '\n';
//         }

//         // cout << endl;

//         vector<bool> visited = vector<bool>(N,0); // 0 = false
//         queue<int> que;
//         for(int i = 0; i < N; i++){ //O(N+lambda) = O(N+E)
//             cout << "i is " << i << "\n";
//             if(visited[i] == true)
//                 continue;
//             if(adj[i].empty()){
//                 cout << "i is empty" << '\n';
//                 visited[i] = true;
//                 continue;
//             }

//             int mul = 2*(int)(dis(gen)*2) -1;

//             que.push(i);
//             visited[i] = true;
//             // BFS
//             while(!que.empty()){
//                 int s = que.front();
//                 que.pop();
                
//                 // cout << s << endl;

//                 model.sc[s] *= mul; // Flip
//                 cout << s << '\n';

//                 auto itr = adj[s].begin();
//                 while(itr != adj[s].end()){
//                     int n = *itr++;
//                     if(!visited[n]){
//                         visited[n] = true;
//                         que.push(n);
//                     }
//                 }
//             }
//         }
//     }
//     Farewell();
// }

// int main(){
//     signal(SIGSEGV, &handler);
//     signal(SIGINT, &handler);
//     Greetings();

//     for(int gg = 0; gg < 1; gg++){
//         Model model = Model(args);
//         model.Initialize(model.BetaV[0]);
//         model.IterateUntilEquilibrium(equil_time);

//         cout << model.J_tot << '\n';

//         double MM, HH;
//         double mcs_i = 1/double(mcs);
//         double kNi = 1/double(kN);

//         cout << "blah" << '\n';
//         int N = kN;
//         int iter_k = 15;

//         vector<short> bond_list = vector<short>(kN*(N-1)/2,0); // hashset을 써야하나?
//         vector<set<int>> adj = vector<set<int>>(kN);
//         for(int it = 0; it < iter_k; it++){ // O(lambda)
//             // Walker's Alias Method
//             int l = ((kN*(kN-1))/2)*dis(gen);    
//             if(dis(gen) > model.Walker_Table_P[l])
//                 l = model.Walker_Table_A[l];

//             int i = model.y_of_bond[l];
//             int j = l- (N*i - i*(i+1)/2 + (-i-1));
//             cout << i << " " << j << " ";
//             cout << model.sc[i] << " " << model.sc[j] << '\n';


//             if(model.sc[i]*model.sc[j] == 1){
//                 adj[i].insert(j);
//                 adj[j].insert(i);
//                 // bond_list[l]++;
//             }
//             // cout << i << ' ' << j << '\n';
//         }

//         // cout << endl;

//         vector<bool> visited = vector<bool>(N,0); // 0 = false
//         queue<int> que;
//         for(int i = 0; i < N; i++){ //O(N+lambda) = O(N+E)
//             cout << "i is " << i << "\n";
//             if(visited[i] == true)
//                 continue;
//             if(adj[i].empty()){
//                 cout << "i is empty" << '\n';
//                 visited[i] = true;
//                 continue;
//             }

//             int mul = 2*(int)(dis(gen)*2) -1;

//             que.push(i);
//             visited[i] = true;
//             // BFS
//             while(!que.empty()){
//                 int s = que.front();
//                 que.pop();
                
//                 // cout << s << endl;

//                 model.sc[s] *= mul; // Flip
//                 cout << s << '\n';

//                 auto itr = adj[s].begin();
//                 while(itr != adj[s].end()){
//                     int n = *itr++;
//                     if(!visited[n]){
//                         visited[n] = true;
//                         que.push(n);
//                     }
//                 }
//             }
//         }
//     }
//     Farewell();
// }

// int main(){
//     signal(SIGSEGV, &handler);
//     signal(SIGINT, &handler);
//     Greetings();

//     for(int gg = 0; gg < 1; gg++){
//         Model model = Model(args);
//         cout << model.J_tot << '\n';

//         double MM, HH;
//         double mcs_i = 1/double(mcs);
//         double kNi = 1/double(kN);

//         cout << "blah" << '\n';
//         for(int i = 0; i < 1; i++){
//             model.Initialize(model.BetaV[i]);
//             double t = 0;

//             for(int j = 0; j < kL; j++){
//             for(int k = 0; k < kL; k++){
//                 cout << model.sc[j*kL+k];
//             }
//             cout << '\n';
//             }

//             for(int mm = 0; mm < kN*(kN-1)/2; mm++){
//                 if(model.Walker_Table_A[mm] != -1)
//                     t += 1;
//                 else
//                     t += model.Walker_Table_P[mm];
//             }
//             cout << "total Prob is " << t << endl;
//             model.Calculate();

//             for(int j = 0; j < kL; j++){
//             for(int k = 0; k < kL; k++){
//                 cout << model.sc[j*kL+k];
//             }
//             cout << '\n';
//             }
//             int div = 1e9;
//             vector<double> simulator = vector<double>(kN*(kN-1),0);
//             for(int kkk = 0; kkk < div; kkk++){
//                 int l = ((kN*(kN-1))/2)*dis(gen);
//                 if(dis(gen) > model.Walker_Table_P[l])
//                     l = model.Walker_Table_A[l];
//                 simulator[l]++;
//             }
//             for(int j = 0; j < kN; j++){
//             for(int k = j+1; k < kN; k++){
//                 cout << kN*j - j*(j+1)/2 + (k-j-1) << '\n';
//                 // cout << simulator[kN*j - j*(j+1)/2 + (k-j-1)]/1e8 << '\n' << pow(model.e2d.dist_ij(j,k),-alpha)/model.J_tot << '\n';
//                 // cout << simulator[kN*j - j*(j+1)/2 + (k-j-1)]/1e8 - pow(model.e2d.dist_ij(j,k),-alpha)/model.J_tot << '\n';
                
//                 cout << simulator[kN*j - j*(j+1)/2 + (k-j-1)]/div << '\n' << model.e2d.pi_ij(j,k)/model.J_tot << '\n';
//                 cout << simulator[kN*j - j*(j+1)/2 + (k-j-1)]/div - model.e2d.pi_ij(j,k)/model.J_tot << '\n';
//             }
//             cout << '\n';
//             }
//         }
//     }
//     Farewell();
// }