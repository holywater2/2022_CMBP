#include "Cluster_Long_range.hpp"
#include "../../headers/Writer.hpp"
#include <iostream>
#include <iomanip>

/***************** Parameters 1 *****************/
int kL          = 5;          /*Parameter: lattice size*/
int kN          = kL*kL;
int kBin        = 40;          /*Parameter: Change binning of temperature*/
double kB       = 0;
double kJ       = 1;
double alpha    = 3;
double Tsrt = 0;
double Tfin = 5;
double isTinf = false;
bool Random = false;
int equil_time_base = 1e4;
int equil_time = equil_time_base;
int mcs = 1e5;
/***************** Parameters 1 *****************/

typedef Cluster_LR_2D Model;

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
vector<string> result_to_file = vector<string>();
string output_path = "";


int main(int argn, char *argv[]){ // Input argument: argv[0]--> file name / argv[1]--> Input parameter)
    signal(SIGSEGV, &handler);
    signal(SIGINT, &handler);
    if(argn >= 2){
        string Input_file = argv[1];
        vector<double> input = Writer::Argument_reader(Input_file,11);
        kL = (int)input[0]; kN = kL*kL; kBin = (int)input[1]; kB = input[2]; kJ = input[3];
        alpha = input[4]; Tsrt = input[5]; Tfin = input[6];
        isTinf = input[7]; Random = input[8];
        equil_time_base = input[9]; mcs = input [10];        
    }
    if(argn >=3){
        output_path = argv[2];
    }
    // arguments list that helps to pass the args to model
    vector<double> args = {kL,kBin,kB,kJ,alpha,Tsrt,Tfin,isTinf,Random,equil_time_base,mcs};    
    // Filename Base: '\Result\(Model Name)_c_(kL)_int[erval]_(kBin) + (blahblah)
    #ifdef _WIN32
    static string kFilename = ".\\Result\\"+output_path+Model::Name()+"_c_"+to_string(kL)+"_int"+to_string(kBin)+"_mcs"+to_string(mcs)+"_a"+to_string(alpha);
    #endif
    #ifdef linux
    static string kFilename = "./Result/"+output_path+Model::Name()+"_c_"+to_string(kL)+"_int"+to_string(kBin)+"_mcs"+to_string(mcs)+"_a"+to_string(alpha);
    #endif
    Greetings();

    for(int gg = 0; gg < 1; gg++){
        Model model = Model(args);

        double MM, HH;
        double mcs_i = 1/double(mcs);
        double kNi = 1/double(kN);

        for(int i = 0; i < kBin; i++){
            model.Initialize(model.BetaV[i]);

            // Output data
            // 0: <m>, 1: <m^2>, 2: <m^4>, 3: <E>/sqrt(N), 4: <E^2>/N
            model.res = vector<double>(5,0);

            equil_time = equil_time_base;
            // if(model.TV[i]<=2.4 || model.TV[i]>=2.0) equil_time =1000;

            // model.IterateUntilEquilibrium(equil_time);

            duo value = model.Measure();
            HH = get<0>(value);
            MM = get<1>(value);

            cout <<"idx: " << left << setw(4) << i << "|| " << left << setw(10) << model.TV[i];
            cout << "|| "  << left << setw(9) << MM/(double)kN << "  " << left << setw(12) << HH << "|| ";
            
            model.IterateUntilEquilibrium(equil_time);


            /***********Monte Carlo Step and Caculate the data***********/
            for(int j = 0; j < mcs; j++){
                model.Calculate();                     //O(N^2)

                value = model.Measure_fast();
                HH = get<0>(value)/(double) kL;        // = E
                MM = abs(get<1>(value))/(double)kN;    // = |M|

                model.res[0] += MM*mcs_i;              // = <m>
                model.res[1] += (MM*mcs_i*MM);         // = <m^2>
                model.res[2] += (MM*mcs_i*MM)*(MM*MM); // = <m^4>
                model.res[3] += HH*mcs_i;              // = <E>/sqrt(N)
                model.res[4] += HH*mcs_i*HH;           // = <E^2>/N
            }
            /***********************************************************/

            /*******Calculate Magnetizaition and Specific Heat.*********/
            model.MV[i] = model.res[0];
            model.CV[i] = (model.BetaV[i]*model.BetaV[i])*(model.res[4]-model.res[3]*model.res[3]);
            /***********************************************************/

            cout << left << setw(13) << model.MV[i] << "  " << right << setw(13) << model.CV[i] << "|| ";
            cout << left << setw(14) << model.Fliped_Step << "  " << left << setw(10) << model.Total_Step << endl;

            string result = to_string(i) + "," + to_string(model.TV[i]) + "," + to_string(model.MV[i]) + "," + to_string(model.CV[i]) + ",";
            result = result + to_string(model.res[0]) + "," + to_string(model.res[1]) + "," + to_string(model.res[2]) + ",";
            result = result + to_string(model.res[3]) + "," + to_string(model.res[4]) + "\n";

            result_to_file.push_back(result);
        }

        /***********Save the result of the Calculation**********/
        Writer modelW = Writer(kFilename+"_Test_");
        modelW.WriteLine("idx,temperture,magnetization,specific heat,abs(mm),mm**2,mm**4,HH/L,HH**2/L\n");
        for(int i = 0; i < kBin; i++)
            modelW.WriteLine(result_to_file.at(i));
        modelW.CloseNewFile();
        /******************************************************/
    }
    Farewell();
}

// Memo
// 3.27 : at kL=64, alpha = 3, equil = 20, mcs = 100
// phase transition of magnetization occurs near 7.2 to 8.4

// 05.07
// 이 모델이 short range에서는 맞는지 먼저 확인하려고 했음. 그래서 alpha = 100으로 두고 5by5계산을 했는데
// SR metropolis 모델이랑 specific heat이 생각보다 크게 차이남. magnetization도 차이가 나긴하는데, 적고, specific heat은
// peak위치도 다르고, 값도 약 1.5배 차이남.
// 그래서 지금 가설이 1) 위 모델의 계산에서 Hamiltonian계산이 잘못 되고 있다. 혹은 mcs를 잘 못 나누고 있다. 2) 그냥 모델 자체의 로직이
// 잘못된거다. 로 생각하고 있음. 아마도 2번인거 같은데.. 어디서 잘못되고 있는지 모르겠음.

// 논문에서 제시한 새로운 partion function을 도입해야하나?