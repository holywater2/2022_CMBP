#include "AA_Metropolis.hpp"
#include "../../../headers/Writer.hpp"
#include <iostream>
#include <iomanip>

/***************** Parameters 1 *****************/
int kLx         = 5;          /*Parameter: lattice size*/
int kLy         = 5;          

int kN          = kLx*kLy;
int kBin        = 40;          /*Parameter: Change binning of temperature*/

double kB       = 0;
double kJx      = 1;
double kJy      = 1;
double alpha    = 100;

double Tsrt = 0;
double Tfin = 5;

double isTinf = false;
bool Random = false;

int equil_time_base = 1e5;
int equil_time = equil_time_base;
int mcs = 1e6;
/***************** Parameters 1 *****************/

typedef AA_Metropolis Model;

// clock used to measure time
clock_t __start__, __finish__;

void Greetings(){
    string Tat = isTinf ? "inf" : "0";

    cout << Model::Name() + "Algorithm\n";
    cout << "Radnomness test(seed): " << seed << '\n';
    cout << "L = " << kLx << "," << kLy << ", " << "bin = " << kBin << ", Start T at " << Tat <<  "\n";
    cout << "N = " << kN << ", " << "alpha = " << alpha << ", mcs " << mcs <<  "\n";
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

int main(int argn, char *argv[]){ // Input argument: argv[0]--> file name / argv[1]--> Input parameter
    signal(SIGSEGV, &handler);
    signal(SIGINT, &handler);
    if(argn >= 2){
        string Input_file = argv[1];
        vector<double> input = Writer::Argument_reader(Input_file,11);
        kLx = (int)input[0]; kLy = (int)input[1]; kN = kLx*kLy;
        kBin = (int)input[2]; kB = input[3]; kJx = input[4]; kJy = input[5];
        alpha = input[6]; Tsrt = input[7]; Tfin = input[8];
        isTinf = input[9]; Random = input[10];
        equil_time_base = input[11]; mcs = input[12];        
    }
    // arguments list that helps to pass the args to model
    vector<double> args = {kLx,kLy,kBin,kB,kJx,kJy,alpha,Tsrt,Tfin,isTinf,Random,equil_time_base,mcs};

    // Filename Base: '\Result\(Model Name)_c_(kL)_int[erval]_(kBin) + (blahblah)
    #ifdef _WIN32
    static string kFilename = ".\\Result\\"+Model::Name()+"_c_"+to_string(kLx)+"_"+to_string(kLy)+"_int"+to_string(kBin)+"_mcs"+to_string(mcs)+"_a"+to_string(alpha);
    #endif
    #ifdef linux
    static string kFilename = "./Result/"+Model::Name()+"_c_"+to_string(kLx)+"_"+to_string(kLy)+"_int"+to_string(kBin)+"_mcs"+to_string(mcs)+"_a"+to_string(alpha);
    #endif

    Greetings();

    for(int gg = 0; gg < 1; gg++){
        Model model = Model(args);

        double MM, HH;
        double mcs_i = 1/double(mcs);
        double kNi  = 1/double(kN);
        double kNir = pow(kN,0.5);
        cout << kNir << '\n';

        for(int i = 0; i < kBin; i++){
            model.Initialize(model.BetaV[i]);

            // Output data
            // 0: <m>, 1: <m^2>, 2: <m^4>, 3: <E>/sqrt(N), 4: <E^2>/N
            model.res = vector<double>(5,0);

            equil_time = equil_time_base;
            // if(model.TV[i]<=2.4 || model.TV[i]>=2.0) equil_time =1000;

            model.IterateUntilEquilibrium(equil_time);

            duo value = model.Measure();
            HH = get<0>(value);
            MM = get<1>(value);

            cout <<"idx: " << left << setw(4) << i << "|| " << left << setw(10) << model.TV[i];
            cout << "|| "  << left << setw(9) << MM/(double)kN << "  " << left << setw(12) << HH << "|| ";

            /***********Monte Carlo Step and Caculate the data***********/
            for(int j = 0; j < mcs; j++){
                model.Calculate();                     //O(N^2)

                value = model.Measure_fast();
                HH = get<0>(value)/(double) kNir;        // = E
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