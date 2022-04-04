#include "Metropolis_Long_range.hpp"
#include "../../headers/Writer.hpp"
#include <iostream>
#include <iomanip>

/***************** Parameters 1 *****************/
const int kL          = 128;          /*Parameter: lattice size*/
const int kN          = kL*kL;
const int kBin        = 25;          /*Parameter: Change binning of temperature*/
const double kB       = 0;
const double kJ       = 1;
const double alpha    = 3;


// const double Tsrt = T_CRIT*(1-0.08);
// const double Tfin = T_CRIT*(1+0.08);

const double Tsrt = 6.4;
const double Tfin = 8.4;
// const double Tsrt = 0;
// const double Tfin = 8;


const double isTinf = false;
const bool Random = false;

const int equil_time_base = 20;
int equil_time = equil_time_base;
int mcs = 1e4;
/***************** Parameters 1 *****************/

typedef Metropolis_LR_2D Model;

// Filename Base: '\Result\(Model Name)_c_(kL)_int[erval]_(kBin) + (blahblah)
#ifdef _WIN32
static string kFilename = ".\\Result\\"+Model::Name()+"_c_"+to_string(kL)+"_int"+to_string(kBin);
#endif
#ifdef linux
static string kFilename = "./Result/"+Model::Name()+"_c_"+to_string(kL)+"_int"+to_string(kBin);
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
        Writer modelW = Writer(kFilename+"_Test_");
        modelW.WriteLine("idx,temperture,magnetization,specific heat,abs(mm),mm**2,mm**4,HH/L,HH**2/L\n");

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

            model.IterateUntilEquilibrium(equil_time);

            duo value = model.Measure();
            HH = get<0>(value);
            MM = get<1>(value);

            cout <<"idx: " << left << setw(4) << i << "|| " << left << setw(10) << model.TV[i];
            cout << "|| "  << left << setw(9) << MM/(double)kN << "  " << left << setw(12) << HH << "|| ";

            // Monte Carlo Step and Caculate the data
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

            // Calculate Magnetizaition and Specific Heat.
            model.MV[i] = model.res[0];
            model.CV[i] = (model.BetaV[i]*model.BetaV[i])*(model.res[4]-model.res[3]*model.res[3]);

            cout << left << setw(13) << model.MV[i] << "  " << right << setw(13) << model.CV[i] << "|| ";
            cout << left << setw(14) << model.Fliped_Step << "  " << left << setw(10) << model.Total_Step << endl;

            string temp = to_string(i) + "," + to_string(model.TV[i]) + "," + to_string(model.MV[i]) + "," + to_string(model.CV[i]) + ",";
            temp = temp + to_string(model.res[0]) + "," + to_string(model.res[1]) + "," + to_string(model.res[2]) + ",";
            temp = temp + to_string(model.res[3]) + "," + to_string(model.res[4]) + "\n";
            modelW.WriteLine(temp);
        }
        modelW.CloseNewFile();
    }
    Farewell();
}

// Memo
// 3.27 : at kL=64, alpha = 3, equil = 20, mcs = 100
// phase transition of magnetization occurs near 7.2 to 8.4