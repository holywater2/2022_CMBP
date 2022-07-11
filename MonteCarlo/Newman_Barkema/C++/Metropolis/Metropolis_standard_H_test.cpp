#include "Metropolis.hpp"
#include "../Writer.hpp"

#include <iostream>
#include <iomanip>

const int kL = 5; /*Parameter: lattice size*/
const int kN = kL*kL;
const int kBin = 3; /*Parametr: Change binning of temperature*/
// const int kB = 0;
double kB = 0;
const int kJ = 1;

const double Tsrt = T_CRIT*(1-0.2);
const double Tfin = T_CRIT*(1+0.2);

double isTinf = false;
const bool Random = true;

#ifdef _WIN32
static string kFilename = ".\\Result\\Metropolis_c_"+to_string(kL)+"_int"+to_string(kBin);
#endif
#ifdef linux
static string kFilename = "./Result/Metropolis_c_"+to_string(kL)+"_int"+to_string(kBin);
#endif

typedef Metropolis_2D Model;

// vector<double> args = {kL,kBin,kB,kJ,Tsrt,Tfin,isTinf};

clock_t __start__, __finish__;

void Greetings(){
    string Tat = isTinf ? "inf" : "0";

    cout << "Metropolis Algorithm\n";
    cout << "Radnomness test(seed): " << seed << '\n';
    cout << "L = " << kL << ", " << "bin = " << kBin << ", Start T at " << Tat << "\n";
    cout << "------------------------------------------------------------------------------------------------------------------" << "\n";
    cout << "--index--||---Temp----||EQ:sig------HH----------||magnetization---specific heat||Fliped Step------Total Step------" << "\n";
    cout << "------------------------------------------------------------------------------------------------------------------" << endl;
    cout << fixed <<setprecision(6);

    __start__ = clock();
}

void Farewell(int N = 0){
    __finish__ = clock();
    if(!N)
        cout << "\nProgram Abonormally Exit. Spent time: " << (double)(__finish__-__start__)/CLOCKS_PER_SEC << "\n";
    else
        cout << "Program Exit. Spent time: " << (double)(__finish__-__start__)/CLOCKS_PER_SEC << "\n";
    cout << "-------------------------------------------------------------------------------------------\n";
}

void handler(int A)
{
    cout << endl;
    Farewell(1);
    exit(A);
}

int main(){
    signal(SIGSEGV, &handler);
    signal(SIGINT, &handler);
    Greetings();
    vector<double> kBL = {-200,-100,-50,0,50,100,200};
    for(int gg = 0; gg < 1; gg++){
        kB = kBL[gg];
        vector<double> args = {kL,kBin,kB,kJ,Tsrt,Tfin,isTinf};

        Model model = Model(args);
        Writer modelW = Writer(kFilename+"_2e4_B"+to_string(kB));
        modelW.WriteLine("idx,temperture,magnetization,specific heat,abs(mm),mm**2,mm**4,HH/L,HH**2/L,<m>\n");

        int equil_time, mcs = 2e4;
        double MM, HH;
        double mcs_i = 1/double(mcs);
        double kNi = 1/double(kN);


        for(int i = 0; i < kBin; i++){
            model.Initialize(model.BetaV[i]);
            model.res = vector<double>(6,0);

            equil_time = 200;

            if(model.TV[i]<=2.4 || model.TV[i]>=2.0) equil_time = 1000;

            model.IterateUntilEquilibrium(2000);

            duo value = model.Measure();
            HH = get<0>(value);
            MM =  get<1>(value);

            cout <<"idx: " << left << setw(4) << i << "|| " << left << setw(10) << model.TV[i];
            cout << "|| "  << left << setw(9) << MM/(double)kN << "  " << left << setw(12) << int(HH) << "|| ";

            for(int j = 0; j < mcs; j++){ // mcs 2e4
                model.Calculate(0,Random);

                value = model.Measure_fast();
                HH = get<0>(value)/(double) kL;        // = E
                MM = abs(get<1>(value))/(double)kN;    // = M
                model.res[0] += MM*mcs_i;              // = <|m|>
                model.res[1] += (MM*mcs_i*MM);         // = <m^2>
                model.res[2] += (MM*mcs_i*MM)*(MM*MM); // = <m^4>
                model.res[3] += HH*mcs_i;              // = <E>/sqrt(N)
                model.res[4] += HH*mcs_i*HH;           // = <E^2>/N
                model.res[5] += get<1>(value)/(double)kN*mcs_i; // = <m>

            }
            model.MV[i] = model.res[0];
            model.CV[i] = (model.BetaV[i]*model.BetaV[i])*(model.res[4]-model.res[3]*model.res[3]);

            cout << left << setw(13) << model.MV[i] << "  " << right << setw(13) << model.CV[i] << "|| ";
            cout << left << setw(14) << model.Fliped_Step << "  " << left << setw(10) << model.Total_Step << "  " << left << setw(10) << model.res[5] << endl;

            string temp = to_string(i) + "," + to_string(model.TV[i]) + "," + to_string(model.MV[i]) + "," + to_string(model.CV[i]) + ",";
            temp = temp + to_string(model.res[0]) + "," + to_string(model.res[1]) + "," + to_string(model.res[2]) + ",";
            temp = temp + to_string(model.res[3]) + "," + to_string(model.res[4]) + "," + to_string(model.res[5]) + "\n";
            modelW.WriteLine(temp);
        }
        modelW.CloseNewFile();
    }
    Farewell();
}