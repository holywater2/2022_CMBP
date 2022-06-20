#include "Wolff.hpp"
#include "../Writer.hpp"

#include <iostream>
#include <iomanip>

#include <signal.h>
#include <cstdlib>


// const int kL = 100; /*Parameter: lattice size*/
// const int kN = kL*kL;
int kL = 0;
int kN = 0;
const int kBin = 25; /*Parametr: Change binning of temperature*/
const int kB = 0;
const int kJ = 1;

// T_crit ~ 2.269
const double Tsrt = T_CRIT*(1-0.1);
const double Tfin = T_CRIT*(1+0.1);

double isTinf = false;

#ifdef _WIN32
static string kFilename = ".\\Result\\Wolff_c_"+to_string(kL)+"_int"+to_string(kBin);
#endif //_WIN32
#ifdef linux
static string kFilename = "./Result/Wolff_c_"+to_string(kL)+"_int"+to_string(kBin);
#endif //linux

typedef Wolff_2D Model;

vector<double> args = {(double) kL,kBin,kB,kJ,Tsrt,Tfin,isTinf};

clock_t __start__, __finish__;

void Greetings(){
    string Tat = isTinf ? "inf" : "0";

    cout << "Wolff Algorithm\n";
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
        cout << "Program Exit Exit. Spent time: " << (double)(__finish__-__start__)/CLOCKS_PER_SEC << "\n";
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

    vector<int> kLL = {24,32,48,64,96};

    for(int gg = 4; gg < 5; gg++){
        kL = kLL[gg];
        kN = kL*kL;
    
        vector<double> args = {(double)kL,kBin,kB,kJ,Tsrt,Tfin,isTinf};
        Model model = Model(args);
        Writer modelW = Writer(kFilename+"_binder1e5_"+to_string(kL));
        modelW.WriteLine("idx,temperture,magnetization,specific heat,abs(MM),MM**2,MM**4,HH,HH**2,c_error,m_error\n");

        int HH, equil_time, mcs = 1e6;
        double MM;
        double mcs_i = 1/double(mcs);
        double kNi = 1/double(kN);

        cout << showpos;

        for(int i = 0; i < kBin; i++){
            model.Initialize(model.BetaV[i]);
            model.res = vector<double>(5,0);

            equil_time = 50;

            if(model.TV[i]<2.4 || model.TV[i]>2.0) equil_time = 100;

            model.IterateUntilEquilibrium(equil_time);

            duo value = model.Measure();
            HH = get<0>(value);
            MM =  get<1>(value);

            cout <<"idx: " << left << setw(4) << i << "|| " << left << setw(10) << model.TV[i];
            cout << "|| "  << left << setw(9) << MM/(double)kN << "  " << left << setw(12) << HH << "|| ";

            int jB = 1000;
            vector<double> Jackknife = vector<double>(mcs/jB,0);
            vector<double> Jackknife_HH = vector<double>(mcs/jB,0);
            vector<double> Jackknife_HH2 = vector<double>(mcs/jB,0);
            double c_error = 0;

            double size, step;
            for(double j = 0; j < mcs;){ // mcs 4000
                size = model.Calculate();
                step = size/(double)kN;
                j += step;
                step /= mcs;

                value = model.Measure();
                HH = get<0>(value);
                MM =  get<1>(value)/(double)kN;
                
                model.res[0] += abs(MM)*step;
                model.res[1] += (MM*step*MM);
                model.res[2] += (MM*step*MM)*(MM*MM);
                model.res[3] += HH*step/(double)kL;
                model.res[4] += HH*step*HH/(double)kN;

                Jackknife_HH[j/jB] += HH;
                Jackknife_HH2[j/jB] += HH*HH;
            }
            model.MV[i] = model.res[0];
            model.CV[i] = (model.BetaV[i]*model.BetaV[i])*(model.res[4]-model.res[3]*model.res[3]);        

            double avg_Step_Size = (model.Fliped_Step/(double) model.Total_Step)/kN;

            for(int j = 0; j < mcs/jB; j++){
                Jackknife_HH[j]  *= avg_Step_Size/kL;
                Jackknife_HH2[j] *= avg_Step_Size/kN;
                Jackknife[j] = (model.BetaV[i]*model.BetaV[i])*((mcs*model.res[4]-Jackknife_HH2[j])/(mcs-jB) \
                            - (mcs*model.res[3]-Jackknife_HH[j])/(mcs-jB)*(mcs*model.res[3]-Jackknife_HH[j])/(mcs-jB));
                c_error += (Jackknife[j]-model.CV[i])*(Jackknife[j]-model.CV[i]);
            }
            c_error = sqrt(c_error);

            cout << left << setw(13) << model.MV[i] << "  " << right << setw(13) << model.CV[i] << "|| ";
            cout << left << setw(14) << model.Fliped_Step << "  " << left << setw(10) << model.Total_Step << "||" << left << setw(10) << c_error <<endl;

            string temp = to_string(i) + "," + to_string(model.TV[i]) + "," + to_string(model.MV[i]) + "," + to_string(model.CV[i]) + ",";
            temp = temp + to_string(model.res[0]) + "," + to_string(model.res[1]) + "," + to_string(model.res[2]) + ",";
            temp = temp + to_string(model.res[3]) + "," + to_string(model.res[4]) + "," + to_string(c_error) + "\n";
            modelW.WriteLine(temp);
        }
        modelW.CloseNewFile();
    }
    Farewell();
}
