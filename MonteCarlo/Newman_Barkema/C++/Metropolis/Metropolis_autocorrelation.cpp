#include "Metropolis.hpp"
#include "../Writer.hpp"

#include <iostream>
#include <iomanip>

const int kL = 100; /*Parameter: lattice size*/
const int kN = kL*kL;
const int kBin = 25; /*Parametr: Change binning of temperature*/
const int kB = 0;
const int kJ = 1;

const double Tsrt = 0;
const double Tfin = 5;

double isTinf = false;

#ifdef _WIN32
static string kFilename = ".\\Result\\Metropolis_c_Error_"+to_string(kL)+"_int"+to_string(kBin);
#endif
#ifdef linux
static string kFilename = "./Result/Metropolis_c_Error_"+to_string(kL)+"_int"+to_string(kBin);
#endif

typedef Metropolis_2D Model;

vector<double> args = {kL,kBin,kB,kJ,Tsrt,Tfin,isTinf};

clock_t __start__, __finish__;

void Greetings(){
    string Tat = isTinf ? "inf" : "0";

    cout << "Metropolis Algorithm\n";
    cout << "Radnomness test(seed): " << seed << endl;
    cout << "L = " << kL << ", " << "bin = " << kBin << ", Start T at " << Tat << "\n";
    cout << "-------------------------------------------------------------------------------------------" << "\n";
    cout << "index---||-Temp-||-sig--HH-||magnetization-specific heat--Fliped--Total Step---------------" << "\n";
    cout << "-------------------------------------------------------------------------------------------" << endl;
    cout << setprecision(6);

    __start__ = clock();
}

void Farewell(){
    __finish__ = clock();
    cout << "Program Exit. Spent time: " << (double)(__finish__-__start__)/CLOCKS_PER_SEC << "\n";
    cout << "-------------------------------------------------------------------------------------------";
}

int main(){
    Greetings();
    
    for(int i = 0; i < 8; i++){
        Model model = Model(args);
        Writer modelW = Writer(kFilename);
        Writer modelW2 = Writer(kFilename + "_auto");
        modelW.WriteLine("idx,temperture,magnetization,specific heat,abs(sigma),sigma**2,sigma**4,HH,HH**2,m_error\n");

        int sigma, HH, equil_time, mcs = 10000;
        double mcs_i = 1/double(mcs);
        double kNi = 1/double(kN);

        cout << showpos;

        for(int i = 0; i < kBin; i++){
            model.Initialize(model.BetaV[i]);
            model.res = vector<double>(5,0);

            duo value = model.Measure();
            HH = get<0>(value);
            sigma =  get<1>(value);

            cout <<"idx: " << i << "\t|| " << model.TV[i] << "\t|| " << sigma << "\t" << HH << setw(3) << "|| ";
            modelW2.WriteLine(model.TV[i]);

            for(int j = 0; j < mcs; j++){
                model.Calculate(0,false);

                value = model.Measure();
                HH = get<0>(value);
                sigma =  get<1>(value);
                
                model.res[0] += abs(sigma)*mcs_i;
                model.res[1] += (sigma*mcs_i*sigma);
                model.res[2] += (sigma*mcs_i*sigma)*(sigma*sigma);
                model.res[3] += HH*mcs_i;
                model.res[4] += HH*mcs_i*HH;
                
                modelW2.WriteLine(",");
                modelW2.WriteLine(sigma*kNi);
            }
            modelW2.WriteLine("\n");

            model.MV[i] = model.res[0]/kN;
            model.CV[i] = (model.BetaV[i]*model.BetaV[i])/kN*(model.res[4]-model.res[3]*model.res[3]);        

            cout << model.MV[i] << "\t" << model.CV[i] << "\t" << model.Fliped_Step << right << setw(10) << "\t" << model.Total_Step << endl;

            string temp = to_string(i) + "," + to_string(model.TV[i]) + "," + to_string(model.MV[i]) + "," + to_string(model.CV[i]) + ",";
            temp = temp + to_string(model.res[0]) + "," + to_string(model.res[1]) + "," + to_string(model.res[2]) + ",";
            temp = temp + to_string(model.res[3]) + "," + to_string(model.res[4]) + "\n";
            modelW.WriteLine(temp);
        }
        modelW.CloseNewFile();
        modelW2.CloseNewFile();
    }
    Farewell();
}