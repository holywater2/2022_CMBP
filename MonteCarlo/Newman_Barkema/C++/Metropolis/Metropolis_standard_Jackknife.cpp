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

const double isTinf = false;
const bool random = true;

#ifdef _WIN32
static string kFilename = ".\\Result\\Metropolis_c_"+to_string(kL)+"_int"+to_string(kBin);
#endif
#ifdef linux
static string kFilename = "./Result/Metropolis_c_"+to_string(kL)+"_int"+to_string(kBin);
#endif

typedef Metropolis_2D Model;

vector<double> args = {kL,kBin,kB,kJ,Tsrt,Tfin,isTinf};

clock_t __start__, __finish__;

void Greetings(){
    string Tat = isTinf ? "inf" : "0";
    string isRandom = random ? "Random" : "Chess board";

    cout << "Metropolis Algorithm\n";
    cout << "Radnomness test(seed): " << seed << '\n';
    cout << "L = " << kL << ", " << "bin = " << kBin << ", Start T at " << Tat <<", selecting sites: " <<  isRandom <<"\n";
    cout << "------------------------------------------------------------------------------------------------------------------" << "\n";
    cout << "--index--||---Temp----||EQ:sig------HH----------||magnetization---specific heat||Fliped Step------Total Step------" << "\n";
    cout << "------------------------------------------------------------------------------------------------------------------" << endl;
    cout << fixed <<setprecision(6);
    cout << showpos;

    __start__ = clock();
}

void Farewell(){
    __finish__ = clock();
    cout << "Program Exit. Spent time: " << (double)(__finish__-__start__)/CLOCKS_PER_SEC << "\n";
    cout << "-------------------------------------------------------------------------------------------";
}

int main(){
    Greetings();
    
    Model model = Model(args);
    Writer modelW = Writer(kFilename);
    modelW.WriteLine("idx,temperture,magnetization,specific heat,abs(sigma),sigma**2,sigma**4,HH,HH**2,c_error\n");

    int equil_time, mcs = 10e6;
    double MM, HH;
    double mcs_i = 1/double(mcs);
    double kNi = 1/double(kN);

    for(int gg = 0; gg < 1; gg++){
        for(int i = 0; i < kBin; i++){
            model.Initialize(model.BetaV[i]);
            model.res = vector<double>(5,0);

            equil_time = 2000;

            if(model.TV[i]<=2.2 || model.TV[i]>=2.0) equil_time =10000;

            model.IterateUntilEquilibrium(equil_time, random);

            duo value = model.Measure();
            HH = get<0>(value);
            MM =  get<1>(value);

            cout <<"idx: " << left << setw(4) << i << "|| " << left << setw(10) << model.TV[i];
            cout << "|| "  << left << setw(9) << MM/(double)kN << "  " << left << setw(12) << HH << "|| ";

            int jB = 100;
            vector<double> Jackknife = vector<double>(mcs/jB,0);
            vector<double> Jackknife_HH = vector<double>(mcs/jB,0);
            vector<double> Jackknife_HH2 = vector<double>(mcs/jB,0);
            double c_error = 0;

            for(int j = 0; j < mcs; j++){ // mcs 15000
                model.Calculate(0,random);

                value = model.Measure_fast();

                HH = get<0>(value)/(double) kL;        // = E
                MM = abs(get<1>(value))/(double)kN;    // = M
                model.res[0] += MM*mcs_i;              // = <m>
                model.res[1] += (MM*mcs_i*MM);         // = <m^2>
                model.res[2] += (MM*mcs_i*MM)*(MM*MM); // = <m^4>
                model.res[3] += HH*mcs_i;              // = <E>/sqrt(N)
                model.res[4] += HH*mcs_i*HH;           // = <E^2>/N

                Jackknife_HH[j/jB] += HH;
                Jackknife_HH2[j/jB] += HH*HH;
            }
            
            model.MV[i] = model.res[0]/kN;
            model.CV[i] = (model.BetaV[i]*model.BetaV[i])*(model.res[4]-model.res[3]*model.res[3]);        

            for(int j = 0; j < mcs/jB; j++){
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