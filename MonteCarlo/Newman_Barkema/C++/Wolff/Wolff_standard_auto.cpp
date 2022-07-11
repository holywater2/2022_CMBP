#include "Wolff.hpp"
#include "../Writer.hpp"

#include <iostream>
#include <iomanip>

const int kL = 100; /*Parameter: lattice size*/
const int kN = kL*kL;
const int kBin = 25; /*Parametr: Change binning of temperature*/
const int kB = 0;
const int kJ = 1;

const double Tsrt = 2.2;
const double Tfin = 2.4;

double isTinf = false;

#ifdef _WIN32
static string kFilename = ".\\Result\\Wolff_c_"+to_string(kL)+"_int"+to_string(kBin);
#endif //_WIN32
#ifdef linux
static string kFilename = "./Result/Wolff_c_"+to_string(kL)+"_int"+to_string(kBin);
#endif //linux

typedef Wolff_2D Model;

vector<double> args = {kL,kBin,kB,kJ,Tsrt,Tfin,isTinf};

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
    Farewell();
    exit(A);
}

int main(){
    signal(SIGSEGV, &handler);
    signal(SIGINT, &handler);
    Greetings();
    
    for(int g = 0; g < 1; g++){
        Model model = Model(args);
        Writer modelW = Writer(kFilename + "_auto3");
        Writer modelW2 = Writer(kFilename + "_auto4");

        modelW.WriteLine("idx,temperture,magnetization,specific heat,abs(sigma),sigma**2,sigma**4,HH,HH**2,m_error\n");

        int HH, equil_time, mcs = 1e5;
        double sigma;
        double mcs_i = 1/double(mcs);
        double kNi = 1/double(kN);

        cout << showpos;

        for(int i = 0; i < kBin; i++){
            model.Initialize(model.BetaV[i]);
            model.res = vector<double>(5,0);
    
            // if(model.TV[i]<2.4 || model.TV[i]>2.0) equil_time =3000;
            // model.IterateUntilEquilibrium(2000);

            duo value = model.Measure();
            HH = get<0>(value);
            sigma =  get<1>(value);

            cout <<"idx: " << i << " \t|| " << model.TV[i] << "\t|| " << sigma << "\t" << HH << setw(3) << "|| ";
            modelW2.WriteLine(model.TV[i]);
           
            double size, step;
            for(double j = 0; j < mcs;){ // mcs 15000
                size = model.Calculate();
                step = size/(double)kN;
                j += step;
                step /= mcs;

                value = model.Measure();
                HH = get<0>(value);
                sigma =  get<1>(value)/(double)kN;
                
                model.res[0] += abs(sigma)*step;
                model.res[1] += (sigma*step*sigma);
                model.res[2] += (sigma*step*sigma)*(sigma*sigma);
                model.res[3] += HH*step/(double)kL;
                model.res[4] += HH*step*HH/(double)kN;

                modelW2.WriteLine(",");
                modelW2.WriteLine(sigma);
            }
            modelW2.WriteLine(",");
            modelW2.WriteLine(model.Fliped_Step/(double)kN/model.Total_Step);
            modelW2.WriteLine("\n");

            model.MV[i] = model.res[0];
            model.CV[i] = (model.BetaV[i]*model.BetaV[i])*(model.res[4]-model.res[3]*model.res[3]);        


            cout << model.MV[i] << "\t" << model.CV[i] << "\t" << model.Fliped_Step << right << setw(10) << "\t" << model.Total_Step << "\t" << (model.Fliped_Step/(double)kN/model.Total_Step) << endl;
            
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
