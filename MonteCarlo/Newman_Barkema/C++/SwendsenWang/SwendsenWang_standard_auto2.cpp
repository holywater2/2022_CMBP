#include "SwendsenWang_HK.hpp"
#include "../Writer.hpp"
#include "../DebugHelper.hpp"

#include <iostream>
#include <iomanip>

// const int kL   = 10;        /*Parameter: lattice size*/
// const int kN   = kL*kL;
const int kBin = 25;        /*Parameter: Change binning of temperature*/
const int kB   = 0;
const int kJ   = 1;

// T_crit ~ 2.269
// const double Tsrt = T_CRIT;
// const double Tfin = T_CRIT;
const double Tsrt = 2.2;
const double Tfin = 2.35;

double isTinf = false;

typedef SwendsenWang_2D Model;

// vector<double> args = {kL,kBin,kB,kJ,Tsrt,Tfin,isTinf};

clock_t __start__, __finish__;

void Greetings(int kL){
    string Tat = isTinf ? "inf" : "0";

    cout << Model::Name() + " Algorithm\n";
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
    if(N)
        cout << "\nProgram Abonormally Exit. Spent time: " << (double)(__finish__-__start__)/CLOCKS_PER_SEC << "\n";
    else
        cout << "Program Exit. Spent time: " << (double)(__finish__-__start__)/CLOCKS_PER_SEC << "\n";
    cout << "-------------------------------------------------------------------------------------------\n";
}

void handler(int A){
    cout << endl;
    Farewell(1);
    exit(A);
}

int main(){
    signal(SIGSEGV, &handler);
    signal(SIGINT, &handler);
    
    // vector<int> kLL = {4,8,16,32,64,128,256,512};
    vector<int> kLL = {24,32,48,64,96};

    for(int gg = 0; gg < kLL.size(); gg++){
        double kL = kLL[gg];
        double kN = kL*kL;
    
        Greetings(kL);

        #ifdef _WIN32
        static string kFilename = ".\\Result\\"+Model::Name()+"_c_"+to_string(kL)+"_int"+to_string(kBin);
        #endif //_WIN32
        #ifdef linux
        static string kFilename = "./Result/"+Model::Name()+"_c_"+to_string(kL)+"_int"+to_string(kBin);
        #endif //linux

        vector<double> args = {kL,kBin,kB,kJ,Tsrt,Tfin,isTinf};

        Model model = Model(args);
        Writer modelW = Writer(kFilename+"1e5_auto3");
        // Writer modelW2 = Writer(kFilename + "1e5_auto4");
        modelW.WriteLine("idx,temperture,magnetization,specific heat,abs(MM),MM**2,MM**4,HH,HH**2,m_error\n");
        
        /* Parameter */
        int equil_time = 100;
        int mcs = 1e5;
        /*************/

        double MM, HH;
        double mcs_i = 1/double(mcs);
        double kNi = 1/double(kN);

        cout << showpos;

        for(int i = 0; i < kBin; i++){
            model.Initialize(model.BetaV[i]);
            model.res = vector<double>(5,0);

            // equil_time = 50;
            // if(model.TV[i]<2.4 || model.TV[i]>2.0) equil_time = 100;

            model.IterateUntilEquilibrium(equil_time);

            duo value = model.Measure();
            HH = get<0>(value);
            MM = get<1>(value);

            cout <<"idx: " << left << setw(4) << i << "|| " << left << setw(10) << model.TV[i];
            cout << "|| "  << left << setw(9) << MM/(double)kN << "  " << left << setw(12) << HH << "|| ";
            // modelW2.WriteLine(model.TV[i]);

            double size;
            for(int j = 0; j < mcs; j++){
                size = model.Calculate();
                
                value = model.Measure();  
                HH = get<0>(value)/(double)kL;          // = E
                MM = get<1>(value)/(double)kN;          // = M
                
                model.res[0] += abs(MM)*mcs_i;               // = <m>
                model.res[1] += (MM*mcs_i*MM);          // = <m^2>
                model.res[2] += (MM*mcs_i*MM)*(MM*MM);  // = <m^4>
                model.res[3] += HH*mcs_i;               // = <E>/sqrt(N)
                model.res[4] += HH*mcs_i*HH;            // = <E^2>/N

                // modelW2.WriteLine(",");
                // modelW2.WriteLine(MM);
            }
            // modelW2.WriteLine(",");
            // modelW2.WriteLine(1);
            // modelW2.WriteLine("\n");

            model.MV[i] = model.res[0];
            model.CV[i] = (model.BetaV[i]*model.BetaV[i])*(model.res[4]-model.res[3]*model.res[3]);        

            cout << left << setw(13) << model.MV[i] << "  " << right << setw(13) << model.CV[i] << "|| ";
            cout << left << setw(14) << model.Fliped_Step << "  " << left << setw(10) << model.Total_Step << endl;

            // cout << ToStringSc(model.sc, kL);

            string temp = to_string(i) + "," + to_string(model.TV[i]) + "," + to_string(model.MV[i]) + "," + to_string(model.CV[i]) + ",";
            temp = temp + to_string(model.res[0]) + "," + to_string(model.res[1]) + "," + to_string(model.res[2]) + ",";
            temp = temp + to_string(model.res[3]) + "," + to_string(model.res[4]) + "\n";
            modelW.WriteLine(temp);
        }
        modelW.CloseNewFile();
        // modelW2.CloseNewFile();
    }
    Farewell();
}
