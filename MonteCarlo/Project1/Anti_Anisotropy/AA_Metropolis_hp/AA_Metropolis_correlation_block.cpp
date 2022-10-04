#include "AA_Metropolis_program_header.hpp"

// arguments list that helps to pass the args to model
vector<string> result_to_file = vector<string>();
vector<string> result_to_file_cor = vector<string>();
vector<string> result_to_file_cor_err = vector<string>();

int main(int argn, char *argv[]){ // Input argument: argv[0]--> file name / argv[1]--> Input parameter
    signal(SIGSEGV, &handler);
    signal(SIGINT, &handler);
    if(argn >= 2){
        string Input_file = argv[1];
        vector<FLOAT1> input = Writer::Argument_readerL(Input_file,13);
        kLx = (INT1)input[0]; kLy = (INT1)input[1]; kN = kLx*kLy;
        kBin = (INT1)input[2]; kB = input[3]; kJx = input[4]; kJy = input[5];
        alpha = input[6]; Tsrt = input[7]; Tfin = input[8];
        isTinf = input[9]; Random = input[10];
        equil_time_base = input[11]; mcs = input[12];
    }
    // arguments list that helps to pass the args to model
    vector<FLOAT1> args = {kLx,kLy,kBin,kB,kJx,kJy,alpha,Tsrt,Tfin,isTinf,Random,equil_time_base,mcs};

    // Filename Base: '\Result\(Model Name)_c_(kL)_int[erval]_(kBin) + (blahblah)
    #ifdef _WIN32
    static string kFilename = ".\\Result\\"+Model::Name()+"_c_"+to_string(kLx)+"_"+to_string(kLy)+"_int"+to_string(kBin)+"_mcs"+to_string(mcs)+"_a"+to_string(alpha);
    #endif
    #ifdef linux
    static string kFilename = "./Result/"+Model::Name()+"_c_"+to_string(kLx)+"_"+to_string(kLy)+"_int"+to_string(kBin)+"_mcs"+to_string(mcs)+"_a"+to_string(alpha);
    #endif

    Greetings();

    Model model = Model(args);

    FLOAT1 MM, HH;
    FLOAT1 mcs_i = 1/FLOAT1(mcs);
    FLOAT1 kNi   = 1/FLOAT1(kN);
    FLOAT1 kNi2  = kNi*kNi;
    FLOAT1 kNir  = 1/pow(kN,0.5);

    for(int cBin = 0; cBin < kBin; cBin++){
        model.Initialize(model.BetaV[cBin]);

        // Output data
        // 0: <m>, 1: <m^2>, 2: <m^4>, 3: <E>/sqrt(N), 4: <E^2>/N
        model.res = vector<FLOAT1>(5,0);
        vector<INT2>   cor_sum(kN);
        vector<FLOAT2> cor_avg(kN);

        equil_time = equil_time_base;
        // if(model.TV[i]<=2.4 || model.TV[i]>=2.0) equil_time =1000;

        model.IterateUntilEquilibrium(equil_time);

        model.Measure_fast();
        HH = model.HH;
        MM = model.staggered;

        cout <<"idx: " << left << setw(4) << cBin << "|| " << left << setw(10) << model.TV[cBin];
        cout << "|| "  << left << setw(9) << MM/(FLOAT1)kN << "  " << left << setw(12) << HH;
        cout << "|| "  << left << setw(9) << model.sigma/(FLOAT1)kN << "|| ";
        cout << '\n';

        /***********Monte Carlo Step and Caculate the data***********/
        int block_size = 1000;
        int blocks = mcs/block_size;
        FLOAT1 bsi = 1/(long double) block_size;

        vector<FLOAT2> Block_MM(blocks);
        vector<FLOAT2> Block_MM2(blocks);
        vector<FLOAT2> Block_MM4(blocks);
        vector<FLOAT2> Block_HH(blocks);
        vector<FLOAT2> Block_HH2(blocks);
        vector<FLOAT2> cor_err(kN);
        vector<INT2>   Block_cor_sum(blocks*kN);

        for(int bidx = 0; bidx < blocks; bidx++){
            INT2 blocksum_MM    = 0;
            INT2 blocksum_MM2   = 0;
            INT2 blocksum_MM4   = 0;
            FLOAT2 blocksum_HH  = 0;
            FLOAT2 blocksum_HH2 = 0;

            for(int k = 0; k < block_size; k++){
                model.Calculate();
                model.Measure_fast();

                HH = model.HH;
                MM = abs(model.staggered);

                blocksum_MM   += MM;
                blocksum_MM2  += MM*MM;
                blocksum_MM4  += MM*MM*MM*MM;
                blocksum_HH   += HH;
                blocksum_HH2  += HH*HH;

                model.CalculateCorrelation();
                // model.CalculateCorrelationStaggered();
                for(int ii = 0; ii < kN; ii++)
                    Block_cor_sum[bidx*kN + ii] += model.cor[ii];
                    // cor_sum[ii] += model.cor_stag[ii];
            }
            Block_MM[bidx]  = blocksum_MM *bsi*kNi;
            Block_MM2[bidx] = blocksum_MM2*bsi*kNi2;
            Block_MM4[bidx] = blocksum_MM4*bsi*kNi2*kNi2;
            Block_HH[bidx]  = blocksum_HH *bsi*kNir;
            Block_HH2[bidx] = blocksum_HH2*bsi*kNi;

            for(int ii = 0; ii < kN; ii++)
                cor_sum[ii] += Block_cor_sum[bidx*kN + ii];

            model.res[0] += Block_MM[bidx];          // = <m>
            model.res[1] += Block_MM2[bidx];         // = <m^2>
            model.res[2] += Block_MM4[bidx];         // = <m^4>
            model.res[3] += Block_HH[bidx];          // = <E>/sqrt(N)
            model.res[4] += Block_HH2[bidx];         // = <E^2>/N
        }

        for(int j = 0; j < 5; j++)
            model.res[j] /= (FLOAT1) blocks;  // average of blocks

        /***********************************************************/

        /*******Calculate Magnetizaition and Specific Heat.*********/
        model.MV[cBin] = model.res[0];
        model.CV[cBin] = (model.BetaV[cBin]*model.BetaV[cBin])*(model.res[4]-model.res[3]*model.res[3]);
        model.BV[cBin] = 0.5*(3-model.res[2]/(model.res[1]*model.res[1]));
        /***********************************************************/

        /*******Calculate Correlation.*********/
        for(int j = 0; j < kN; j++)
            cor_avg[j] = (FLOAT2)(2*cor_sum[j]-mcs)/mcs;
        /***********************************************************/

        /*******Calculate Block Error.*********/
        FLOAT2 MMerr = 0;
        FLOAT2 CCerr = 0;
        FLOAT2 BBerr = 0;
        for(int bidx = 0; bidx < blocks; bidx++){
            MMerr += (Block_MM[bidx]-model.MV[cBin])*(Block_MM[bidx]-model.MV[cBin]);

            FLOAT2 CCtemp = model.BetaV[cBin]*model.BetaV[cBin]*(Block_HH2[bidx]-Block_HH[bidx]*Block_HH[bidx]);
            CCerr += (CCtemp-model.CV[cBin])*(CCtemp-model.CV[cBin]);

            FLOAT2 BBtemp = 0.5*(3-Block_MM4[bidx]/(Block_MM2[bidx]*Block_MM2[bidx]));
            BBerr += (BBtemp-model.BV[cBin])*(BBtemp-model.BV[cBin]);
            for(int ii = 0; ii < kN; ii++){
                FLOAT2 cor_temp = (FLOAT2)(2*Block_cor_sum[bidx*kN + ii]-block_size)/block_size;
                cor_err[ii] += (cor_avg[ii]-cor_temp) * (cor_avg[ii]-cor_temp);
            }
        }
        MMerr = pow(MMerr/(blocks-1),0.5);
        CCerr = pow(CCerr/(blocks-1),0.5);
        BBerr = pow(BBerr/(blocks-1),0.5);
        for(int ii = 0; ii < kN; ii++){
            cor_err[ii] = pow(cor_err[ii]/(blocks-1),0.5);
        }
        /***********************************************************/


        cout <<"                         ";
        cout << left << setw(8) << model.MV[cBin] << "  " << left << setw(12) << model.CV[cBin] << "|| ";
        cout << left << setw(14) << model.Fliped_Step << "  " << left << setw(10) << model.Total_Step << endl;

        // Standard result save
        string result = to_string(cBin) + "," + to_string(model.TV[cBin]) + "," + to_string(model.MV[cBin]) + "," + to_string(model.CV[cBin]) + ",";
        result = result + to_string(model.res[0]) + "," + to_string(model.res[1]) + "," + to_string(model.res[2]) + ",";
        result = result + to_string(model.res[3]) + "," + to_string(model.res[4]) + "," + to_string(model.BV[cBin]) + ",";
        result = result + to_string(MMerr) + "," + to_string(CCerr) + "," + to_string(BBerr) + "\n";

        result_to_file.push_back(result);

        // correlation result save
        string result2 = to_string(cBin) + "," + to_string(model.TV[cBin]) + ",";
        for(int j = 0; j < kN; j++){
            result2 = result2 + to_string(cor_avg[j]) + ",";
        }
        result2.pop_back();
        result_to_file_cor.push_back(result2 + "\n");

        // correlation err result save
        string result3 = to_string(cBin) + "," + to_string(model.TV[cBin]) + ",";
        for(int j = 0; j < kN; j++){
            result3 = result3 + to_string(cor_err[j]) + ",";
        }
        result3.pop_back();
        result_to_file_cor_err.push_back(result3 + "\n");
    }

    kFilename += "_corBlock";
    /***********Save the result of the Calculation**********/
    Writer modelW = Writer(kFilename);
    modelW.WriteLine("idx,temperture,(staggered)magnetization,specific heat,abs(mm),mm**2,mm**4,HH/L,HH**2/L,Binder,MMerr,CCerr,BBerr\n");
    for(int i = 0; i < kBin; i++)
        modelW.WriteLine(result_to_file.at(i));
    modelW.CloseNewFile();
    /******************************************************/

    string cor_idx = "idx,temperture,cor\n";
    for(int i = 0; i < kN; i++){
        cor_idx += to_string(i) + ',';
    }
    cor_idx.pop_back();
    cor_idx += '\n';

    /***********Save the result of the Correlation**********/
    Writer modelW2 = Writer(kFilename+"_corres");
    modelW2.WriteLine(cor_idx);
    for(int i = 0; i < kBin; i++)
        modelW2.WriteLine(result_to_file_cor.at(i));
    modelW2.CloseNewFile();
    /******************************************************/

    /***********Save the result of the Correlation**********/
    Writer modelW3 = Writer(kFilename+"_corerr");
    modelW3.WriteLine(cor_idx);
    for(int i = 0; i < kBin; i++)
        modelW3.WriteLine(result_to_file_cor_err.at(i));
    modelW3.CloseNewFile();
    /******************************************************/
    Farewell();
}