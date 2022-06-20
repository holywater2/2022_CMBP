#ifndef ____DebugHelper____
#define ____DebugHelper____

#include <iostream>
#include <string>

using namespace std;

string ToStringSc(short* array, int L){
    string res = "";
    for(int i = 0; i < L; i++){
        for(int j = 0; j < L; j++){
            if(array[i*L+j] == 1) res += "■";
            else res += "□";	
        }
        res += "\n";
    }
    return res;   
}

#endif // ____DebugHelper____