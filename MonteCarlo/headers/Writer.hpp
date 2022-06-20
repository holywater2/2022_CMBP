#ifndef ____Writer____
#define ____Writer____

// File Stream and etc. 
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string>

// Error catch
#include <signal.h>
#include <cstdlib>

#ifdef _WIN32
#include <windows.h>
#endif // _WIN32
#ifdef linux
#include <sys/stat.h>
#endif // linux




using namespace std;

class Writer {
    public:
        ofstream myfile;
        string file_name;
        string new_file_name;
        string directory_name;
        void MakeDirectory();
        // Wirter(string filename);
        // ~Writer();
        Writer(string filename, string directory_name);
        // Writer(bool Foo){/*Intentialnally empty for FooWriter*/};
        void FindNextFileNum();
        void OpenNewFile();
        void WriteLine(string contents);
        void WriteLine(double contents);
        void WriteLine(int contents);
        void CloseNewFile();;
};

// Writer::Writer(string filename){
//     this-> new_file_name = filename;
//     this-> directory_name = "Result";
// }

Writer::Writer(string filename, string directory_name = "Result"){
    this-> file_name = filename;
    this-> directory_name = directory_name;

    this-> MakeDirectory();
    this-> FindNextFileNum();
    this-> OpenNewFile();
    cout << "Saving Start: " << this-> new_file_name << "\n";
}
// Writer::~Writer(){
//     this-> CloseNewFile();
// }

void Writer::MakeDirectory(){
    #ifdef _WIN32
    CreateDirectory(this->directory_name.c_str(), NULL);
    #endif // _WIN32

    #ifdef linux
    mkdir(this->directory_name.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    #endif // linux
}

void Writer::FindNextFileNum(){
    ifstream f(this-> new_file_name);
    int filenum = 1;
    while(true){
        this-> new_file_name = this->file_name + "_" + to_string(filenum++) + ".csv";
        ifstream f(this-> new_file_name);
        if(!f.good()) break;
    }
    f.close();
}

void Writer::OpenNewFile(){
    this->myfile.open(this->new_file_name);
}

void Writer::WriteLine(string contents){
    // This assumes that there is no line break.
    myfile << contents;
}

void Writer::WriteLine(double contents){
    // This assumes that there is no line break.
    myfile << contents;
}

void Writer::WriteLine(int contents){
    // This assumes that there is no line break.
    myfile << contents;
}

void Writer::CloseNewFile(){
    this-> myfile.close();
    cout << "Save Completed: " << this-> new_file_name << "\n";
}

class FooWriter : Writer {
    public:
        ofstream myfile;
        string file_name;
        string new_file_name;
        string directory_name;
        void MakeDirectory(){}
        // FooWriter(string filename, string directory_name = "Result") : Writer(true) {}
        FooWriter(string filename, string directory_name = "Result") : Writer(filename,directory_name) {}
        void FindNextFileNum(){}
        void OpenNewFile(){}
        void WriteLine(string contents){}
        void WriteLine(double contents){}
        void WriteLine(int contents){}
        void CloseNewFile(){}
};

#endif // ____Writer____
