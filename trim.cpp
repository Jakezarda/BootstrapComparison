#include<iostream>
#include<iomanip>
#include<fstream>
#include<sstream>
#include<vector>
#include<string>

std::string filename(std::string base, int digits, int num, std::string ext){
    std::stringstream file;
    file << base << std::setw(digits) << std::setfill('0') << num << "." << ext;       //creates files with names such that filename0001.extension
    return file.str();
}

int main(){
    int N_s = 1000;
    int N_trim = 24;                                                                      //came to 23 in class, from own data 0's stop after 16 and before 150, N_trim = 24?
    
    for (int sample = 0; sample < N_s; ++sample){
        std::string file = filename("Pis_hist-" , 4, sample + 1, "dat");
        
        std::vector<double> bins;
        std::vector<int> hist;
        
        
        std::ifstream fin(file);
        while(!fin.eof()){
            double bin;
            int val;
            fin >> bin >> val;
            if (!fin.eof()){
                bins.push_back(bin);
                hist.push_back(val);
            }
            
        }
        fin.close();
        
        std::ofstream fout(file);
        fout.precision(15);
        for (int i = N_trim; i < bins.size() - N_trim; ++i){
            fout <<bins[i] << " " << hist[i] << "\n";
            
        }
        fout.close();
    }
    std::vector<double> bins;
    std::vector<int> hist;
    
    std::ifstream fin("Pis_hist.dat");
    while(!fin.eof()){
        double bin;
        int val;
        fin >> bin >> val;
        if (!fin.eof()){
            bins.push_back(bin);
            hist.push_back(val);
        }
    }
    fin.close();
    
    std::ofstream fout("Pis_hist.dat");
    fout.precision(15);
    for (int i = N_trim; i < bins.size() - N_trim; ++i){
        fout << bins[i] << " " << hist[i] << "\n";
    }
    fout.close();
         
}






















