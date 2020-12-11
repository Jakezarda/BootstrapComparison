//can be applied to many codes, this one will be more tuned to the pi problem
#include<iostream>
#include<iomanip>
#include<fstream>
#include<sstream>
#include<vector>
#include<string>
#include<cmath>

std::string filename(std::string base, int digits, int num, std::string ext){
    std::stringstream file;
    file << base << std::setw(digits) << std::setfill('0') << num << "." << ext;       //creates files with names such that filename0001.extension
    return file.str();
}


int main(){
    int N = 0;
    int N_s = 1000;
    std::ifstream fin("Pk_LNKNLogs_0001.dat");
    //std::ofstream cov_err1("Pis_cov_err1.dat");
    //cov_err1.precision(15);
    while(!fin.eof()){
        double bin;
        int val;
        fin >> bin >> val;
    
        if (!fin.eof()){
            N += 1;   
        }
        
    }
    N /= 2;
    int N_dof = (N/2)*(N+1);
      
    fin.close();
    std::cout << "initial file read in" << std::endl;
    std::vector<std::vector<double>> C(N,std::vector<double>(N));
    std::vector<double> avg(N);
    std::vector<std::vector<double>> C_var(N,std::vector<double>(N));
    
    for (int sample = 0; sample < N_s; ++sample){
        std::string file = filename("Pk_LNKNLogs_", 4, sample + 1, "dat");
        
        fin.open(file);
        for (int i = 0; i < N; ++i){
            double bin, val;
            fin >> bin >> val;
            avg[i] += val/double(N_s);
        }   
        fin.close();
    }
    std::cout << "averages found in first file read for loop" << std::endl;
    
    double var = 0;
    double var2 = 0;
    for (int sample = 0; sample < N_s; ++sample){
        std::string file = filename("Pk_LNKNLogs_", 4, sample + 1, "dat");
        std::vector<int> hist(N);
        
        fin.open(file);
        for (int i = 0; i < N; ++i){
            double bin, val;
            fin >> bin >> val;
            hist[i] = val;
            
        }
        fin.close();
        
        for (int i = 0; i < N; ++i){
            for (int j = 0; j < N; ++j){
                C[i][j] += ((hist[i]-avg[i])*(hist[j]-avg[j]))/(N_s-1.0);
                C_var[i][j] += pow((((hist[i]-avg[i])*(hist[j]-avg[j]))-C[i][j]),2)/(N_s-1.0);
            }
        }
        //std::cout << ((hist[36]-avg[36])*(hist[36]-avg[36]))/(999) << std::endl;
        //var += (hist[0] - avg[0])*(hist[0]-avg[0]);
        //var2 += (hist[0] - avg[0])*(hist[1]-avg[1]);
        //std::cout << avg[0] << std::endl;
    }

    
    std::cout << "C[i][j] matrix made" << std::endl;
    //std::cout << C[36][36] << std::endl;
    
    std::ofstream cov_out("Pk_mock_cov.dat");
    std::ofstream cor_out("Pk_mock_cor.dat");
    //std::ofstream stdev_out("AAAA.dat");
    cov_out.precision(15);
    cor_out.precision(15);
    //stdev_out.precision(15);
    
    for (int i = 0; i < N; ++i){
        for (int j = 0; j < N; ++j){
            cov_out.width(25);
            cov_out << C[i][j];
            cor_out.width(25);
            cor_out << C[i][j]/std::sqrt(C[i][i]*C[j][j]);
            //stdev_out.width(25);
            //stdev_out << C_stdev[i][j];
        }
        cov_out << "\n";
        cor_out << "\n";
    }
    cov_out.precision(15);
    cor_out.precision(15);
    cov_out.close();
    cor_out.close();
    
    //This section is to attempt to read in the bootstrapped version and compute its covariance in the same function as the mock spectrum.
    
    std::vector<std::vector<double>> C_boot(N,std::vector<double>(N));
    std::vector<double> avg_boot(N);
    
        for (int sample = 0; sample < N_s; ++sample){
            std::string file_boot = filename("Pk_bootstrap", 4, sample + 1, "dat");
        
        fin.open(file_boot);
        for (int i = 0; i < N; ++i){
            double bin, val;
            fin >> bin >> val;
            avg_boot[i] += val/double(N_s);
        }   
        fin.close();
    }
    std::cout << "secondary file read in" << std::endl;
        for (int sample = 0; sample < N_s; ++sample){
        std::string file_boot = filename("Pk_bootstrap", 4, sample + 1, "dat");
        std::vector<int> hist_boot(N);
        
        fin.open(file_boot);
        for (int i = 0; i < N; ++i){
            double bin, val;
            fin >> bin >> val;
            hist_boot[i] = val;
            
        }
        fin.close();
        
        for (int i = 0; i < N; ++i){
            for (int j = 0; j < N; ++j){
                C_boot[i][j] += ((hist_boot[i]-avg_boot[i])*(hist_boot[j]-avg_boot[j]))/(N_s-1.0);
            }
        }
    
    }
    
    std::ofstream cov_out_boot("Pk_boot_cov.dat");
    std::ofstream cor_out_boot("Pk_boot_cor.dat");
    //std::ofstream stdev_out("AAAA.dat");
    cov_out_boot.precision(15);
    cor_out_boot.precision(15);
    //stdev_out.precision(15);
    
    for (int i = 0; i < N; ++i){
        for (int j = 0; j < N; ++j){
            cov_out_boot.width(25);
            cov_out_boot << C_boot[i][j];
            cor_out_boot.width(25);
            cor_out_boot << C_boot[i][j]/std::sqrt(C_boot[i][i]*C_boot[j][j]);
            
        }
        cov_out_boot << "\n";
        cor_out_boot << "\n";
    }
    cov_out_boot.precision(15);
    cor_out_boot.precision(15);
    cov_out_boot.close();
    cor_out_boot.close();
    
    double chi_squared = 0;
    
    for (int i = 0; i < N; ++i){
        for (int j = i; j < N; ++j) {
            chi_squared += (C_boot[i][j] - C[i][j])/(C_var[i][j]);
        }
        
    }
    
    chi_squared /= N_dof;
    
    std::cout << chi_squared << std::endl;
    
    return 0;
}
//to change directory in gnu plot; cd "directory path"




