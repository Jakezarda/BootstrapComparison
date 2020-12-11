#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <fftw3.h>

#include "binner.hpp"
#include "galaxy.hpp"

std::string filename(std::string base, int digits, int num, std::string ext){
    std::stringstream file;
    file << base << std::setw(digits) << std::setfill('0') << num << "." << ext;       //creates files with names such that filename0001.extension
    return file.str();
}


/*std::vector<galaxy> readRandoms(std::string filename) {
    std::vector<galaxy> randoms;
    std::ifstream fin(filename);
    
    while (!fin.eof()) {
        std::vector<double> ranPos(6);
      
        
        fin >> ranPos[0] >> ranPos[1] >> ranPos[2] >> ranPos[3] >> ranPos[4] >> ranPos[5];
        if (!fin.eof()) {
            galaxy rans(ranPos);      
            randoms.push_back(rans);
        }
    }
    fin.close();
    
    return randoms;
}*/

std::vector<galaxy> fileRead(std::string filename) {
    std::vector<galaxy> gals;
    std::ifstream fin(filename);
    
    while (!fin.eof()) {
        std::vector<double> pos(3), vel(3), pos2(3);
        vel[0] = 0;
        vel[1] = 0;
        vel[2] = 0;
        pos2[0] = 0;
        pos2[1] = 0;
        pos2[2] = 0;
        
        fin >> pos[0] >> pos[1] >> pos[2] >> pos2[0] >> pos2[1] >> pos2[2];
        if (!fin.eof()) {
            galaxy gal(pos, vel, 1.0, 0.0004, "Cartesian");      
            gals.push_back(gal);
        }
    }
    fin.close();
    
    return gals;
}

std::vector<double> fftFreq(int N, double L) {
    std::vector<double> k;
    double dk = 2.0*M_PI/L;
    for (int i = 0; i <= N/2; ++i)
        k.push_back(i*dk);
    for (int i = N/2 + 1; i < N; ++i)
        k.push_back((i - N)*dk);
    return k;
}

int main() {
    //std::cout << "no memory allocated yet" << std::endl;
    std::vector<int> N_grid = {512, 512, 512};
    std::vector<double> L = {1024.0, 1024.0, 1024.0};
    double nbar = 0.0004;
    int N_mocks = 1000;                                            //temporarily changed to 1, will have to make 1000 by the end of this
    
    //std::cout << "creating galaxy arrays" << std::endl;
    //std::vector<galaxy> gals = fileRead("LNKNLogsVel_01.dat");
    
    std::vector<galaxy> rans = fileRead("LNKNLogsVel_Random.dat");
    
    for (int sample = 0; sample < N_mocks; ++sample){
      //  std::cout << "entering for loop" << std::endl;
        std::string inputfile = filename("LNKNLogs_" , 4, sample + 1, "dat");
        std::vector<galaxy> gals = fileRead(inputfile);
        //std::cout << "mock catalog read" << std::endl;
        
        binner binGals;
        std::vector<double> n_gals(N_grid[0]*N_grid[1]*N_grid[2]);
        std::vector<double> n_rans(N_grid[0]*N_grid[1]*N_grid[2]);
        std::vector<double> nbw_gal(3);
        std::vector<double> nbw_ran(3);
        
        std::vector<fftw_complex> dk(N_grid[0]*N_grid[1]*(N_grid[2]/2 + 1));
        //std::cout << "arrays created" << std::endl;
        
        //std::cout << "entering while statement" << std::endl;
        std::ifstream fin(inputfile);
        while(!fin.eof()){
            double n_gal;
            double n_ran;
            double nbw_gal_val;
            double nbw_ran_val;
            double dk_val;
            //fftw_complex dk_val;
            fin >> n_gal >> n_ran >> nbw_gal_val >> nbw_ran_val >> dk_val;
            if (!fin.eof()){
                n_gals.push_back(n_gal);
                n_rans.push_back(n_ran);
                nbw_gal.push_back(nbw_gal_val);
                nbw_ran.push_back(nbw_ran_val);
            }
        }
        fin.close();
        
        //std::cout << "left while statement" << std::endl;
        
        fftw_import_wisdom_from_filename("wisdom");
        fftw_plan dr2dk = fftw_plan_dft_r2c_3d(N_grid[0], N_grid[1], N_grid[2], n_gals.data(), dk.data(), FFTW_ESTIMATE);
        fftw_export_wisdom_to_filename("wisdom");
        
        binGals.binNGP(gals, N_grid, L, n_gals, nbw_gal);
        binGals.binNGP(rans, N_grid, L, n_rans, nbw_ran);
        
        double alpha = nbw_gal[0]/nbw_ran[0];
        double shotnoise = nbw_gal[1] + alpha*alpha*nbw_ran[1];
        nbw_gal[2] = (gals.size()*nbw_gal[1])/(L[0]*L[1]*L[2]);
        
        for (size_t i = 0; i < n_gals.size(); ++i) {
            n_gals[i] -= alpha*n_rans[i];
        }
        
        fftw_execute(dr2dk);
        
        double delta_k = 0.01;
        double k_min = 0.01;
        double k_max = 0.2;
        int N_bins = int((k_max - k_min)/delta_k + 1);
        std::vector<double> Pk(N_bins);
        std::vector<int> Nk(N_bins);
        
        std::vector<double> k_x = fftFreq(N_grid[0], L[0]);
        std::vector<double> k_y = fftFreq(N_grid[1], L[1]);
        std::vector<double> k_z = fftFreq(N_grid[2], L[2]);    
        
        
        for (int i = 0; i < N_grid[0]; ++i) {
            for (int j = 0; j < N_grid[1]; ++j) {
                for (int k = 0; k <= N_grid[2]/2; ++k) {
                    double k_mag = std::sqrt(k_x[i]*k_x[i] + k_y[j]*k_y[j] + k_z[k]*k_z[k]);
                    
                    if (k_mag >= k_min and k_mag < k_max) {
                        int index = k + (N_grid[2]/2 + 1)*(j + N_grid[1]*i);
                        int bin = (k_mag - k_min)/delta_k;
                        
                        Pk[bin] += (dk[index][0]*dk[index][0] + dk[index][1]*dk[index][1]) - shotnoise;
                        Nk[bin]++;
                    }
                }
            }
        }    
        
        std::cout << nbw_gal[2] << "\n";
        //std::ofstream
        
        std::string outfile = filename("Pk_LNKNLogs_", 4, sample + 1, "dat");     
        std::ofstream fout(outfile);
        fout.precision(15);
        for (int i = 0; i < N_bins; ++i) {
            if (Nk[i] > 0) {
                double k = k_min + (i + 0.5)*delta_k;
                Pk[i] /= (Nk[i]*nbw_gal[2]);
                fout << k << " " << Pk[i] << "\n";
            }
        }
    fout.close();
    }
    
}
