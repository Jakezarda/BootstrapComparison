#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <random>
#include <fftw3.h>
#include <sstream>

#include "binner.hpp"
#include "galaxy.hpp"

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

std::string filename(std::string base, int digits, int num, std::string ext){
    std::stringstream file;
    file << base << std::setw(digits) << std::setfill('0') << num << "." << ext;       //creates files with names such that filename0001.extension
    return file.str();
}

int main() {
    
    std::vector<int> N_grid = {512, 512, 512};
    std::vector<double> L = {1024.0, 1024.0, 1024.0};
    double nbar = 0.0004;
    
    std::vector<galaxy> gals = fileRead("LNKNLogs_0001.dat");
    std::vector<galaxy> rans = fileRead("LNKNLogsVel_Random.dat");
    
    binner binGals;
    std::vector<double> n_gals(N_grid[0]*N_grid[1]*N_grid[2]);
    std::vector<double> n_rans(N_grid[0]*N_grid[1]*N_grid[2]);
    std::vector<double> nbw_gal(3);
    std::vector<double> nbw_ran(3);
    
    std::vector<fftw_complex> dk(N_grid[0]*N_grid[1]*(N_grid[2]/2 + 1));
    
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
    
    std::vector<std::vector<double>> Pk_bins(N_bins);               //vector of vectors to store individual values
    
    for (int i = 0; i < N_grid[0]; ++i) {
        for (int j = 0; j < N_grid[1]; ++j) {
            for (int k = 0; k <= N_grid[2]/2; ++k) {
                double k_mag = std::sqrt(k_x[i]*k_x[i] + k_y[j]*k_y[j] + k_z[k]*k_z[k]);
                
                if (k_mag >= k_min and k_mag < k_max) {
                    int index = k + (N_grid[2]/2 + 1)*(j + N_grid[1]*i);
                    int bin = (k_mag - k_min)/delta_k;
                    
                    Pk[bin] += (dk[index][0]*dk[index][0] + dk[index][1]*dk[index][1]) - shotnoise;
                    Nk[bin]++;
                    
                    Pk_bins[bin].push_back(dk[index][0]*dk[index][0] + dk[index][1]*dk[index][1] - shotnoise);
                }
            }
        }
    }
    
    
    std::random_device seeder;
    std::mt19937_64 gen(seeder());
    
    for (int bootstrap = 0; bootstrap < 1000; ++bootstrap) {
        std::vector<double> Pk_boot(N_bins);
        
        for (int i = 0; i < N_bins; ++i) {
            std::uniform_int_distribution<int> dist(0, Pk_bins[i].size());
            
            for (int j = 0; j < Pk_bins[i].size(); ++j) {
                Pk_boot[i] += Pk_bins[i][dist(gen)]/double(Pk_bins[i].size());
            }
            
        }
        
        std::string outfile = filename("Pk_bootstrap", 4, bootstrap + 1, "dat");     
        std::ofstream fout(outfile);
        fout.precision(15);
        
        for (int l = 0; l < N_bins; ++l) {
            if (Nk[l] > 0) {
                double k = k_min + (l + 0.5)*delta_k;
                //Pk_boot[l] /= (Nk[l]*nbw_gal[2]);
                //std::cout << Nk[l] << std::endl;
                fout << k << " " << Pk_boot[l] << "\n";
            }
            
        }

    }
    
    
    std::cout << nbw_gal[2] << "\n";
    std::ofstream fout("Bootstrap_single_FFT.dat");                   
    fout.precision(15);
    std::cout << Pk[1] << " " << Nk[1] << " " << nbw_gal[2] << std::endl;
    for (int i = 0; i < N_bins; ++i) {
        //std::cout << Nk[i] << "\n";
        if (Nk[i] > 0) {
            std::cout << Pk[i] / (Nk[i]*nbw_gal[2]) << std::endl;
            double k = k_min + (i + 0.5)*delta_k;
            Pk[i] /= (Nk[i]*nbw_gal[2]);
            fout << k << " " << Pk[i] << "\n";
        }
    }
    fout.close();    
}
