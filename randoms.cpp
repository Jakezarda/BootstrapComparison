#include <iostream>
#include <fstream>
#include <random>

int main() {
    double L = 1024.0;
    double nbar = 0.0004;
    double timesRan = 10.0;
    
    size_t N_ran = nbar*L*L*L*timesRan;
    
    std::random_device seeder;
    std::mt19937_64 gen(seeder());
    std::uniform_real_distribution<double> dist (0.0, 1.0);
    
    std::ofstream fout("LNKNLogsVel_Random.dat");
    fout.precision(15);
    for (size_t i = 0; i < N_ran; ++i) {
        double x = L*dist(gen);
        double y = L*dist(gen);
        double z = L*dist(gen);
        
        fout << x << " " << y << " " << " " << z << " " << 0.0 << " " << 0.0 << " " << 0.0 << "\n";
        
    }
    fout.close();
    
    return 0;
}
