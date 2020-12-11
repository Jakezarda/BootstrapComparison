#ifndef _BINNER_HPP_
#define _BINNER_HPP_

#include <vector>
#include "galaxy.hpp"

class binner{
public:
    binner();
    
    void binNGP(std::vector<galaxy> &gals, std::vector<int> N_grid, std::vector<double> L, std::vector<double> &n, std::vector<double> &nbw);

};

binner::binner() {
    
}

void binner::binNGP(std::vector<galaxy> &gals, std::vector<int> N_grid, std::vector<double> L, std::vector<double> &n, std::vector<double> &nbw) {
    std::vector<double> dr = {L[0]/N_grid[0], L[1]/N_grid[1], L[2]/N_grid[2]};
    
    for (size_t i = 0; i < gals.size(); ++i) {
        std::vector<double> pos = gals[i].getCartPos();
        double w = gals[i].getWeight();
        double nbar = gals[i].getNbar();
        int x = pos[0]/dr[0];
        int y = pos[1]/dr[1];
        int z = pos[2]/dr[2];
        
        if (x >= N_grid[0]) {
            x -= 1;
        }
        if (y >= N_grid[1]) {
            y -= 1;
        }
        if (z >= N_grid[2]) {
            z -= 1;
        }
        
        int index = z + N_grid[2]*(y + N_grid[1]*x);
        n[index] += w;
        
        nbw[0] += w;
        nbw[1] += w*w;
        nbw[2] += nbar*w*w;
    }
}

#endif
