#ifndef _GALAXY_HPP_
#define _GALAXY_HPP_

#include <iostream>
#include <vector>
#include <string>

class galaxy{
    std::vector<double> r_cart, r_astro, vel;
    double w, nbar;
    
public:
    galaxy(std::vector<double> &pos, double w, double nbar, std::string posType);
    
    galaxy(std::vector<double> &pos, std::vector<double> &vel, double w, double nbar, std::string posType);
    
    std::vector<double> getCartPos();
    
    std::vector<double> getAstroPos();
    
    std::vector<double> getVel();
    
    double getWeight();
    
    double getNbar();
    
    void setCartPos(std::vector<double> &pos);
    
    void setAstroPos(std::vector<double> &pos);
    
    void setVel(std::vector<double> &vel);
    
    void setWeight(double w);
};

galaxy::galaxy(std::vector<double> &pos, double w, double nbar, std::string posType) {
    this->w = w;
    this->nbar = nbar;
    
    if (posType == "cartesian" or posType == "Cartesian") {
        this->r_cart = pos;
    } else if (posType == "Astro" or posType == "astro") {
        this->r_astro = pos;
    } else {
        std::cout << "Invalid position type.\n";
        throw std::runtime_error("Invalid position type"); // Have this throw and error.
    }
}

galaxy::galaxy(std::vector<double> &pos, std::vector<double> &vel, double w, double nbar, std::string posType) {
    this->vel = std::vector<double>(3);
    this->r_cart = std::vector<double>(3);
    this->r_astro = std::vector<double>(3);
    this->w = w;
    this->vel = vel;
    this->nbar = nbar;
    
    if (posType == "cartesian" or posType == "Cartesian") {
        this->r_cart = pos;
    } else if (posType == "Astro" or posType == "astro") {
        this->r_astro = pos;
    } else {
        std::cout << "Invalid position type.\n";
        throw std::runtime_error("Invalid position type");// Have this throw and error.
    }
}

std::vector<double> galaxy::getCartPos() {
    return this->r_cart;
}

std::vector<double> galaxy::getAstroPos() {
    return this->r_astro;
}

std::vector<double> galaxy::getVel() {
    return this->vel;
}

double galaxy::getWeight() {
    return this->w;
}

double galaxy::getNbar() {
    return this->nbar;
}

void galaxy::setCartPos(std::vector<double> &pos) {
    this->r_cart = pos;
}

void galaxy::setAstroPos(std::vector<double> &pos) {
    this->r_astro = pos;
}

void galaxy::setVel(std::vector<double> &vel) {
    this->vel = vel;
}

void galaxy::setWeight(double w) {
    this->w = w;
}

#endif
