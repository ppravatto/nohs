#include <iostream>
#include <fstream>
#include <cmath>
#include "nohs.h"

double Quartic(double x, void* pvoid){
    double delta = *(double *) pvoid;
    return delta*std::pow(x*x - 1., 2.);
}

int main(){

    double barrier = 40.;
    
    int N_min = 20;
    int N_max = 40;
    double alpha_min = 4.;
    double alpha_max = 4.;

    int N = 2*N_min + N_max;

    nohs::Solver System(N, &Quartic, &barrier);

    for(int i=0; i<N_max; i++){
        nohs::Hermite funct(i, alpha_max, 0.);
        System.add(funct);
    }

    for(int i=0; i<N_min; i++){
        nohs::Hermite funct_L(i, alpha_min, 1.);
        System.add(funct_L);
        nohs::Hermite funct_R(i, alpha_min, -1.);
        System.add(funct_R);
    }

    System.solve(1e-8);

    std::cout << "Barrier: " << barrier << std::endl;
    std::cout << "Basis functions: " << N << std::endl;
    std::cout << "Effective Basis functions: " << System.get_N_reduced() << std::endl << std::endl;

    std::cout << "n" << '\t' << "Energy" << std::endl;
    std::cout << "-------------------------------" << std::endl;
    for(int i=0; i<10; i++){
        std::cout << i << ")\t" << System.energy(i) << std::endl;
    }
    std::cout << std::endl;

    double TS = System.energy(1)-System.energy(0);

    std::cout << "Ground-state tunneling splitting: " << TS << std::endl << std::endl;;

    int npt_plot = 10000;
    int max_psi = 4;
    std::ofstream file("quartic.txt");
    for(int i=0; i<=npt_plot; i++){
        double x = -2.5 + 5.*i/npt_plot;
        file << x << '\t' << Quartic(x, &barrier);
        for(int j=0; j<=max_psi; j++) file << '\t' << System.psi(j, x);
        file << std::endl;
    }
    file.close();

    return 0.;
}