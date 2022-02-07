#include <iostream>
#include <fstream>
#include <cmath>
#include "nohs.h"
#include "nohs_optimizer.h"

static double Quartic(double x, void* pvoid){
    double delta = *(double *) pvoid;
    return delta*std::pow(x*x - 1., 2.);
}

static double TS_Garg(double barrier){
    return (32./std::sqrt(2.*M_PI))*std::pow(barrier, 3./4.)*std::exp(-(4./3.)*std::sqrt(barrier));
}

int main(){

    double barrier = 40.;
    
    int N_min = 20;
    int N_max = 40;

    double alpha_guess = 3.;

    int N = 2*N_min + N_max;

    std::cout << "Barrier: " << barrier << std::endl;
    std::cout << "Basis functions: " << N << std::endl << std::endl;;
    std::cout << "Basis-set optimization:" << std::endl;
    std::cout << "---------------------------------------" << std::endl;

    nohs::Optimizer MyOptimizer(&Quartic, &barrier);
    MyOptimizer.add(0., 10, alpha_guess, 0);
    MyOptimizer.add(1., 5, alpha_guess, 1);
    MyOptimizer.add(-1., 5, alpha_guess, 1);
    MyOptimizer.optimize(true);
    std::cout << std::endl;


    std::vector<int> order = {N_max, N_min, N_min};
    std::vector<nohs::Hermite> BasisSet = MyOptimizer.generate_basis_set(order);
    nohs::Solver System(BasisSet, &Quartic, &barrier);
    System.solve(1e-8);

    std::cout << "Effective Basis functions: " << System.get_N_reduced() << std::endl << std::endl;

    std::cout << "n" << '\t' << "Energy" << std::endl;
    std::cout << "---------------------------------------" << std::endl;
    for(int i=0; i<10; i++){
        std::cout << i << ")\t" << System.energy(i) << std::endl;
    }
    std::cout << std::endl;

    double TS = System.energy(1)-System.energy(0);

    std::cout << "Ground-state tunneling splitting:" << std::endl;
    std::cout << "Hermite: " << TS << std::endl;
    std::cout << "WKB (Garg): " << TS_Garg(barrier) << std::endl << std::endl;

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