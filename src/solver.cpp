#include <cmath>
#include <armadillo>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>

#ifdef _OPENMP
    #include <omp.h>
#endif

#include "nohs.h"
#include "nohs_exceptions.h"

namespace nohs{

    static double QAGI_integrator(double (*f)(double, void*), void * pvoid, int npt, double abs, double rel){
        double result, error;
        gsl_function integrand;
        integrand.function = f;
        integrand.params = pvoid;
        gsl_integration_workspace* w = gsl_integration_workspace_alloc(npt);
        gsl_integration_qagi(&integrand, abs, rel, npt, w, &result, &error);
        gsl_integration_workspace_free(w);
        return result;
    }

    double Solver::overlap_integrand(int row, int col, double x){
        return BasisSet[row].f(x)*BasisSet[col].f(x);
    }

    double Solver::hamiltonian_integrand(int row, int col, double x){
        return BasisSet[row].f(x)*(-BasisSet[col].d2f(x) + V(x, parameters)*BasisSet[col].f(x));
    }

    Solver::Solver(unsigned int N_, double (*V_)(double, void*), void* parameters_) : solved(false), N(N_), N_red(-1), npt(10000), abs(1e-10), rel(1e-10), V(V_), parameters(parameters_) {
        H = arma::mat(N, N, arma::fill::zeros);
        S = arma::mat(N, N, arma::fill::zeros);
        BasisSet.reserve(N);
    }

    Solver::Solver(std::vector<Hermite> BasisSet_, double (*V_)(double, void*), void* parameters_) : Solver(BasisSet_.size(), V_, parameters_) {
        BasisSet = BasisSet_;
    }

    void Solver::set_integration_parameters(unsigned int npt_, double abs_, double rel_){
        npt = npt_; abs = abs_, rel = rel_;
    }

    void Solver::add(Hermite function_){
        if(int(BasisSet.size()) >= N) throw exceptions::MaxDimensionError();
        BasisSet.push_back(function_);
    }

    void Solver::solve(double threshold_){

        #ifdef _OPENMP
            #pragma omp parallel for collapse(2)
        #endif
        for(int row=0; row<N; row++){
            for(int col=0; col<N; col++){
                utils::Carrier data(this, &Solver::overlap_integrand, row , col);
                S(row, col) = QAGI_integrator(&utils::aux_func_gsl, &data, npt, abs, rel);
                data.fpt = &Solver::hamiltonian_integrand;
                H(row, col) = QAGI_integrator(&utils::aux_func_gsl, &data, npt, abs, rel);
            }
        }

        #ifdef _OPENMP
            #pragma omp parallel for
        #endif
        for(int row=0; row<N; row++){
            for(int col=0; col<row; col++){
                double S_average = 0.5*(S(row, col) + S(col, row));
                S(row, col) = S_average;
                S(col, row) = S_average;
                double H_average = 0.5*(H(row, col) + H(col, row));
                H(row, col) = H_average;
                H(col, row) = H_average;
            }
        }

        arma::mat Q = arma::mat(N, N, arma::fill::zeros);
        arma::vec S_eval = arma::vec (N, arma::fill::zeros);
        arma::eig_sym(S_eval, Q, S, "std");

        int stop_index = N-1;
        while(stop_index > 0){
            if(S_eval(stop_index-1) < threshold_){
                break;
            }
            stop_index--;
        }

        N_red = N - stop_index;
        arma::mat Q_red = arma::mat(N, N_red, arma::fill::zeros);

        #ifdef _OPENMP
            #pragma omp parallel for
        #endif
        for(int col=0; col<N_red; col++){
            double factor = 1./std::sqrt(S_eval(stop_index + col));
            for(int row=0; row<N; row++){
                Q_red(row, col) = Q(row, stop_index+col)*factor;
            }
        }

        arma::mat H_red = arma::mat(N_red, N_red, arma::fill::zeros);
        H_red = Q_red.t() * H * Q_red;

        #ifdef _OPENMP
            #pragma omp parallel for
        #endif
        for(int row=0; row<N_red; row++){
            for(int col=0; col<=row; col++){
                double average = 0.5*(H_red(row, col) + H_red(col, row));
                H_red(row, col) = average;
                H_red(col, row) = average;
            }
        }

        arma::mat C_red = arma::mat(N_red, N_red, arma::fill::zeros);
        E = arma::vec(N_red, arma::fill::zeros);
        arma::eig_sym(E, C_red, H_red, "std");
        C = arma::mat(N, N_red, arma::fill::zeros);
        C = Q_red * C_red;

        solved = true;
    }

    int Solver::get_N_reduced(){
        if(solved == false) throw exceptions::SolverError();
        return N_red;
    }

    double Solver::energy(int index_){
        if(solved == false) throw exceptions::SolverError();
        if(index_ >= N_red) throw exceptions::BoundError();
        return E(index_);
    }

    double Solver::psi(int index_, double x_){
        if(solved == false) throw exceptions::SolverError();
        if(index_ >= N_red) throw exceptions::BoundError();
        double value = 0.;
        for(int i=0; i<N; i++){
            value += C(i, index_)*BasisSet[i].f(x_);
        }
        return value;
    }

}