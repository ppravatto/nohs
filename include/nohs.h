#ifndef NOHS_H
#define NOHS_H

#include <vector>
#include <armadillo>

namespace nohs{

    class Hermite;
    class Solver;

    namespace utils{

        struct Carrier{
            int r, c;
            Solver* cpt;
            double (Solver::* fpt)(int, int, double);
            Carrier() : r(0), c(0), cpt(nullptr), fpt(nullptr) {}
            Carrier(Solver* cpt_, double (Solver::* fpt_)(int, int, double), int r_, int c_) : r(r_), c(c_), cpt(cpt_), fpt(fpt_) {}
        };

        inline double aux_func_gsl(double x, void* pvoid){
            Carrier data = *(Carrier*) pvoid;
            return ((*(data.cpt)).*(data.fpt))(data.r, data.c, x);
        }

    }

    class Hermite{
        private:
            bool init;
            int order;
            double alpha, center;
        public:
            Hermite();
            Hermite(int order_, double alpha_, double center_);
            double f(double x_);
            double d1f(double x_);
            double d2f(double x_);
    };

    class Solver{
        private:
            bool solved;
            int N, N_red, npt;
            double abs, rel;
            std::vector<Hermite> BasisSet;
            arma::mat H, S, C;
            arma::vec E;
            double (*V)(double, void*);
            void* parameters;
            
            double overlap_integrand(int row, int col, double x);
            double hamiltonian_integrand(int row, int col, double x);

        public:
            Solver(unsigned int N_, double (*V_)(double, void*), void* parameters_);
            Solver(std::vector<Hermite> BasisSet_, double (*V_)(double, void*), void* parameters_);           
            void set_integration_parameters(unsigned int npt_, double abs_, double rel_);
            void add(Hermite function_);
            void solve(double threshold_);
            int get_N_reduced();
            double energy(int index_);
            double psi(int index_, double x_);

            friend double utils::aux_func_gsl(double x, void* pvoid);

    };

}

#endif