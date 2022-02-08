#ifndef NOHS_OPTIMIZATION_H
#define NOHS_OPTIMIZATION_H

#include <vector>
#include <gsl/gsl_multimin.h>

#include "nohs.h"
#include "nohs_exceptions.h"

namespace nohs{

    double minimization_target(const gsl_vector* v, void* pvoid);

    class Optimizer{
        private:
            bool optimized;
            int N_labels, npt;
            double abs, rel;
            std::vector<int> max_order, label;
            std::vector<double> center, guess, optimized_alpha;
            double (*V)(double, void*);
            void* parameters;
        public:
            Optimizer(double (*V_)(double, void*), void* parameters_);
            void add(double center_, int max_order_, double guess_, int label_);
            void set_integration_parameters(unsigned int npt_, double abs_, double rel_);
            void optimize(size_t max_iter_, double stop_size_, bool verbose_);
            std::vector<Hermite> generate_basis_set(std::vector<int> max_order_list_);

            friend double minimization_target(const gsl_vector* v, void* pvoid);
    };

}

#endif