#ifndef NOHS_OPTIMIZATION_H
#define NOHS_OPTIMIZATION_H

#include <vector>
#include "nohs.h"
#include "nohs_exceptions.h"

namespace nohs{

    class Optimizer{
        private:
            bool optimized;
            int N_labels;
            std::vector<int> max_order, label;
            std::vector<double> center, guess, optimized_alpha;
            double (*V)(double, void*);
            void* parameters;
        public:
            Optimizer(double (*V_)(double, void*), void* parameters_);
            void add(double center_, int max_order_, double guess_, int label_);
            void optimize(bool verbose_);
            std::vector<Hermite> generate_basis_set(std::vector<int> max_order_list_);
    };

}

#endif