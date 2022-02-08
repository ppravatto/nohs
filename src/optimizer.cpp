#include <cmath>
#include <vector>
#include <gsl/gsl_multimin.h>

#include "nohs_optimizer.h"
#include "nohs.h"

template <typename T>
static int find_index(std::vector<T> vector, T object){
    int index = -1.;
    for(int i=0; i<int(vector.size()); i++){
        if(vector[i] == object){
            index = i;
            break;
        }
    }
    return index;
}

namespace nohs{

    struct OptCarrier{
        Optimizer* opt_ptr;
        std::vector<int> assignment;
        OptCarrier(Optimizer* opt_ptr_, std::vector<int> assignment_) : opt_ptr(opt_ptr_), assignment(assignment_) {}
    };

    double minimization_target(const gsl_vector* v, void* pvoid){
        OptCarrier p = *(OptCarrier*) pvoid;
        Optimizer* origin = p.opt_ptr;

        int N = int(origin->max_order.size());

        std::vector<Hermite> BasisSet;
        for(int i=0; i<N; i++){
            double alpha = gsl_vector_get(v, p.assignment[i]);
            for(int order=0; order<=origin->max_order[i]; order++){
                Hermite function(order, alpha, origin->center[i]);
                BasisSet.push_back(function);
            }
        }

        Solver System(BasisSet, origin->V, origin->parameters);
        System.set_integration_parameters(origin->npt, origin->abs, origin->rel);
        System.solve(1e-8);

        return System.energy(0);
    }


    
    Optimizer::Optimizer(double (*V_)(double, void*), void* parameters_) : optimized(false), N_labels(0), npt(10000), abs(1e-10), rel(1e-10), V(V_), parameters(parameters_) {}

    void Optimizer::add(double center_, int max_order_, double guess_, int label_){
        int idx = find_index<int>(label, label_);
        if(idx == -1) N_labels++;
        else if(guess[idx] != guess_) throw exceptions::InvalidError();
        center.push_back(center_);
        max_order.push_back(max_order_);
        guess.push_back(guess_);
        label.push_back(label_);
    }

    void Optimizer::set_integration_parameters(unsigned int npt_, double abs_, double rel_){
        npt = npt_; abs = abs_, rel = rel_;
    }

    void Optimizer::optimize(size_t max_iter_, double stop_size_, bool verbose_){

        const gsl_multimin_fminimizer_type* T = gsl_multimin_fminimizer_nmsimplex2;
        gsl_multimin_fminimizer* s = gsl_multimin_fminimizer_alloc (T, N_labels);
        
        gsl_vector* step = gsl_vector_alloc(N_labels);
        gsl_vector_set_all (step, 0.5);
        
        gsl_vector* x = gsl_vector_alloc(N_labels);

        int idx = 0;
        std::vector<int> skip, assignment(label.size(), -1);
        for(int i=0; i<int(label.size()); i++){
            if(find_index(skip, label[i]) == -1){
                gsl_vector_set(x, idx, guess[i]);
                skip.push_back(label[i]);
                for(int j=0; j<int(label.size()); j++){
                    if(label[j] == label[i]) assignment[j] = idx;
                }
                idx++;
            } 
        }

        gsl_multimin_function target;
        target.n = N_labels;
        target.f = &minimization_target;
        OptCarrier data(this, assignment);
        target.params = &data;
        
        gsl_multimin_fminimizer_set(s, &target, x, step);

        size_t iter = 0;
        int status;
        double size;
        
        do{
            iter++;
            status = gsl_multimin_fminimizer_iterate(s);
            if (status) break;
            size = gsl_multimin_fminimizer_size(s);
            status = gsl_multimin_test_size (size, stop_size_);
            if(status == GSL_SUCCESS && verbose_==true){
                std::cout << "Optimization converged:" << std::endl;
                std::cout << "-> Iterations: " << iter << std::endl;
                std::cout << "-> Alpha values:" << std::endl;
                for(int i=0; i<N_labels; i++) std::cout << "        " << i << " -> " << gsl_vector_get(s->x, i) << std::endl;
                std::cout << "-> Function value: " << s->fval << std::endl;
                std::cout << "-> Simplex size: " << size << std::endl;
            }
        }while (status == GSL_CONTINUE && iter < max_iter_);

        optimized_alpha.clear();
        for(std::vector<int>::iterator ptr = assignment.begin(); ptr < assignment.end(); ptr++){
            optimized_alpha.push_back(gsl_vector_get(s->x, *ptr));
        }
        
        gsl_vector_free(x);
        gsl_vector_free(step);
        gsl_multimin_fminimizer_free(s);

        optimized = true;

    }

    std::vector<Hermite> Optimizer::generate_basis_set(std::vector<int> max_order_list_){
        if(optimized == false) throw exceptions::OptimizeError();
        if(max_order_list_.size() != max_order.size()) throw exceptions::InvalidError();
        std::vector<Hermite> BasisSet;
        for(int i=0; i<int(max_order_list_.size()); i++){
            for(int order=0; order<=max_order_list_[i]; order++){
                Hermite function(order, optimized_alpha[i], center[i]);
                BasisSet.push_back(function);
            }
        }
        return BasisSet;
    }

}