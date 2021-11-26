#ifndef NOHS_EXCEPTIONS_H
#define NOHS_EXCEPTIONS_H

#include <exception>

namespace nohs{
    namespace exceptions{

        class InitError: public std::exception {
            virtual const char* what() const throw(){
                return "Cannot operate on an uninitilized object";
            }
        };

        class MaxDimensionError: public std::exception {
            virtual const char* what() const throw(){
                return "Index exceeds the maximum object size";
            }
        };

        class SolverError: public std::exception {
            virtual const char* what() const throw(){
                return "Cannot access requested data without a prior call to solve()";
            }
        };

        class BoundError: public std::exception {
            virtual const char* what() const throw(){
                return "Index out of bounds";
            }
        };      

    }
}

#endif