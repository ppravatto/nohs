#include <cmath>
#include "nohs.h"
#include "nohs_exceptions.h"

namespace nohs{

    static double hermite(double x, double alpha, int order){
        if(order == 0) return sqrt(alpha/sqrt(M_PI))*exp(-0.5*pow(alpha*x, 2.));						    //0-th order Hermite function
        else if(order == 1) return sqrt(alpha/(2.*sqrt(M_PI)))*2.*x*alpha*exp(-0.5*pow(alpha*x, 2.));		//1-st order hermite function
        else{
            double a = hermite(x, alpha, 0);												                //Set the (n-1)-th term to the value of the 0-th Hermite function
            double b = hermite(x, alpha, 1);												                //Set the n-th term to the value of the 1-st Hermite function
            for(int i=2; i<=order; i++){                                                                    //Iterate until the order-th function of the recursion series
                double var = alpha*x*sqrt(2./i)*b - sqrt((i-1.)/i)*a;			                            //Apply the recursion relation to compute the (n+1)-th term
                a = b;															                            //Set the current n-th term to the (n-1)-th	position
                b = var;																                    //Copy the obtained result as the n-th value (needed for the next iteration)
            }
            return b;																	                    //Return the function value
        }
    }

    Hermite::Hermite() : init(false), order(0), alpha(0.) {}

    Hermite::Hermite(int order_, double alpha_, double center_) : init(true), order(order_), alpha(alpha_), center(center_) {}

    double Hermite::f(double x_) {
        if(init == false) throw exceptions::InitError();
        return hermite(x_-center, alpha, order);
    }

    double Hermite::d1f(double x_) {
        if(init == false) throw exceptions::InitError();
        double y = x_-center;
        double b = -std::sqrt((order+1.)/2.)*hermite(y, alpha, order+1);
        if(order==0) return alpha*b;
        else return alpha*(std::sqrt(order/2.)*hermite(y, alpha, order-1) + b);
    }

    double Hermite::d2f(double x_) {
        if(init == false) throw exceptions::InitError();
        double y = x_-center;
        double c = std::sqrt((order+2.)*(order+1.))*hermite(y, alpha, order+2);
        double b = (2.*order+1.)*hermite(y, alpha, order);
        if(order==0) return -0.5*std::pow(alpha, 2.)*((order+1.)*hermite(y, alpha, order) - c);
        else if(order==1) return 0.5*std::pow(alpha, 2.)*(c-b);
        else{
            double a = std::sqrt(order*(order-1))*hermite(y, alpha, order-2);
            return 0.5*std::pow(alpha, 2.)*(a-b+c);
        }
    }

}