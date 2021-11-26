#include <fstream>
#include <sstream>
#include "nohs.h"

int main(int argc, char* argv[]){

    const int npt = 10000;
    const double exc = 10.;

    if(argc != 2){
        std::cout << "ERROR: the program expects a single argument (Hermite function order)" << std::endl;
        return -1;
    }

    int order;
    std::stringstream(argv[1]) >> order;

    nohs::Hermite function(order, 1., 0.);

    std::ofstream file("Herm.txt");

    for(int i=0; i<=npt; i++){
        double x = -exc + 2.*exc*i/npt;
        file << x << '\t' 
            << function.f(x) << '\t' 
            << function.d1f(x) << '\t' 
            << function.d2f(x) << std::endl;
    }

    file.close();

    return 0;
}