# nohs
The "Non-Orthogonal Hermite Solver" library (NOHS) is a simple C++ library that implements the solution of the eigenvalue problem associated with a non-periodic one-dimensional quantum Hamiltonian by expansion of the problem over a set of non-orthogonal Hermite functions. The library automatically applies the canonical orthogonalization procedure to ensure numerical stability.

## Prerequisites
* C++ compiler compatible with C++11 or higher (possibly g++ ver. 10.3 or newer)
* cmake 3.12 or newer
* [armadillo](http://arma.sourceforge.net/) (possibly ver. 10.6 or newer)
* [gsl](https://www.gnu.org/software/gsl/) (possibly ver. 2.6 or newer)
* [OpenMP](https://www.openmp.org/) (needed for parallel execution)

NOTE: all the version requirements listed with "(possibly ver. [...])" are not strictly mandatory. The version indicated represents the one used for development and testing. Older library versions can possibly work as well.

## How to install
After downloading this repository you can enter the main folder and use the following commands to compile and install the library:
```
mkdir build
cd build
cmake ..
make
sudo make install
```
This will compile and install the library in the `/usr/local` system folder (requires `sudo` privileges). You can change the installation location setting the `CMAKE_INSTALL_PREFIX` variable using the command:
```
cmake -DCMAKE_INSTALL_PREFIX=<your path> ..
```
If you want to compile the provided examples you can set the `COMPILE_EXAMPLES` variable to `ON` and set the installation path setting the correspondent `EXAMPLES_INSTALL_PATH` variable.

If you want, you can use also the `ccmake` interactive mode to set the required variables using the GUI.


## The `Hermite` function calss
The `Hermite` class, contained in the `nohs` namespace, implements the definition of the Hermite function:

<div align="center">
    <img src="https://render.githubusercontent.com/render/math?math=\psi_n(x):=\sqrt{\frac{\alpha}{2^nn!\sqrt{\pi}}}H_n(\alpha x)e^{-\frac{1}{2}\alpha^2x^2}">
</div>

where `H_n(q)` represents the [physicist's Hermite polynomial](https://en.wikipedia.org/wiki/Hermite_polynomials#Definition) of order `n`.

An instance of the `Hermite` class can be created using both the default unparameterized constructor or the parameterized constructor:

```
Hermite(int order_, double alpha_, double center_)
```

the latter sets both the order `n` and the amplitude coefficient `alpha`. The `center` parameter can be used to specify an arbitrary offset of the function center from the origin.

The `Hermite` class provides to the user the following public member functions:

| Function |      Description      |
|:----------:|:----------:|
| `double f(double x_)` |  returns the value of the function computed at the point `x_` |
| `double d1f(double x_)` |  returns the value of the function first derivative computed at the point `x_` |
| `double d2f(double x_)` |  returns the value of the function second-derivative computed at the point `x_` |

In order to ensure numerical stability the value of each Hermite function is computed recursively, starting from the function definitions for orders 0 and 1. The value of the first and second derivative terms is computed by the composition of Hermite function values.

WARNING: Acting upon a non-initialized Hermite object will result in a `nohs::exceptions::InitError` type exception.

## The `Solver` class
The `Solver` class, defined in the `nohs` namespace, implements the routines required to solve the Hamiltonian eigenvalue problem associated with a given potential function `V(x)`. The potential function  must be implemented according to the form:
```
double V(double x, void* pvoid)
```
in which the `pvoid` variable represents a pointer to some form of parameter required for the potential parametrization. 

An instance of the `Solver` class can be created specifying the pointer to the potential function considered, a pointer to the parameters required by the potential function definition and an indication of the basis set to be used in the computation. The basis set can be specified directly on construction, (as a `std::vector<nohs::Hermite>` vector), using the constructor:
```
Solver(std::vector<nohs::Hermite> BasisSet_, double (*V_)(double, void*), void* parameters_)        
```
or can be simply pre-allocated using the constructor:
```
Solver(unsigned int N_, double (*V_)(double, void*), void* parameters_)
```
If the second constructor (pre-allocation) is employed, the basis-set can be specified at a later time by adding the preselected number of basis functions through the function:
```
void add(nohs::Hermite function_)
```

The computation of the overlap and Hamiltonian matrix elements is performed by the [QAGI integration routine](https://www.gnu.org/software/gsl/doc/html/integration.html#qagi-adaptive-integration-on-infinite-intervals) of the GNU-GSL mathematical library. The integration parameters can be set using the function:
```
void set_integration_parameters(unsigned int npt_, double abs_, double rel_)
```
If the `set_integration_parameters` is not called the following set of integration parameters will be considered as defaults:
* Maximum number of integration points (`npt_`): `10000`
* Relative error threshold (`rel_`): `1e-10`
* Absolute error threshold (`abs_`): `1e-10`

Once the wanted solver configuration is obtained the solution of the problem can be computed calling the function:
```
void solve(double threshold_)
``` 
in which the `threshold_` parameter represents the overlap eigenvalue threshold adopted during the canonical orthogonalization process.

Once the problem solution has been completed the following set of public member functions can be interrogated:
| Function |      Description      |
|:----------:|:----------:|
| `int get_N_reduced()` |  returns the number of overlap-matrix eigenvalues not discarded during the canonical orthogonalization process |
| `double energy(int index_)` |  returns the `index_`-th energy eigenvalue |
| `double psi(int index_, double x_)` |  returns the value of the `index_`-th eigenfunction of the system computed at the point `x_` |

Invoking the previous functions without a previous call to solve will result in a `nohs::exceptions::SolverError` exception. If an invalid `index_` is specified in accessing the computed data a `nohs::exceptions::BoundError` exception will be raised.

## The `Optimizer` class
The `Optimizer` class, defined in the `nohs` namespace, is a simple auxiliary class that allows the user to pre-optimize the parameters of a small basis-set. An object of the `Optimizer` class can be initialized using the constructor:
```
Optimizer(double (*V_)(double, void*), void* parameters_)
```
in which, in perfect analogy with the `Solver` class, the target potential function is passed through the `V_` function pointer while the required parameters are encoded into the `parameters_` pointer.

The basis set to be optimized is specified in terms of groups of Hermite basis functions. Each group of functions is specified using the function:
```
void add(double center_, int max_order_, double guess_, int label_);
```
in which `center_` encodes the position of the function center, the `max_orede_` variable set the maximum order of the Hermite functions associated with the function group while `guess_` encodes a starting value for the group amplitude parameter `alpha`. The `label_` parameter can be used to link groups of functions that must share the same amplitude parameter (e.g. functions group related by symmetry). If the `label_` associated with each group is different from the one associated with the others, no link between parameters will be considered.

The optimization procedure generates small instance of the `Solver` class during its execution. The integration parameters considered during the optimization procedure can be set by calling the function:
```
void set_integration_parameters(unsigned int npt_, double abs_, double rel_)
```
whose definition and default parameters match the one considered in the `Solver` class definition.

Once the Hermite function groups have been defined, the optimization can be performed by invoking the:
```
void optimize(size_t max_iter_, double stop_size_, bool verbose_)
```
The optimization procedure is performed by the `gsl_multimin_fminimizer_nmsimplex2` [GNU_GSL routine](https://www.gnu.org/software/gsl/doc/html/multimin.html#algorithms-without-derivatives). The maximum number of iterations is set by the `max_iter_` variable while the convergence criterion (`gsl_multimin_test_size`) is set by the `stop_size_` variable. If the `verbose_` flag is set to `true` some information about the progress of the optimization process will be printed in the `iostream`.

Once the optimization procedure is completed the function:
```
std::vector<Hermite> generate_basis_set(std::vector<int> max_order_list_)
```
can be invoked in order to obtain a pre-optimized basis-set that can be used directly in the initialization of a `Solver` class object. The `max_order_list_` vector must be used to indicate the number fo functions required for each group.