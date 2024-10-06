#include <numbers>
#include <fstream>
#include <Eigen/core>
#include <string>

#if !defined(MYLIB_CONSTANTS_H)
#define MYLIB_CONSTANTS_H 1

const int N_x      = 6;
const int L	       = 100;
//const double U     = 0;

const double t     = 1;
const double alpha = 3;
const double a     = 1;
const double pi    = std::numbers::pi;


const double R     = L * a / (2 * pi);
const std::string destination = "C:/Users/DzJas/Desktop/Code_Prjcts/Masters_project_py/Flexo_electric/cpp_results_ploting/data_storage/";


#endif