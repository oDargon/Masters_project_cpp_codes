#include <Eigen/core>
using namespace Eigen;


#if !defined(MYLIB_HELPER_STRUCTS_H)
#define MYLIB_HELPER_STRUCTS_H 1

struct ee_pair {
	double                     energy1;
	double                     energy2;
};

struct ed_pair {
	double                     energy;
	Matrix<double, Dynamic, 1> density;
};

struct dd_pair {
	Matrix<double, Dynamic, 1> density1;
	Matrix<double, Dynamic, 1> density2;
};

struct edd_pair {
	Matrix<double, Dynamic, 1> density1;
	Matrix<double, Dynamic, 1> density2;
	double                     energy_tot;
	
};

struct eedd_pair {
	Matrix<double, Dynamic, 1> density1;
	Matrix<double, Dynamic, 1> density2;
	std::vector<double>        energy1;
	std::vector<double>        energy2;
};

#endif