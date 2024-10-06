#include <Eigen/core>
using namespace Eigen;


#if !defined(MYLIB_HELPER_STRUCTS_H)
#define MYLIB_HELPER_STRUCTS_H 1

struct ed_pair {
	double                     energy;
	Matrix<double, Dynamic, 1> density;
};

struct dd_pair {
	Matrix<double, Dynamic, 1> density1;
	Matrix<double, Dynamic, 1> density2;
};

struct eedd_pair {
	Matrix<double, Dynamic, 1> density1;
	Matrix<double, Dynamic, 1> density2;
	std::vector<double>        energy1;
	std::vector<double>        energy2;
};

#endif