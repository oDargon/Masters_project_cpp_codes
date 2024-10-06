#include <Eigen/core>
#include "Constants.h"
#include <vector>
#include <array>


#if !defined(MYLIB_SETUP_H)
#define MYLIB_SETUP_H 1

using namespace Eigen;

extern Matrix<double, 3, 1>                  tube_center;
extern Matrix<double, N_x, 3>                cell_basis;
extern std::array<Matrix<double, N_x, 3>, L> positions;
extern Matrix<double, 1, L>				     brillouin_zone;
extern Matrix<double, 1, N_x>				 distance_change;
extern Matrix<double, 1, L>				     cos_arr;

void set_up();

#endif