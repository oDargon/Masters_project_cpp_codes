#include <iostream>
#include <Eigen/core>
#include "Constants.h"
#include <vector>
#include <array>
#include <cmath>

using namespace Eigen;

Matrix<double, 3, 1>                  tube_center;
Matrix<double, N_x, 3>                cell_basis;
std::array<Matrix<double, N_x, 3>, L> positions;
Matrix<double, 1, L>				  brillouin_zone;
Matrix<double, 1, N_x>				  distance_change;
Matrix<double, 1, L>				  cos_arr;

void set_up() {
	tube_center << R, 0, 0;

	for (int i{ 0 }; i < N_x; ++i) {
		RowVector3d temp_vec;
		temp_vec << std::cos(2 * pi * i / N_x), 0, std::sin(2 * pi * i / N_x);

		cell_basis.row(i) = a * temp_vec;
	}

	brillouin_zone.setLinSpaced(L, -pi, pi - (2 * pi / L));

	for (int i{ 0 }; i < L; ++i) {

		double rot_cos = std::cos(2 * pi * i / L);
		double rot_sin = std::sin(2 * pi * i / L);

		Matrix3d z_rot_mat;
		z_rot_mat << rot_cos, rot_sin, 0, -rot_sin, rot_cos, 0, 0, 0, 1;

		for (int j{ 0 }; j < N_x; ++j) {
			positions[i].row(j) = z_rot_mat * (tube_center + cell_basis.row(j).transpose());

		}
	}

	for (int i{ 0 }; i < N_x; ++i) {
		distance_change(i) = (positions[0].row(i) - positions[1].row(i)).norm() - a;
	}

	for (int i{ 0 }; i < L; ++i) {
		cos_arr(i) = std::cos(brillouin_zone(i));
	}

	std::cout << distance_change << "\n";
}