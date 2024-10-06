#include "model.h"
#include "constants.h"
#include <array>
#include <iostream>


void Model::setup() {
	
	tube_center << R, 0, 0;
	cell_basis.resize(N_x, 3);
	brillouin_zone.resize(1, L);
	distance_change.resize(1, N_x);
	cos_arr.resize(1, L);


	std::array<Matrix<double, Dynamic, 3>,2> positions;


	for (int i{ 0 }; i < N_x; ++i) {
		RowVector3d temp_vec;
		temp_vec << std::cos(2 * pi * i / N_x), 0, std::sin(2 * pi * i / N_x);

		cell_basis.row(i) = a * temp_vec;
	}

	brillouin_zone.setLinSpaced(L, -pi, pi - (2 * pi / L));

	for (int i{ 0 }; i < 2; ++i) {

		positions[i].resize(N_x, 3);

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

	tk_mat.resize(N_x, N_x);
	indices.resize(1, L * N_x);
	eigen_vals.resize(1, L * N_x);
	eigen_vecs.resize(N_x, L * N_x);
	

	
}

