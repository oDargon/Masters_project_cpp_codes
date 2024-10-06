#include <Eigen/core>
#include <Eigen/Eigenvalues>
#include "model.h"
#include "constants.h"
#include <iostream>
#include <string>
#include <algorithm>
#include <fstream>
#include <string>
#include <format>
#include <cmath>
#include <iomanip>   

const static IOFormat CSVFormat(StreamPrecision, DontAlignCols, ", ", "\n");

Model::Model(int ring_N = 6, int length_N = 100, double hopping_t = 1, double alpha = 3, double spacing_a = 1)
	: N_x{ ring_N }, L{ length_N }, t{ hopping_t }, alpha{ alpha }, a{ spacing_a }, R{ L * a / (2 * pi)  } {

	setup();
}

std::vector<int> Model::gen_filling_factors(int end, int step, bool reverse) {

	std::vector<int> return_arr;
	int start = 2;

	while (start <= 2 * end) {
		return_arr.push_back(start);
		start += step * 2;
	}

	if (reverse) {
		std::reverse(return_arr.begin(), return_arr.end());
	}

	return return_arr;
}

void Model::tk_hubbard(Matrix<double, Dynamic, 1> rho_vec, int k, double U) {

	if (tkh_count == 0) {
		tk_mat.setZero();

		for (int i{ 0 }; i < N_x; ++i) {
			int neighbour = 0;
			if (i + 1 < N_x) { neighbour = i + 1; }

			tk_mat(i, neighbour) = -t;
			tk_mat(neighbour, i) = -t;
		}
	}

	for (int j{ 0 }; j < N_x; ++j) {
		tk_mat(j, j) = -t * (1 - alpha * distance_change(j)) * cos_arr(k);
	}

	tk_mat.diagonal() += U * rho_vec;

	tkh_count++;

}

ed_pair Model::extract_vecs(Matrix<double, Dynamic, 1> rho_vec, int fill_m, double U) {

	for (int i{ 0 }; i < L; ++i) {
		tk_hubbard(rho_vec, i, U);
		solver.compute(tk_mat);

		eigen_vecs.block(0, N_x * i, N_x, N_x) = solver.eigenvectors();
		eigen_vals.segment(i * N_x, N_x) = solver.eigenvalues().transpose();

	}

	eigen_vecs /= sqrt(L);
	indices.setLinSpaced(L * N_x, 0, L * N_x);
	std::sort(indices.data(), indices.data() + L * N_x, [&](int a, int b) { return eigen_vals(b) > eigen_vals(a); });

	for (int i = 0; i < L * N_x; ++i) {

		if (indices[i] == -1) continue;
		double value = eigen_vals(i);
		Matrix<double, Dynamic, 1 > value2 = eigen_vecs.col(i);
		int x = i;
		int y = indices(i);

		while (y != i) {
			indices(x) = -1;
			eigen_vals(x) = eigen_vals(y);
			eigen_vecs.col(x) = eigen_vecs.col(y);
			x = y;
			y = indices(x);
		}
		eigen_vals(x) = value;
		eigen_vecs.col(x) = value2;
		indices(x) = -1;

	}

	eigen_vecs = eigen_vecs.array().square();

	ed_pair pair;
	/*pair.energy  = eigen_vals.head(fill_m).sum() / fill_m;*/
	pair.energy = eigen_vals(0);
	pair.density = eigen_vecs.leftCols(fill_m).rowwise().sum();

	return pair;

}

dd_pair Model::hubbard_iteration(int electron_count, Matrix<double, Dynamic, 1> initial_density, double U) {

	int he_count = static_cast<int>(std::round(electron_count / 2));

	std::vector<double> up_energy;
	std::vector<double> down_energy;

	Matrix<double, Dynamic, 1> up_density;
	Matrix<double, Dynamic, 1> down_density;

	up_density = initial_density;
	ed_pair temp;

	for (int i{ 0 }; i < 100; ++i) {

		temp = extract_vecs(up_density, he_count, U);
		down_energy.push_back(temp.energy);
		down_density = temp.density;

		temp = extract_vecs(down_density, he_count, U);
		up_energy.push_back(temp.energy);
		up_density = temp.density;

	}

	dd_pair return_pair;
	return_pair.density1 = up_density;
	return_pair.density2 = down_density;

	std::cout << "done with: " << counter << "\n";
	counter++;
	return return_pair;
    
}

eedd_pair Model::hubbard_iteration_with_energy(int electron_count, Matrix<double, Dynamic, 1> initial_density, double U) {

	int he_count = static_cast<int>(std::round(electron_count / 2));

	std::vector<double> up_energy;
	std::vector<double> down_energy;

	Matrix<double, Dynamic, 1> up_density;
	Matrix<double, Dynamic, 1> down_density;

	up_density = initial_density;
	ed_pair temp;

	for (int i{ 0 }; i < 500; ++i) {

		temp = extract_vecs(up_density, he_count, U);
		down_energy.push_back(temp.energy);
		down_density = temp.density;

		temp = extract_vecs(down_density, he_count, U);
		up_energy.push_back(temp.energy);
		up_density = temp.density;

	}

	eedd_pair return_pair;
	return_pair.density1 = up_density;
	return_pair.density2 = down_density;
	return_pair.energy1  = up_energy;
	return_pair.energy2  = down_energy;

	std::cout << "done with: " << counter << "\n";
	counter++;
	return return_pair;

}

void Model::scan_over_filling(double U, std::string mode, bool reverse) {

	Matrix<double, Dynamic, 1> initial_density;
	std::vector<int> filling_arr = gen_filling_factors(N_x * L, 1, reverse);

	initial_density.resize(N_x, 1);
	initial_density.setZero();

	

	std::string file_name;
	if (reverse) {
		file_name = "f_scan_reverse_" + mode + "_U_=" + std::format("{:.2f}", U) + "_L = " + std::to_string(L) + "_N_x = " + std::to_string(N_x) + ".csv";

	}



	else {
		file_name = "f_scan_" + mode + "_U_=" + std::format("{:.2f}", U) + "_L = " + std::to_string(L) + "_N_x = " + std::to_string(N_x) + ".csv";
	}


	dd_pair sol_pair;

	std::ofstream file(destination + file_name);

	for (int n{ 0 }; n < filling_arr.size(); ++n) {

		if (mode == "even") {
			initial_density.setOnes();
			initial_density *= (double(filling_arr[n]) / (L * N_x));
		}

		sol_pair = hubbard_iteration(filling_arr[n], initial_density, U);
		file << sol_pair.density1.transpose().format(CSVFormat) << ',' << sol_pair.density2.transpose().format(CSVFormat) << ',';

		if (mode == "follow" && (n < filling_arr.size() - 1)) {

			initial_density = sol_pair.density1 * (double(filling_arr[n + 1]) / filling_arr[n]);
		}
	}

	file << U << ',' << N_x << ',' << L << '\n';

	for (int n{ 0 }; n < filling_arr.size(); ++n) {
		file << filling_arr[n];

		if (n < filling_arr.size() - 1) {
			file << ',';
		}
	}


	file.close();

}

void Model::scan_over_filling_and_energy(double U, std::string mode, bool reverse) {

	Matrix<double, Dynamic, 1> initial_density;
	std::vector<int> filling_arr = gen_filling_factors(N_x * L, 1);

	initial_density.resize(N_x, 1);
	initial_density.setZero();



	std::string file_name;
	if (reverse) {
		file_name = "f&E_scan_reverse_" + mode + "_U_=" + std::format("{:.2f}", U) + "_L = " + std::to_string(L) + "_N_x = " + std::to_string(N_x) + ".csv";
	}

	else {
		file_name = "f&E_scan_" + mode + "_U_=" + std::format("{:.2f}", U) + "_L = " + std::to_string(L) + "_N_x = " + std::to_string(N_x) + ".csv";
	}


	eedd_pair sol_pair;

	std::ofstream file(destination + file_name);

	file << std::setprecision(16);

	for (int n{ 0 }; n < filling_arr.size(); ++n) {

		if (mode == "even") {
			initial_density.setOnes();
			initial_density *= (double(filling_arr[n]) / (L * N_x));
		}

		sol_pair = hubbard_iteration_with_energy(filling_arr[n], initial_density, U);

		file << sol_pair.density1.transpose().format(CSVFormat) << ',' << sol_pair.density2.transpose().format(CSVFormat) << ',';

		for (int i{ 0 }; i < sol_pair.energy1.size(); ++i) {
			file << sol_pair.energy2[i] + sol_pair.energy1[i] << ',';
		}
		
		file << U << ',' << N_x << ',' << L << '\n';
			

		if (mode == "follow" && (n < filling_arr.size() - 1)) {

			initial_density = sol_pair.density1 * (double(filling_arr[n + 1]) / filling_arr[n]);
		}
	}

	for (int n{ 0 }; n < filling_arr.size(); ++n) {
		file << filling_arr[n];

		if (n < filling_arr.size() - 1) {
			file << ',';
		}
	}


	file.close();

}

void Model::scan_over_filling_and_U(double start_U, double end_U, int steps_U, std::string mode, bool reverse) {

	Matrix<double, Dynamic, 1> initial_density;
	std::vector<int> filling_arr = gen_filling_factors(N_x * L, 1, reverse);

	initial_density.resize(N_x, 1);
	initial_density.setZero();

	std::string file_name;
	if (reverse) {
		file_name = "f&U_scan_reverse" + mode + "_U=" + std::format("{:.2f}", start_U) + "-" + std::format("{:.2f}", end_U) +
			"-" + std::to_string(steps_U) + "_L=" + std::to_string(L) + "_N_x=" + std::to_string(N_x) + ".csv";
	}
	
	else {
		file_name = "f&U_scan_" + mode + "_U=" + std::format("{:.2f}", start_U) + "-" + std::format("{:.2f}", end_U) +
			"-" + std::to_string(steps_U) + "_L=" + std::to_string(L) + "_N_x=" + std::to_string(N_x) + ".csv";
	}





	std::ofstream file(destination + file_name);

	Matrix<double, 1, Dynamic> array_U;

	array_U.setLinSpaced(steps_U, start_U, end_U);

	dd_pair sol_pair;
	for (int i{ 0 }; i < steps_U; ++i) {
		double U = array_U(i);

		for (int n{ 0 }; n < filling_arr.size(); ++n) {

			//std::cout << filling_arr[n] << "\n";

			if (mode == "even") {
				initial_density.setOnes();
				initial_density *= (double(filling_arr[n]) / (L * N_x));
			}

			sol_pair = hubbard_iteration(filling_arr[n], initial_density, U);
			file << sol_pair.density1.transpose().format(CSVFormat) << ',' << sol_pair.density2.transpose().format(CSVFormat) << ',';

			if (mode == "follow" && (n < filling_arr.size() - 1)) {
				initial_density = sol_pair.density1 * (double(filling_arr[n + 1]) / filling_arr[n]);
			}
		}

		file << U << ',' << N_x << ',' << L << '\n';
	}

	for (int n{ 0 }; n < filling_arr.size(); ++n) {
		file << filling_arr[n];

		if (n < filling_arr.size() - 1) {
			file << ',';
		}
	}

	file.close();

}


