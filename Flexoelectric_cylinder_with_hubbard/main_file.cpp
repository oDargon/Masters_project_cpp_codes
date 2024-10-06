#include <iostream>
#include <Eigen/core>
#include "Constants.h"
#include "Setup.h"
#include <vector>
#include <array>
#include <cmath>
#include <Eigen/Eigenvalues>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <string>
#include <format>



const static IOFormat CSVFormat(StreamPrecision, DontAlignCols, ", ", "\n");
using namespace Eigen;

struct ed_pair {
	double                 energy;
	Matrix<double, N_x, 1> density;
};

struct dd_pair {
	Matrix<double, N_x, 1> density1;
	Matrix<double, N_x, 1> density2;
};

std::vector<int> gen_filling_factors(int end, int step) {

	std::vector<int> return_arr;
	int start = 2;

	while (start <= 2*end) {
		return_arr.push_back(start);
		start += step * 2;
	}

	return return_arr;
}

Matrix<double, N_x, N_x> tk_hubbard(Matrix<double, N_x, 1> rho_vec, int k, double U) {

	static Matrix<double, N_x, N_x> tk_mat;
	static int						tkh_count = 0;

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

	tk_mat.diagonal() += -U * rho_vec;

	return tk_mat;
}

ed_pair extract_vecs(Matrix<double, N_x, 1> rho_vec, int fill_m, double U) {

	static SelfAdjointEigenSolver<Matrix<double, N_x, N_x>> solver;

	static Matrix<int, 1, L* N_x>      indices;
	static Matrix<double, 1, L* N_x>   eigen_vals;
	static Matrix<double, N_x, L* N_x> eigen_vecs;

	for (int i{ 0 }; i < L; ++i) {
		solver.compute(tk_hubbard(rho_vec, i, U));

		eigen_vecs.block<N_x, N_x>(0, N_x * i) = solver.eigenvectors();
		eigen_vals.segment(i * N_x, N_x) = solver.eigenvalues().transpose();
		

	}

	eigen_vecs /= sqrt(L);
	indices.setLinSpaced(L * N_x, 0, L * N_x);
	std::sort(indices.data(), indices.data() + L * N_x, [&](int a, int b) { return eigen_vals(b) > eigen_vals(a); } );

	for (int i = 0; i < L * N_x; ++i) {
		
		if (indices[i] == -1) continue;
		double value = eigen_vals(i);
		Matrix<double, N_x, 1 > value2 = eigen_vecs.col(i);
		int x     = i;
		int y     = indices(i);
		
		while (y != i) {
			indices(x)        = -1;
			eigen_vals(x)     = eigen_vals(y);
			eigen_vecs.col(x) = eigen_vecs.col(y);
			x                 = y;
			y			      = indices(x);
		}
		eigen_vals(x)     = value;
		eigen_vecs.col(x) = value2;
		indices(x)        = -1;

	}

	eigen_vecs = eigen_vecs.array().square();

	ed_pair pair;
	pair.energy  = eigen_vals.head(fill_m).sum()/fill_m;
	pair.density = eigen_vecs.leftCols(fill_m).rowwise().sum();

	return pair;

}

dd_pair hubbard_iteration( int electron_count, Matrix<double, N_x, 1> initial_density, double U ) {

	static int counter = 0;
	int he_count = static_cast<int>(std::round(electron_count/2));

	std::vector<double> up_energy;
	std::vector<double> down_energy;

	Matrix<double, N_x, 1> up_density;
	Matrix<double, N_x, 1> down_density;

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

void scan_over_filling(double U) {

	Matrix<double, N_x, 1> initial_density;
	std::vector<int> filling_arr = gen_filling_factors(N_x*L, 1);

	initial_density << 0, 0, 0, 0, 0, 0;

	dd_pair sol_pair;


	std::string name = "f_scan_U=" + std::format("{:.2f}", U ) + "_L=" + std::to_string(L) + "_N_x=" + std::to_string(N_x) + ".csv";
	std::ofstream file(destination + name );

	for (int n{ 0 }; n < filling_arr.size(); ++n) {
		sol_pair = hubbard_iteration(filling_arr[n] , initial_density, U);
		file << sol_pair.density1.transpose().format(CSVFormat) << ',' << sol_pair.density2.transpose().format(CSVFormat) << ',';
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

void scan_over_filling_and_U(double start_U, double end_U, int steps_U, std::string mode = "zero") {

	Matrix<double, N_x, 1> initial_density;
	std::vector<int> filling_arr = gen_filling_factors(N_x * L, 1);

	initial_density << 0, 0, 0, 0, 0, 0;

	dd_pair sol_pair;
	
	std::string file_name = "f&U_scan_" + mode + "_U=" + std::format("{:.2f}", start_U) + "-" + std::format("{:.2f}", end_U) +
		"-" + std::to_string(steps_U) + "_L=" + std::to_string(L) + "_N_x=" + std::to_string(N_x) + ".csv";
	 
	std::ofstream file(destination + file_name);

	Matrix<double, 1, Dynamic> array_U;

	array_U.setLinSpaced(steps_U, start_U, end_U);

	for (int i{ 0 }; i < steps_U; ++i) {
		double U = array_U(i);

		for (int n{ 0 }; n < filling_arr.size(); ++n) {

			if (mode == "even") {
				initial_density.setOnes();
				initial_density *= (double(filling_arr[n]) / (L * N_x));
			}

			sol_pair = hubbard_iteration(filling_arr[n], initial_density, U);
			file << sol_pair.density1.transpose().format(CSVFormat) << ',' << sol_pair.density2.transpose().format(CSVFormat) << ',';

			if (mode == "follow" && (n < filling_arr.size()-1)) {
				initial_density = sol_pair.density1*(filling_arr[n+1]/ filling_arr[n]);
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


int main() {

	set_up();

	scan_over_filling(2);

	//scan_over_filling_and_U(0, 5, 101, "zero");
	//scan_over_filling_and_U(0, 5, 101, "even");
	//scan_over_filling_and_U(0, 5, 101, "follow");
	//scan_over_filling_and_U(2.5, 3, 101, "follow");

	//scan_over_filling_and_U(-5, 0, 101, "zero");
	//scan_over_filling_and_U(-5, 0, 101, "follow");

	std::cout << "Done\n";

	return 0;
}