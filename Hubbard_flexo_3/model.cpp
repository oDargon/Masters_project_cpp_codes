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

Model::Model(int ring_N = 6, int length_N = 100, double hopping_t = 1, double alpha = 3, double spacing_a = 1, double beta = 100)
	: N_x{ ring_N }, L{ length_N }, t{ hopping_t }, alpha{ alpha }, a{ spacing_a }, beta{ beta }, R{ L * a / (2 * pi) } {

	setup();
}

double Model::fermi_dirac(const double E, const double mu) {
    return 1 / (1 + exp((E - mu) * beta));
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
		tk_mat(j, j) = -2 * t * (1 - alpha * distance_change(j)) * cos_arr(k);
	}

	tk_mat.diagonal() += U * rho_vec;

	tkh_count++;

}

ed_pair Model::extract_vecs(Matrix<double, Dynamic, 1> rho_vec, double fill_mu, double U) {

	for (int i{ 0 }; i < L; ++i) {
		tk_hubbard(rho_vec, i, U);
		solver.compute(tk_mat);

		eigen_vecs.block(0, N_x * i, N_x, N_x) = solver.eigenvectors();
		eigen_vals.segment(i * N_x, N_x) = solver.eigenvalues().transpose();

	}

	//std::cout << eigen_vals;


	eigen_vecs = eigen_vecs.array().square();

	double                     new_energy = 0;
	Matrix<double, Dynamic, 1> new_density;
	new_density.resize(N_x, 1);
	new_density.setZero();
	
	new_density.resize(N_x, 1);

	for (int i{ 0 }; i < L * N_x; ++i) {
		new_energy  += fermi_dirac(eigen_vals(i), fill_mu) * (eigen_vals(i) - fill_mu);
		new_density += fermi_dirac(eigen_vals(i), fill_mu) * eigen_vecs.col(i);
	}

	ed_pair pair;

	pair.energy  = new_energy / L;
	pair.density = new_density / L;

	return pair;

}

void Model::extract_base_energy() {
	
	Matrix<double, Dynamic, 1> rho_vec;
	rho_vec.setLinSpaced(N_x, 1, 10);

	for (int i{ 0 }; i < L; ++i) {
		tk_hubbard(rho_vec, i, 0);
		solver.compute(tk_mat);

		eigen_vecs.block(0, N_x * i, N_x, N_x) = solver.eigenvectors();
		eigen_vals.segment(i * N_x, N_x) = solver.eigenvalues().transpose();

	}

	base_energy.energy1 = eigen_vals.minCoeff();
	base_energy.energy2 = eigen_vals.maxCoeff();

}

edd_pair Model::hubbard_iteration(double mu, dd_pair initial_density, double U) {

	std::vector<double> up_energy;
	std::vector<double> down_energy;
	std::vector<double> tot_energy;

	Matrix<double, Dynamic, 1> up_density;
	Matrix<double, Dynamic, 1> down_density;
	Matrix<double, Dynamic, 1> new_up_density;
	Matrix<double, Dynamic, 1> new_down_density;

	up_density = initial_density.density1;
	down_density = initial_density.density2;

	ed_pair temp;


	int mini_count = 0;
	double differ = 0;

	/*while (mini_count < 2 || (abs((1-(tot_energy[tot_energy.size()-2]/ tot_energy.back()))) > epsilon && mini_count <= MAX_ITERS ) ) {*/
	/*for (int i{ 0 }; i < 10; ++i ) {*/
	while (mini_count < 2 || (differ > epsilon && mini_count <= MAX_ITERS)) {

		//std::cout << up_density.transpose() << "; " << down_density.transpose() << "\n";

		temp = extract_vecs(down_density, mu, U);
		up_energy.push_back(temp.energy);
		new_up_density = temp.density;

		temp = extract_vecs(up_density, mu, U);
		down_energy.push_back(temp.energy);
		new_down_density = temp.density;

		double i_correction_0 = U * (new_up_density.array() * new_down_density.array()).sum();
		double i_correction_1 = U * (up_density.array() * new_down_density.array()).sum();
		double i_correction_2 = U * (new_up_density.array() * down_density.array()).sum();

		tot_energy.push_back(up_energy.back() + down_energy.back() + i_correction_0 - i_correction_1 - i_correction_2);


		differ = (new_up_density - up_density).cwiseAbs().sum() + (new_down_density - down_density).cwiseAbs().sum();

		up_density = new_up_density * mix + up_density * (1 - mix);
		down_density = new_down_density * mix + down_density * (1 - mix);

		++mini_count;
	}

	edd_pair return_pair;
	return_pair.density1   = up_density;
	return_pair.density2   = down_density;
	return_pair.energy_tot = tot_energy.back();

	std::cout << "done with: " << counter << " mini count: " << mini_count << "\n";
	counter++;

	if (mini_count == MAX_ITERS + 1) {
		++max_count;
	}

	return return_pair;

}

void Model::scan_over_filling(double U, int mu_granularity, std::string mode, bool reverse) {
	
	Matrix<double, Dynamic, 1> mu_arr;
	mu_arr.setLinSpaced(mu_granularity, base_energy.energy1, base_energy.energy2 + U);

	Matrix<double, Dynamic, 1> initial_density_part;
	initial_density_part.resize(N_x, 1);
	initial_density_part.setZero();


	dd_pair start_pair;
	start_pair.density1 = initial_density_part;
	start_pair.density2 = initial_density_part;
	
	
	std::string file_name;
	if (reverse) {
		file_name = "f_scan_reverse_" + mode + "_U_=" + std::format("{:.2f}", U) + "_L = " + std::to_string(L) + "_N_x = " + std::to_string(N_x) + ".csv";

	}
	else {
		file_name = "f_scan_" + mode + "_U_=" + std::format("{:.2f}", U) + "_L = " + std::to_string(L) + "_N_x = " + std::to_string(N_x) + ".csv";
	}


	edd_pair sol_pair;

	std::ofstream file(destination + file_name);

	for (int n{ 0 }; n < mu_arr.size(); ++n) {

		if (mode == "even") {
			start_pair.density1.setOnes();
			start_pair.density2.setOnes();

			start_pair.density1 *= 0.51;		
		    start_pair.density2 *= 0.49;
		}

		sol_pair = hubbard_iteration(mu_arr( n ), start_pair, U);
		file << sol_pair.density1.transpose().format(CSVFormat) << ',' << sol_pair.density2.transpose().format(CSVFormat) << ',';

		if (mode == "follow" && (n < mu_arr.size() - 1)) {

			start_pair.density1 = sol_pair.density1;
			start_pair.density2 = sol_pair.density2;
		}
	}

	file << U << ',' << N_x << ',' << L << '\n';

	for (int n{ 0 }; n < mu_arr.size(); ++n) {
		file << mu_arr(n);

		if (n < mu_arr.size() - 1) {
			file << ',';
		}
	}


	file.close();

}

void Model::scan_over_filling_and_U(double start_U, double end_U, int steps_U, int mu_granularity, std::string mode, bool reverse) {

	Matrix<double, 1, Dynamic> array_U;

	array_U.setLinSpaced(steps_U, start_U, end_U);

	Matrix<double, Dynamic, 1> initial_density_part;
	initial_density_part.resize(N_x, 1);
	initial_density_part.setZero();

	dd_pair start_pair;
	start_pair.density1 = initial_density_part;
	start_pair.density2 = initial_density_part;


	std::string file_name;
	if (reverse) {
		file_name = "f&U_scan_reverse" + mode + "_U=" + std::format("{:.2f}", start_U) + "-" + std::format("{:.2f}", end_U) +
			"-" + std::to_string(steps_U) + "_L=" + std::to_string(L) + "_N_x=" + std::to_string(N_x) + ".csv";
	}

	else {
		file_name = "f&U_scan_" + mode + "_U=" + std::format("{:.2f}", start_U) + "-" + std::format("{:.2f}", end_U) +
			"-" + std::to_string(steps_U) + "_L=" + std::to_string(L) + "_N_x=" + std::to_string(N_x) + ".csv";
	}


	edd_pair sol_pair;

	std::ofstream file(destination + file_name);

	

	for (int i{ 0 }; i < steps_U; ++i) {
		double U = array_U(i);

		Matrix<double, Dynamic, 1> mu_arr;
		mu_arr.setLinSpaced(mu_granularity, base_energy.energy1, base_energy.energy2 + U);

		if (reverse) {
			mu_arr.reverse();
		}

		for (int n{ 0 }; n < mu_arr.size(); ++n) {

			if (mode == "even") {
				start_pair.density1.setOnes();
				start_pair.density2.setOnes();

				start_pair.density1 *= 0.51;
				start_pair.density2 *= 0.49;
			}

			sol_pair = hubbard_iteration(mu_arr(n), start_pair, U);

			file << sol_pair.density1.transpose().format(CSVFormat) << ',' << sol_pair.density2.transpose().format(CSVFormat) << ',';

			if (mode == "follow" && (n < mu_arr.size() - 1)) {

				start_pair.density1 = sol_pair.density1*0.99;
				start_pair.density2 = sol_pair.density2*1.01;

				//std::cout << sol_pair.density1.transpose() << ", " << sol_pair.density2.transpose() << "\n";

			}
		}

		file << U << ',' << N_x << ',' << L << '\n';

		for (int n{ 0 }; n < mu_arr.size(); ++n) {
			file << mu_arr(n);

			if (n < mu_arr.size() - 1) {
				file << ',';
			}
		}

		file << '\n';
	}


	file.close();

	std::cout << "Number of hits on maxiters: " << max_count << "\n";

}

void Model::scan_over_filling_and_U_E(double start_U, double end_U, int steps_U, int mu_granularity, std::string mode, bool reverse) {

	Matrix<double, 1, Dynamic> array_U;

	array_U.setLinSpaced(steps_U, start_U, end_U);

	Matrix<double, Dynamic, 1> initial_density_part;
	initial_density_part.resize(N_x, 1);
	initial_density_part.setZero();

	dd_pair start_pair;
	start_pair.density1 = initial_density_part;
	start_pair.density2 = initial_density_part;


	std::string file_name;
	if (reverse) {
		file_name = "f&U_scan_reverse" + mode + "_U=" + std::format("{:.2f}", start_U) + "-" + std::format("{:.2f}", end_U) +
			"-" + std::to_string(steps_U) + "_L=" + std::to_string(L) + "_N_x=" + std::to_string(N_x) + "_" + std::to_string(MAX_ITERS) + ".csv";
	}

	else {
		file_name = "f&U_scan_" + mode + "_U=" + std::format("{:.2f}", start_U) + "-" + std::format("{:.2f}", end_U) +
			"-" + std::to_string(steps_U) + "_L=" + std::to_string(L) + "_N_x=" + std::to_string(N_x) + "_Al=" + std::format("{:.2f}", alpha) + "_" + std::to_string(MAX_ITERS) + ".csv";
	}


	edd_pair sol_pair;

	std::ofstream file(destination + file_name);

	std::ofstream file_E(destination + "Energy_"  + file_name  );



	for (int i{ 0 }; i < steps_U; ++i) {
		double U = array_U(i);

		Matrix<double, Dynamic, 1> mu_arr;
		mu_arr.setLinSpaced(mu_granularity, base_energy.energy1, base_energy.energy2 + U);

		if (reverse) {
			mu_arr.reverseInPlace();	
		}

		for (int n{ 0 }; n < mu_arr.size(); ++n) {

			sol_pair = hubbard_iteration(mu_arr(n), start_pair, U);

			file << sol_pair.density1.transpose().format(CSVFormat) << ',' << sol_pair.density2.transpose().format(CSVFormat) << ',';

			file_E << sol_pair.energy_tot << ',';

			if (mode == "even") {
				start_pair.density1.setOnes();
				start_pair.density2.setOnes();

				double avg = (sol_pair.density1.mean() + sol_pair.density2.mean()) / 2;


				start_pair.density1 *= (1.01*(avg));
				start_pair.density2 *= (0.99)*(avg);

				/*std::cout << start_pair.density1 << "\n";*/
			}

			if (mode == "follow" && (n < mu_arr.size() - 1)) {

				start_pair.density1 = sol_pair.density1 * 0.99;
				start_pair.density2 = sol_pair.density2 * 1.01;

				//std::cout << sol_pair.density1.transpose() << ", " << sol_pair.density2.transpose() << "\n";

			}

			if (mode == "para") {
				start_pair.density1.setOnes();
				start_pair.density2.setOnes();

				double avg = (sol_pair.density1.mean() + sol_pair.density2.mean()) / 2;

				start_pair.density1 *= (avg);
				start_pair.density2 *= (avg);

			}

			if (mode == "FM") {
				start_pair.density1.setOnes();
				start_pair.density2.setOnes();

				double avg = (sol_pair.density1.mean() + sol_pair.density2.mean()) / 2;

				start_pair.density1 *= (1.99 * (avg));
				start_pair.density2 *= (0.01) * (avg);

			}

			if (mode == "AFM") {
				start_pair.density1.setOnes();
				start_pair.density2.setOnes();

				double avg = (sol_pair.density1.mean() + sol_pair.density2.mean()) / 2;

				for (int q{0}; q < N_x; ++q) {
					if (q % 2 == 0) {
						start_pair.density1[q] *= (1.99 * (avg));
						start_pair.density2[q] *= (0.01) * (avg);
					}
					else {
						start_pair.density1[q] *= (0.01 * (avg));
						start_pair.density2[q] *= (1.99) * (avg);
					}
				}

			}

			if (mode == "AFM2") {
				start_pair.density1.setOnes();
				start_pair.density2.setOnes();

				double avg = (sol_pair.density1.mean() + sol_pair.density2.mean()) / 2;

				for (int q{ 0 }; q < N_x; ++q) {
					if (q % 2 == 0) {
						start_pair.density1[q] *= (0.01 * (avg));
						start_pair.density2[q] *= (1.99) * (avg);
					}
					else {
						start_pair.density1[q] *= (1.99 * (avg));
						start_pair.density2[q] *= (0.01) * (avg);
					}
				}

			}

			if (mode == "polar") {
				start_pair.density1.setOnes();
				start_pair.density2.setOnes();

				double avg = (sol_pair.density1.mean() + sol_pair.density2.mean()) / 2;
				
				for (int q{ 0 }; q < N_x; ++q) {
					
					if (q == 0) {
						start_pair.density1[q] = avg * 1.99;
						start_pair.density2[q] = avg * 0.01;
					}
					else {
						start_pair.density1[q] = avg * 0.01;
						start_pair.density2[q] = avg * 1.99;
					}
				}

			}
		
			if (mode == "polar2") {
				start_pair.density1.setOnes();
				start_pair.density2.setOnes();

				double avg = (sol_pair.density1.mean() + sol_pair.density2.mean()) / 2;

				for (int q{ 0 }; q < N_x; ++q) {

					if (q == 0) {
						start_pair.density1[q] = avg * 0.01;
						start_pair.density2[q] = avg * 1.99;
					}
					else {
						start_pair.density1[q] = avg * 1.99;
						start_pair.density2[q] = avg * 0.01;
					}
				}

			}
		
		
		}

		file << U << ',' << N_x << ',' << L << '\n';

		file_E << U << '\n';

		for (int n{ 0 }; n < mu_arr.size(); ++n) {
			file << mu_arr(n);

			file_E << mu_arr(n);

			if (n < mu_arr.size() - 1) {
				file << ',';
				file_E << ',';
			}
		}

		file << '\n';
		file_E << '\n';
	}


	file.close();

	std::cout << "Number of hits on maxiters: " << max_count << "\n";
}