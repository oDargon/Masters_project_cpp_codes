#include <Eigen/core>
#include <Eigen/Eigenvalues>
#include <vector>
#include "helper_structs.h"

using namespace Eigen;

#if !defined(MYLIB_MODEL_H)
#define MYLIB_MODEL_H 1

class Model {

private:

	int N_x;
	int L;

	double t;
	double alpha;
	double a;
	double R;
    int tkh_count = 0;
	int counter   = 0;


	Matrix<double, 3, 1>                      tube_center;
	Matrix<double, Dynamic, 3>                cell_basis;
	Matrix<double, 1, Dynamic>				  brillouin_zone;
	Matrix<double, 1, Dynamic>				  distance_change;
	Matrix<double, 1, Dynamic>				  cos_arr;


	Matrix<double, Dynamic, Dynamic>                         tk_mat;
	SelfAdjointEigenSolver<MatrixXd>                         solver;
	Matrix<int, 1, Dynamic>                                  indices;
	Matrix<double, 1, Dynamic>                               eigen_vals;
	Matrix<double, Dynamic, Dynamic>                         eigen_vecs;


public:

	Model(int ring_N, int length_N, double hopping_t, double alpha, double spacing_a);

	void setup();
	
	std::vector<int> gen_filling_factors(int end, int step, bool reverse = false);

	void tk_hubbard(Matrix<double, Dynamic, 1> rho_vec, int k, double U);

	ed_pair extract_vecs(Matrix<double, Dynamic, 1> rho_vec, int fill_m, double U);

	dd_pair hubbard_iteration(int electron_count, Matrix<double, Dynamic, 1> initial_density, double U);

	eedd_pair hubbard_iteration_with_energy(int electron_count, Matrix<double, Dynamic, 1> initial_density, double U);



	void scan_over_filling(double U, std::string mode, bool reverse);

	void scan_over_filling_and_energy(double U, std::string mode, bool reverse);

	void scan_over_filling_and_U(double start_U, double end_U, int steps_U, std::string mode = "zero", bool reverse = false);

};




#endif