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

	double	t;
	double	alpha;
	double	a;
	double	beta;
	double	R;
	ee_pair base_energy;


    int tkh_count = 0;
	int counter   = 0;
	int max_count = 0;


	Matrix<double, 3, 1>                      tube_center;
	Matrix<double, Dynamic, 3>                cell_basis;
	Matrix<double, 1, Dynamic>				  brillouin_zone;
	Matrix<double, 1, Dynamic>				  distance_change;
	Matrix<double, 1, Dynamic>				  cos_arr;


	Matrix<double, Dynamic, Dynamic>          tk_mat;
	SelfAdjointEigenSolver<MatrixXd>          solver;
	Matrix<double, 1, Dynamic>                eigen_vals;
	Matrix<double, Dynamic, Dynamic>          eigen_vecs;


public:

	Model(int ring_N, int length_N, double hopping_t, double alpha, double spacing_a, double beta);

	void setup();

	double fermi_dirac(const double E, const double mu);
	
	std::vector<int> gen_filling_factors(int end, int step, bool reverse = false);

	void tk_hubbard(Matrix<double, Dynamic, 1> rho_vec, int k, double U);

	void extract_base_energy();

	ed_pair extract_vecs(Matrix<double, Dynamic, 1> rho_vec, double fill_mu, double U);

	edd_pair hubbard_iteration(double mu, dd_pair initial_density, double U);

	void scan_over_filling(double U, int mu_granularity, std::string mode, bool reverse);

	void scan_over_filling_and_U(double start_U, double end_U, int steps_U, int mu_granularity, std::string mode, bool reverse);

	void scan_over_filling_and_U_E(double start_U, double end_U, int steps_U, int mu_granularity, std::string mode, bool reverse);

};




#endif