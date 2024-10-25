#include "constants.h"
#include "model.h"
#include <iostream>
#include <vector>
#include <Eigen/core>

using namespace Eigen;

int main() {

	//Model dog1{ 6, 100, 1, 3, 1, 100 };

	/*dog1.scan_over_filling(7, 500, "even", false);
	dog1.scan_over_filling(7, 500, "follow", false);*/


	Model dog1{ 6, 100, 1, 3, 1, 100 };
	dog1.scan_over_filling_and_U(0, 10, 101, 500, "follow", false);

	/*Model dog2{ 5, 100, 1, 3, 1, 100 };
	dog2.scan_over_filling_and_U(0, 10, 101, 500, "even", false);

	Model dog3{ 4, 100, 1, 3, 1, 100 };
	dog3.scan_over_filling_and_U(0, 10, 101, 500, "even", false);*/


	/*Matrix<double, Dynamic, 1> new_density;
	new_density.resize(6, 1); 
	new_density.setOnes();
	
	dd_pair start_pair;
	start_pair.density1 = new_density*0.51;
	start_pair.density2 = new_density*0.49;

	dog.hubbard_iteration(3, start_pair, 0);*/


	/*dog.scan_over_filling(3.5, 600, "even", false);*/
	//dog.scan_over_filling_and_U(0,10,101, 600, "even", false);

	/*dog.scan_over_filling(2, 600, "zero", false);
	dog.scan_over_filling(2, 600, "follow", false);*/

	return 0;
}