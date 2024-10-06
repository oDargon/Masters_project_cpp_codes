#include "constants.h"
#include "model.h"
#include <iostream>
#include <vector>
#include <Eigen/core>

using namespace Eigen;

int main() {

	Model dog{ 6, 100, 1, 3, 1 };

	//dog.scan_over_filling(2, "zero", false);
	//dog.scan_over_filling(2, "zero", true);

	/*dog.scan_over_filling(2.5, "follow", false);
	dog.scan_over_filling(2.5, "follow", true);*/

	//dog.scan_over_filling_and_U(1, 2, 3, "follow", false);
	//dog.scan_over_filling_and_U(1, 2, 3, "follow", true);

	//dog.scan_over_filling(2, "follow", true);
	//dog.scan_over_filling(2, "zero", true);
	//dog.scan_over_filling(2, "even", true);

	//dog.scan_over_filling_and_energy(2, "follow", false);

	/*dog.scan_over_filling(2, "zero", false);
	dog.scan_over_filling(2, "follow", false);*/

	//Model dog1{ 6, 25, 1, 3, 1 };
	//Model dog2{ 6, 50, 1, 3, 1 };
	//Model dog3{ 6, 75, 1, 3, 1 };
	//Model dog4{ 6, 100, 1, 3, 1 };
	//Model dog5{ 6, 125, 1, 3, 1 };
	//Model dog6{ 6, 150, 1, 3, 1 };
	//Model dog7{ 6, 175, 1, 3, 1 };
	//Model dog8{ 6, 200, 1, 3, 1 };
	//Model dog9{ 6, 225, 1, 3, 1 };

	/*dog1.scan_over_filling(2, "follow", false);
	dog2.scan_over_filling(2, "follow", false);
	dog3.scan_over_filling(2, "follow", false);
	dog4.scan_over_filling(2, "follow", false);
	dog5.scan_over_filling(2, "follow", false);
	dog6.scan_over_filling(2, "follow", false);*/
	//dog7.scan_over_filling(2, "follow", false);
	//dog8.scan_over_filling(2, "follow", false);
	//dog9.scan_over_filling(2, "follow", false);

	/*dog.scan_over_filling_and_energy(2, "follow", false);*/

	/*dog.scan_over_filling(2, "zero", false);
	dog.scan_over_filling(2, "follow", false);*/


	//dog.scan_over_filling_and_U(0, 5, 20, "zero", false);
	/*dog.scan_over_filling_and_U(4,6,21, "follow", false);*/

	dog.scan_over_filling(2, "zero", false);
	dog.scan_over_filling(2, "follow", false);

	return 0;
}