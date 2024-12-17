#include "constants.h"
#include "model.h"
#include <iostream>
#include <vector>
#include <Eigen/core>

using namespace Eigen;

int main() {


	/*Model dog1{ 4, 100, 1, 3, 1, 100 };

	dog1.scan_over_filling_and_U_E(0, 0, 1, 1000, "follow", false);

	Model dog2{ 5, 100, 1, 3, 1, 100 };

	dog2.scan_over_filling_and_U_E(0, 0, 1, 1000, "follow", false);

	Model dog3{ 6, 100, 1, 3, 1, 100 };

	dog3.scan_over_filling_and_U_E(0, 0, 1, 1000, "follow", false);*/



	//Model dog1{ 6, 100, 1, 3, 1, 100 };

	/*dog1.scan_over_filling(7, 500, "even", false);
	dog1.scan_over_filling(7, 500, "follow", false);*/


	

	//Model dog1{ 6, 100, 1, 10, 1, 100 };

	for (int i = 0; i < 100; ++i) {

		Model dog1{ 6, 100, 1, i*0.3 , 1, 100 };

		dog1.scan_over_filling_and_U_E(0, 0, 1, 1000, "FM", false);

		dog1.scan_over_filling_and_U_E(0, 0, 1, 1000, "AFM", false);

		dog1.scan_over_filling_and_U_E(0, 0, 1, 1000, "polar", false);



		dog1.scan_over_filling_and_U_E(0, 3, 1, 1000, "FM", false);

		dog1.scan_over_filling_and_U_E(0, 3, 1, 1000, "AFM", false);

		dog1.scan_over_filling_and_U_E(0, 3, 1, 1000, "polar", false);



		dog1.scan_over_filling_and_U_E(0, 6, 1, 1000, "FM", false);

		dog1.scan_over_filling_and_U_E(0, 6, 1, 1000, "AFM", false);

		dog1.scan_over_filling_and_U_E(0, 6, 1, 1000, "polar", false);



		dog1.scan_over_filling_and_U_E(0, 9, 1, 1000, "FM", false);

		dog1.scan_over_filling_and_U_E(0, 9, 1, 1000, "AFM", false);

		dog1.scan_over_filling_and_U_E(0, 9, 1, 1000, "polar", false);



	}

	/*Model dog1{ 5, 100, 1, 0, 1, 100 };

	dog1.scan_over_filling_and_U_E(0, 20, 201, 1000, "follow", false);

	dog1.scan_over_filling_and_U_E(0, 20, 201, 1000, "FM", false);

	dog1.scan_over_filling_and_U_E(0, 20, 201, 1000, "AFM", false);

	dog1.scan_over_filling_and_U_E(0, 20, 201, 1000, "polar", false);

	dog1.scan_over_filling_and_U_E(0, 20, 201, 1000, "para", false);*/



	/*Model dog2{ 5, 100, 1, 1, 1, 100 };

	dog2.scan_over_filling_and_U_E(0, 20, 201, 1000, "follow", false);

	dog2.scan_over_filling_and_U_E(0, 20, 201, 1000, "FM", false);

	dog2.scan_over_filling_and_U_E(0, 20, 201, 1000, "AFM", false);

	dog2.scan_over_filling_and_U_E(0, 20, 201, 1000, "polar", false);

	dog2.scan_over_filling_and_U_E(0, 20, 201, 1000, "para", false);



	Model dog3{ 5, 100, 1, 10, 1, 100 };

	dog3.scan_over_filling_and_U_E(0, 20, 201, 1000, "follow", false);

	dog3.scan_over_filling_and_U_E(0, 20, 201, 1000, "FM", false);

	dog3.scan_over_filling_and_U_E(0, 20, 201, 1000, "AFM", false);

	dog3.scan_over_filling_and_U_E(0, 20, 201, 1000, "polar", false);

	dog3.scan_over_filling_and_U_E(0, 20, 201, 1000, "para", false);*/



	/*Model dog2{ 4, 100, 1, 0.1, 1, 100 };

	dog2.scan_over_filling_and_U_E(0, 20, 101, 100, "follow", false);

	dog2.scan_over_filling_and_U_E(0, 20, 101, 100, "FM", false);

	dog2.scan_over_filling_and_U_E(0, 20, 101, 100, "AFM", false);

	dog2.scan_over_filling_and_U_E(0, 20, 101, 100, "polar", false);

	Model dog3{ 4, 100, 1, 1, 1, 100 };

	dog3.scan_over_filling_and_U_E(0, 20, 101, 100, "follow", false);

	dog3.scan_over_filling_and_U_E(0, 20, 101, 100, "FM", false);

	dog3.scan_over_filling_and_U_E(0, 20, 101, 100, "AFM", false);

	dog3.scan_over_filling_and_U_E(0, 20, 101, 100, "polar", false);

	Model dog4{ 4, 100, 1, 10, 1, 100 };

	dog4.scan_over_filling_and_U_E(0, 20, 101, 100, "follow", false);

	dog4.scan_over_filling_and_U_E(0, 20, 101, 100, "FM", false);

	dog4.scan_over_filling_and_U_E(0, 20, 101, 100, "AFM", false);

	dog4.scan_over_filling_and_U_E(0, 20, 101, 100, "polar", false);





	Model dog5{ 5, 100, 1, 0.1, 1, 100 };

	dog5.scan_over_filling_and_U_E(0, 20, 101, 100, "follow", false);

	dog5.scan_over_filling_and_U_E(0, 20, 101, 100, "FM", false);

	dog5.scan_over_filling_and_U_E(0, 20, 101, 100, "AFM", false);

	dog5.scan_over_filling_and_U_E(0, 20, 101, 100, "polar", false);

	Model dog6{ 5, 100, 1, 1, 1, 100 };

	dog6.scan_over_filling_and_U_E(0, 20, 101, 100, "follow", false);

	dog6.scan_over_filling_and_U_E(0, 20, 101, 100, "FM", false);

	dog6.scan_over_filling_and_U_E(0, 20, 101, 100, "AFM", false);

	dog6.scan_over_filling_and_U_E(0, 20, 101, 100, "polar", false);

	Model dog7{ 5, 100, 1, 10, 1, 100 };

	dog7.scan_over_filling_and_U_E(0, 20, 101, 100, "follow", false);

	dog7.scan_over_filling_and_U_E(0, 20, 101, 100, "FM", false);

	dog7.scan_over_filling_and_U_E(0, 20, 101, 100, "AFM", false);

	dog7.scan_over_filling_and_U_E(0, 20, 101, 100, "polar", false);





	Model dog8{ 7, 100, 1, 0.1, 1, 100 };
	
	dog8.scan_over_filling_and_U_E(0, 20, 101, 100, "follow", false);

	dog8.scan_over_filling_and_U_E(0, 20, 101, 100, "FM", false);

	dog8.scan_over_filling_and_U_E(0, 20, 101, 100, "AFM", false);

	dog8.scan_over_filling_and_U_E(0, 20, 101, 100, "polar", false);

	Model dog9{ 7, 100, 1, 1, 1, 100 };

	dog9.scan_over_filling_and_U_E(0, 20, 101, 100, "follow", false);

	dog9.scan_over_filling_and_U_E(0, 20, 101, 100, "FM", false);

	dog9.scan_over_filling_and_U_E(0, 20, 101, 100, "AFM", false);

	dog9.scan_over_filling_and_U_E(0, 20, 101, 100, "polar", false);

	Model dog10{ 7, 100, 1, 10, 1, 100 };

	dog10.scan_over_filling_and_U_E(0, 20, 101, 100, "follow", false);

	dog10.scan_over_filling_and_U_E(0, 20, 101, 100, "FM", false);

	dog10.scan_over_filling_and_U_E(0, 20, 101, 100, "AFM", false);

	dog10.scan_over_filling_and_U_E(0, 20, 101, 100, "polar", false);
	





	Model dog11{ 8, 100, 1, 0.1, 1, 100 };

	dog11.scan_over_filling_and_U_E(0, 20, 101, 100, "follow", false);

	dog11.scan_over_filling_and_U_E(0, 20, 101, 100, "FM", false);

	dog11.scan_over_filling_and_U_E(0, 20, 101, 100, "AFM", false);

	dog11.scan_over_filling_and_U_E(0, 20, 101, 100, "polar", false);

	Model dog12{ 8, 100, 1, 1, 1, 100 };

	dog12.scan_over_filling_and_U_E(0, 20, 101, 100, "follow", false);

	dog12.scan_over_filling_and_U_E(0, 20, 101, 100, "FM", false);

	dog12.scan_over_filling_and_U_E(0, 20, 101, 100, "AFM", false);

	dog12.scan_over_filling_and_U_E(0, 20, 101, 100, "polar", false);

	Model dog13{ 8, 100, 1, 10, 1, 100 };

	dog13.scan_over_filling_and_U_E(0, 20, 101, 100, "follow", false);

	dog13.scan_over_filling_and_U_E(0, 20, 101, 100, "FM", false);

	dog13.scan_over_filling_and_U_E(0, 20, 101, 100, "AFM", false);

	dog13.scan_over_filling_and_U_E(0, 20, 101, 100, "polar", false);*/

	

	/*Model dog2{ 6, 100, 1, 2, 1, 100 };

	dog2.scan_over_filling_and_U_E(0, 10, 101, 500, "AFM", false);

	Model dog3{ 6, 100, 1, 1, 1, 100 };

	dog3.scan_over_filling_and_U_E(0, 10, 101, 500, "AFM", false);

	Model dog4{ 6, 100, 1, 0.1, 1, 100 };

	dog4.scan_over_filling_and_U_E(0, 10, 101, 500, "AFM", false);

	Model dog5{ 6, 100, 1, 0.01, 1, 100 };

	dog5.scan_over_filling_and_U_E(0, 10, 101, 500, "AFM", false);

	Model dog6{ 6, 100, 1, 0.001, 1, 100 };

	dog6.scan_over_filling_and_U_E(0, 10, 101, 500, "AFM", false);*/

	/*Model dog2{ 5, 100, 1, 3, 1, 100 };

	dog2.scan_over_filling_and_U_E(0, 10, 101, 500, "follow", false);*/

	/*Model dog3{ 6, 100, 1, 3, 1, 100 };

	dog3.scan_over_filling_and_U_E(7.5, 7.8, 81, 500, "follow", false);*/

	/*Model dog4{ 7, 100, 1, 3, 1, 100 };

	dog4.scan_over_filling_and_U_E(0, 10, 101, 500, "follow", false);

	Model dog5{ 8, 100, 1, 3, 1, 100 };

	dog5.scan_over_filling_and_U_E(0, 10, 101, 500, "follow", false);*/


	//dog1.scan_over_filling_and_U_E(0, 10, 11, 500, "even", false);

	//dog1.scan_over_filling_and_U(0, 10, 11, 500, "even", false);

	/*Model dog2{ 5, 100, 1, 3, 1, 100 };
	dog2.scan_over_filling_and_U(0, 20, 101, 500, "follow", false);

	Model dog3{ 4, 100, 1, 3, 1, 100 };
	dog3.scan_over_filling_and_U(0, 20, 101, 500, "follow", false);*/


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