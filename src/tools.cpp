#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */
	VectorXd rmse(4);
	rmse << 0,0,0,0;


	return rmse;

}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
	MatrixXd H(3,4);

	H << 0,0,0,0,
		0,0,0,0,
		0,0,0,0;


	return H; 
  /**
  TODO:
    * Calculate a Jacobian here.
  */
}
