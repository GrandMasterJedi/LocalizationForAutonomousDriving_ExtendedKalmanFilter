#include "kalman_filter.h"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;

// Please note that the Eigen library does not initialize 
// VectorXd or MatrixXd objects with zeros upon creation.

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  // *Predict the state

	// state transition prediction, innovation are gaussian with zero mean
	//x_ = F_ * x_  + u;
	x_ = F_ * x_;

	// Covariance P of the states, Q is innovation
	P_ = F_ * P_ * F_.transpose() + Q_;

}

void KalmanFilter::Update(const VectorXd &z) {
  // * update the state by using Kalman Filter equations


	long nx = x_.size();
	MatrixXd I = MatrixXd::Identity(nx, nx);
	MatrixXd Ht = H_.transpose();

	// Measurement update (correction step)
	
	
	VectorXd z_pred = H_ * x_;
	VectorXd y = z - z_pred;
	MatrixXd S = H_*P_*Ht + R_;
	MatrixXd K = P_*Ht*S.inverse();

	
	
	// updated new state
	x_ = x_ + (K*y);
	P_ = (I-K*H_) * P_;


}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  // * update the state by using Extended Kalman Filter equations
	// in EKF the difference is the process innovation is not gaussian, Tailor expansion approximation of Hx
	// double process, position (p) and velocity (v)

	  double px = x_(0);
	  double py = x_(1);
	  double vx = x_(2);
	  double vy = x_(3);

	  // polar coordinates
	  double rho = sqrt(px*px + py*py);  // radial coord from Cartesian 	  
	  double theta = atan2(py, px);		// angular coord from Cartesian
	  double rho_dot = (px*vx + py*vy) / rho; // rate of change
 
	  if (rho < 0.0001) {
        cout << "problem with distance rho = " << rho << endl;
        rho_dot = 0;
      }

  	  //Check angle
	  while ( theta > M_PI || theta < -M_PI ) {
        cout << "problem with theta__ = " << theta << endl;   
		if ( theta > M_PI ) {
		  theta -= M_PI;
		} else {
		  theta += M_PI;
		}
	  }
//   		while (theta < -M_PI)	theta += M_PI;
//   		while (theta > M_PI) theta -= M_PI;

	  VectorXd h = VectorXd(3);
	  h << rho, theta, rho_dot;

  
	  // VectorXd y = z - z_pred;
	  // VectorXd z_pred = H_ * x_;
	  VectorXd y = z - h;
	  
	  // Check prediction value
	  while ( y(1) > M_PI || y(1) < -M_PI ) {
        cout << "problem with theta predicted gap = " << y(1) << endl;   
		if ( y(1) > M_PI ) {
		  y(1) -= M_PI;
		} else {
		  y(1) += M_PI;
		}
	  }

	  MatrixXd S = H_*P_*H_.transpose() + R_;
	  MatrixXd K = P_*H_.transpose()*S.inverse();

	  int nx = x_.size();
	  MatrixXd I = MatrixXd::Identity(nx, nx);
	
	  // updated
	  x_ = x_ + K*y;
	  P_ = (I-K*H_) * P_;

}
