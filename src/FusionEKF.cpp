#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

// Constructor.

FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0, 
			 0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
	          0, 0.0009, 0,
			  0, 0, 0.09;

  // **DO:* Finish initializing the FusionEKF. * Set the process and measurement noises
  //
  //// Processes covariance
  //ekf_.P_ = MatrixXd(4, 4);
  //ekf_.P_ << 1, 0, 0, 0,
		//	0, 1, 0, 0,
		//	0, 0, 1000, 0,
		//	0, 0, 0, 1000;

  H_laser_ << 1.0, 0.0, 0.0, 0.0,
              0.0, 1.0, 0.0, 0.0;

  Hj_ << 1.0, 0.0, 0.0, 0.0,
         1.0, 1.0, 0.0, 0.0,
         1.0, 1.0, 1.0, 1.0;

  
}

//Destructor.
FusionEKF::~FusionEKF() {}


void FusionEKF::ProcessMeasurement(const MeasurementPackage &mp) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    /**DO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;

    if (mp.sensor_type_ == MeasurementPackage::RADAR) {
      // Convert radar from polar to cartesian coordinates and initialize state.
		const double rho = mp.raw_measurements_[0]; // Range - radial distance from origin in polar coordinate
		const double phi = mp.raw_measurements_[1]; // Bearing - angle between rho and axis
		const double rho_dot = mp.raw_measurements_[2]; // Radial Velocity - rho rate

		double c1 =  rho * cos(phi) < 1e-6 ? 1e-6 : rho * cos(phi);
		double c2 =  rho * sin(phi) < 1e-6 ? 1e-6 : rho * sin(phi);
		
		double v1 = rho_dot * cos(phi);
  		double v2 = rho_dot * sin(phi);

		ekf_.x_ << c1, c2, v1 , v2;
		
    }
    else if (mp.sensor_type_ == MeasurementPackage::LASER) {
      // Initialize state.
		// The coordinates are cartesian
		ekf_.x_ <<  mp.raw_measurements_[0], mp.raw_measurements_[1], 0.0, 0.0;
      
    }

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/
  /* DO:
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */

  
  double dt = (mp.timestamp_ - previous_timestamp_) / 1e6;
  previous_timestamp_ = mp.timestamp_;

  if (dt<1e-9) { return; }  // simultaneous measurement: no not predict

  // State transition matrix update
  ekf_.F_ = MatrixXd(4, 4);
  ekf_.F_ << 1, 0, dt, 0,
             0, 1, 0, dt,
             0, 0, 1, 0,
             0, 0, 0, 1;

  // Noise covariance matrix computation
  // Set suggested noise value
  double sig1 = 9.0;
  double sig2 = 9.0;
  
  double dt2 = dt * dt; //dt^2
  double dt3 = dt2 * dt; //dt^3
  double dt4 = dt3 * dt; //dt^4
  
  double dt4_4 = dt4 / 4; // dt^4 / 4
  double dt3_2 = dt3 / 2; // dt^3 / 2
  ekf_.Q_ = MatrixXd(4, 4);
  ekf_.Q_ << dt4_4*sig1,	0,				dt3_2 * sig1,	0,
	         0,				dt4_4*sig2,		0,				dt3_2*sig2,
	         dt3_2*sig1,	0,				dt2*sig1,		0,
 	         0,				dt3_2*sig2,		0,				dt2*sig2;

  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /*DO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */


  if (mp.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
	  ekf_.H_ = tools.CalculateJacobian(ekf_.x_); 
	  ekf_.R_ = R_radar_;
	  ekf_.UpdateEKF(mp.raw_measurements_);
  } else {
    // Laser updates
	  ekf_.H_ = H_laser_;
	  ekf_.R_ = R_laser_;
	  ekf_.Update(mp.raw_measurements_);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
