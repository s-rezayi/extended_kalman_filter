#include "FusionEKF.h"
#include <iostream>
#include "Eigen/Dense"
#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;
using std::vector;

/**
 * Constructor.
 */
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

  /**
   * TODO: Finish initializing the FusionEKF.
   * TODO: Set the process and measurement noises
   */

  // ekf_.x_ = VectorXd(4);



  H_laser_ << 1.0, 0, 0, 0,
              0, 1.0, 0, 0;

//   Hj_<< 0.5, 0.5, 0, 0,
//         0.5, 0.5, 0, 0,
//         0.5, 0.5, 0.5, 0.5;

  ekf_.F_ = MatrixXd(4, 4);
  ekf_.F_ << 1.0, 0.0, 1.0, 0.0,
             0.0, 1.0, 0.0, 1.0,
             0.0, 0.0, 1.0, 0.0,
             0.0, 0.0, 0.0, 1.0;

  // ekf_.Q_ = MatrixXd(4, 4);
}

/**
 * Destructor.
 */
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
  /**
   * Initialization
   */
  if (!is_initialized_) {
    /**
     * TODO: Initialize the state ekf_.x_ with the first measurement.
     * TODO: Create the covariance matrix.
     * You'll need to convert radar from polar to cartesian coordinates.
     */

    // first measurement
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1.0, 1.0, 1.0, 1.0;

    ekf_.P_ = MatrixXd(4, 4);
    ekf_.P_ << 1.0, 0.0, 0.0, 0.0,
        	     0.0, 1.0, 0.0, 0.0,
        	     0.0, 0.0, 1000.0, 0.0,
        	     0.0, 0.0, 0.0, 1000.0;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      float rho = measurement_pack.raw_measurements_[0];
      float phi = measurement_pack.raw_measurements_[1];
      float rho_dot = measurement_pack.raw_measurements_[2];

      //normalizing phi
      // while(phi > M_PI){
      //   phi -= 2 * M_PI;
      // }
      // while(phi < -M_PI){
      //   phi += 2 * M_PI;
      // }
      
      float px_radar = rho * cos(phi);
      float py_radar = rho * sin(phi);
      float vx_radar = rho_dot * cos(phi);
      float vy_radar = rho_dot * sin(phi);
      
      // if(px_radar < 0.0001)
      // {
      //   px_radar = 0.0001;
      // }
      
      // if(py_radar < 0.0001)
      // {
      //   py_radar = 0.0001;
      // }
      
      // if(vx_radar < 0.0001)
      // {
      //   vx_radar = 0.0001;
      // }
      
      // if(vy_radar < 0.0001)
      // {
      //   vy_radar = 0.0001;
      // }

      ekf_.x_ << px_radar, py_radar, vx_radar, vy_radar;
      previous_timestamp_ = measurement_pack.timestamp_;

    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      float px_laser = measurement_pack.raw_measurements_[0];
      float py_laser = measurement_pack.raw_measurements_[1];
      
      // if(px_laser < 0.0001)
      // {
      //   px_laser = 0.0001;
      // }
      
      // if(py_laser < 0.0001)
      // {
      //   py_laser = 0.0001;
      // }
      
      ekf_.x_ << px_laser, py_laser, 0, 0;
      previous_timestamp_ = measurement_pack.timestamp_;
    }

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /**
   * Prediction
   */

  /**
   * TODO: Update the state transition matrix F according to the new elapsed time.
   * Time is measured in seconds.
   * TODO: Update the process noise covariance matrix.
   * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */

  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
  previous_timestamp_ = measurement_pack.timestamp_;
  
//   if (dt < 0.0001)
//   {
//     dt = 0.0001;
//   }

  float dt_2 = dt * dt;
  float dt_3 = dt_2 * dt;
  float dt_4 = dt_3 * dt;

  // Modify the F matrix so that the time is integrated
  ekf_.F_(0, 2) = dt;
  ekf_.F_(1, 3) = dt;

  // set the process covariance matrix Q
  float noise_ax = 9.0;
  float noise_ay = 9.0;

  ekf_.Q_ = MatrixXd(4, 4);
  ekf_.Q_ <<  dt_4/4*noise_ax, 0, dt_3/2*noise_ax, 0,
         	  0, dt_4/4*noise_ay, 0, dt_3/2*noise_ay,
         	  dt_3/2*noise_ax, 0, dt_2*noise_ax, 0,
         	  0, dt_3/2*noise_ay, 0, dt_2*noise_ay;

  ekf_.Predict();

  /**
   * Update
   */

  /**
   * TODO:
   * - Use the sensor type to perform the update step.
   * - Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    ekf_.H_ = tools.CalculateJacobian(ekf_.x_);
    ekf_.R_ = R_radar_;
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);

  } else {
    ekf_.H_ = H_laser_;
    ekf_.R_ = R_laser_;
    ekf_.Update(measurement_pack.raw_measurements_);

  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
