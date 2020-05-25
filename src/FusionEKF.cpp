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
FusionEKF::FusionEKF()
{
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  // Initializing P
  ekf_.P_ = MatrixXd(4, 4);
  ekf_.P_ << 1, 0, 0, 0,
      0, 1, 0, 0,
      0, 0, 1000, 0,
      0, 0, 0, 1000;

  H_laser_ << 1, 0, 0, 0,
      0, 1, 0, 0;

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
      0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
      0, 0.0009, 0,
      0, 0, 0.09;

  // Measurement noise
  noise_ax_ = 9.0;
  noise_ay_ = 9.0;
}

/**
 * Destructor.
 */
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack)
{
  if (!is_initialized_)
  {
    // Init
    // first measurement
    VectorXd x = VectorXd(4);

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR)
    {
      // Read out the Radara measurement components
      float rho = measurement_pack.raw_measurements_[0];
      float phi = measurement_pack.raw_measurements_[1];

      float px = rho * cos(phi);
      float py = rho * sin(phi);

      // Do not use rho dot in init phase since
      // radar measurement does not contain enough information to determin
      // the state variable velocities vx, vy
      x << px,
          py,
          0.0,
          0.0;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER)
    {
      x << measurement_pack.raw_measurements_[0],
          measurement_pack.raw_measurements_[1],
          0.0,
          0.0;
    }

    previous_timestamp_ = measurement_pack.timestamp_;

    MatrixXd F(4, 4);
    // delta t is zero, since this is first measurement
    F << 1.0, 0.0, 0.0, 0.0,
        0.0, 1.0, 0.0, 0.0,
        0.0, 0.0, 1.0, 0.0,
        0.0, 0.0, 0.0, 1.0;

    MatrixXd P(4, 4);
    // We are very uncertain about velocity
    P << 1.0, 0.0, 0.0, 0.0,
        0.0, 1.0, 0.0, 0.0,
        0.0, 0.0, 1000.0, 0.0,
        0.0, 0.0, 0.0, 1000.0;

    MatrixXd Q(4, 4);
    Q << 0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0;

    ekf_.Init(x, P, F, Q);

    // done initializing, no need to predict or update
    is_initialized_ = true;

    return;
  }

  // Prediction
  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
  previous_timestamp_ = measurement_pack.timestamp_;

  ekf_.F_(0, 2) = dt;
  ekf_.F_(1, 3) = dt;

  float dt_2 = dt * dt;
  float dt_3 = dt_2 * dt;
  float dt_4 = dt_3 * dt;
  float dt_4_div_4 = dt_4 / 4.0;
  float dt_3_div_2 = dt_3 / 2.0;

  ekf_.Q_ = MatrixXd(4, 4);
  ekf_.Q_ << dt_4_div_4 * noise_ax_, 0.0, dt_3_div_2 * noise_ax_, 0.0,
      0.0, dt_4_div_4 * noise_ay_, 0.0, dt_3_div_2 * noise_ay_,
      dt_3_div_2 * noise_ax_, 0.0, dt_2 * noise_ax_, 0.0,
      0.0, dt_3_div_2 * noise_ay_, 0.0, dt_2 * noise_ay_;

  ekf_.Predict();

  // Update
  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR)
  {
    ekf_.R_ = R_radar_;
    ekf_.H_ = Tools::CalculateJacobian(ekf_.x_);

    // Do extended kalman filter update here (Non-linear)
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  }
  else
  {
    ekf_.R_ = R_laser_;
    ekf_.H_ = H_laser_;

    // Do kalman filter update here (Linear)
    ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  std::cout << "x_ = " << ekf_.x_ << std::endl;
  std::cout << "P_ = " << ekf_.P_ << std::endl;
}
