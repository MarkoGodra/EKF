#include "kalman_filter.h"
#include "tools.h"
#include <math.h>

using Eigen::MatrixXd;
using Eigen::VectorXd;

/* 
 * Please note that the Eigen library does not initialize 
 *   VectorXd or MatrixXd objects with zeros upon creation.
 */

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_lidar_in, MatrixXd &R_radar_in, MatrixXd &Q_in)
{
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_lidar_ = R_lidar_in;
  R_radar_ = R_radar_in;
  Q_ = Q_in;
  I_ = MatrixXd::Identity(x_in.size(), x_in.size());
}

void KalmanFilter::Predict()
{
  x_ = F_ * x_;
  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z)
{
  VectorXd z_pred = H_ * x_;
  VectorXd y = z - z_pred;
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_lidar_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

  //new estimate
  x_ = x_ + (K * y);
  P_ = (I_ - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z)
{
  // Calculation constants
  MatrixXd Hj = Tools::CalculateJacobian(x_);
  MatrixXd Ht = Hj.transpose();
  MatrixXd S = Hj * P_ * Ht + R_radar_;
  MatrixXd Si = S.inverse();
  MatrixXd K = P_ * Ht * Si;

  // h(x')
  VectorXd polar_state = Tools::ConvertCartesianToPolar(x_);
  // Get new measurement error
  VectorXd y = z - polar_state;

  //new estimate
  x_ = x_ + (K * y);
  P_ = (I_ - K * Hj) * P_;
}