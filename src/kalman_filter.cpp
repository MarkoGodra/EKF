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

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in, MatrixXd &Q_in)
{
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
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

  CommonUpdate(y);
}

void KalmanFilter::UpdateEKF(const VectorXd &z)
{
  // h(x')
  VectorXd polar_state = Tools::ConvertCartesianToPolar(x_);
  // Get new measurement error
  VectorXd y = z - polar_state;
  Tools::NormalizeAngle(y(1));

  CommonUpdate(y);
}

void KalmanFilter::CommonUpdate(const Eigen::VectorXd &y)
{
  // Calculation constants
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd K = P_ * Ht * Si;

  //new estimate
  x_ = x_ + (K * y);
  P_ = (I_ - K * H_) * P_;
}