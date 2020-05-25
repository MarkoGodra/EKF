#include "tools.h"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth)
{
   VectorXd rmse(4);
   rmse << 0, 0, 0, 0;

   if (estimations.size() != ground_truth.size() || estimations.size() == 0)
   {
      std::cout << "Invalid estimation or ground_truth data" << std::endl;
      return rmse;
   }

   // accumulate squared residuals
   for (unsigned int i = 0; i < estimations.size(); ++i)
   {

      VectorXd residual = estimations[i] - ground_truth[i];

      // coefficient-wise multiplication
      residual = residual.array() * residual.array();
      rmse += residual;
   }

   // calculate the mean
   rmse = rmse / estimations.size();

   // calculate the squared root
   rmse = rmse.array().sqrt();

   // return the result
   return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd &x_state)
{
   MatrixXd Hj(3, 4);
   // recover state parameters
   double px = x_state(0);
   double py = x_state(1);
   double vx = x_state(2);
   double vy = x_state(3);

   // pre-compute a set of terms to avoid repeated calculation
   double c1 = px * px + py * py;
   double c2 = sqrt(c1);
   double c3 = (c1 * c2);

   // check division by zero
   if (fabs(c1) < 0.0001)
   {
      std::cout << "CalculateJacobian () - Error - Division by Zero" << std::endl;
      return Hj;
   }

   // compute the Jacobian matrix
   Hj << (px / c2), (py / c2), 0, 0,
       -(py / c1), (px / c1), 0, 0,
       py * (vx * py - vy * px) / c3, px * (px * vy - py * vx) / c3, px / c2, py / c2;

   return Hj;
}

// Actually represents h(x')
VectorXd Tools::ConvertCartesianToPolar(const VectorXd &state)
{
   const double px = state(0);
   const double py = state(1);
   const double vx = state(2);
   const double vy = state(3);

   double rho;
   double phi;
   double rho_dot;

   rho = sqrt(px * px + py * py);
   phi = std::atan2(py, px);

   // Avoid zero divison
   // Rho cannot be negative value
   if (rho < 0.000001)
   {
      rho = 0.000001;
   }

   rho_dot = (px * vx + py * vy) / rho;

   VectorXd return_value = VectorXd(3);
   return_value << rho, phi, rho_dot;

   return return_value;
}

void Tools::NormalizeAngle(double &angle)
{
   while (angle > M_PI || angle < -M_PI)
   {
      if (angle > M_PI)
      {
         angle -= M_PI;
      }
      else
      {
         angle += M_PI;
      }
   }
}
