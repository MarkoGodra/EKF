#ifndef TOOLS_H_
#define TOOLS_H_

#include <vector>
#include "Eigen/Dense"

// TODO: Remove this
#include <iostream>

class Tools
{
public:
  /**
   * Constructor.
   */
  Tools();

  /**
   * Destructor.
   */
  virtual ~Tools();

  /**
   * A helper method to calculate RMSE.
   */
  static Eigen::VectorXd CalculateRMSE(const std::vector<Eigen::VectorXd> &estimations,
                                       const std::vector<Eigen::VectorXd> &ground_truth);

  /**
   * A helper method to calculate Jacobians.
   */
  static Eigen::MatrixXd CalculateJacobian(const Eigen::VectorXd &x_state);

    /**
   * Utility function that converts state from cartesian to polar coordinate sys
   * @param state state matrix
   */
  static Eigen::VectorXd ConvertCartesianToPolar(const Eigen::VectorXd &state);

  /**
   * Utlity function that normalizes angle to {-pi, pi} range
   * @param angle Angle that shall be normalized
   */
  static void NormalizeAngle(float &angle);
};

#endif // TOOLS_H_
