#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

/* 
 * Please note that the Eigen library does not initialize 
 *   VectorXd or MatrixXd objects with zeros upon creation.
 */

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
  /**
   * TODO: predict the state
   */
   x_ = F_ * x_;
   MatrixXd Ft = F_.transpose();
   P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
   * TODO: update the state by using Kalman Filter equations
   */
  const VectorXd y = z - H_ * x_;
  UpdateCommon(y);
}

void KalmanFilter::UpdateCommon(const VectorXd &y)
{
  const MatrixXd Ht = H_.transpose();
  const MatrixXd pht = P_ * Ht;
  const MatrixXd S = H_ * pht + R_;
  const MatrixXd Si = S.inverse();
  const MatrixXd K =  pht * Si;
  // New state
  x_ = x_ + (K * y);
  int x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
   * TODO: update the state by using Extended Kalman Filter equations
   */
   const double rho = sqrt(x_(0)*x_(0)+x_(1)*x_(1));
   const double rhoDot = (x_(0)*x_(2) + x_(1)*x_(3)) / rho;
   const double theta = atan2(x_(1),x_(0));
   VectorXd h = VectorXd(3);
   h << rho, theta, rhoDot;
   VectorXd y = z - h;
  // Normalization
  while (y(1)>M_PI)
    y(1)-=M_PI;
  while (y(1)<-M_PI)
    y(1)+=M_PI;
  UpdateCommon(y);
}
