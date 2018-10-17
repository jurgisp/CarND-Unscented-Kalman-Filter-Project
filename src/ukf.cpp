#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = false;

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 0.2;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.2;
  
  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.

  n_x_ = 5;
  n_aug_ = 7;
  lambda_ = 3 - n_aug_;

  is_initialized_ = false;
  x_ = VectorXd(5);
  P_ = MatrixXd(5, 5);
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {

  if (!is_initialized_) {
    if (use_laser_ && meas_package.sensor_type_ == MeasurementPackage::LASER) {
      ///* state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate]
      x_ = VectorXd(5);
      x_ << meas_package.raw_measurements_, 0, 0, 0;

      P_ = MatrixXd(5, 5);
      // Position uncertainty is given by measurement uncertainty
      P_(0,0) = std_laspx_;
      P_(1,1) = std_laspy_;
      // Use 10. for "maximum uncertainty"
      P_(2,2) = 10.;
      P_(3,3) = 10.;
      P_(4,4) = 10.;

      time_us_ = meas_package.timestamp_;
      is_initialized_ = true;
    }

    if (use_radar_ && meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      double r = meas_package.raw_measurements_(0);
      double phi = meas_package.raw_measurements_(1);
//      ekf_.x_ = VectorXd(4);
//      ekf_.x_ << r*cos(phi), r*sin(phi), 0, 0;
//      ekf_.P_ = MatrixXd(4, 4);
//      ekf_.P_ << 0.09, 0, 0, 0, 0, 0.09, 0, 0, 0, 0, 10., 0, 0, 0, 0, 10.;
//      previous_timestamp_ = measurement_pack.timestamp_;
//      is_initialized_ = true;
    }
  }
  else {
    auto dt = (meas_package.timestamp_ - time_us_) / 1e6;

    Prediction(dt);

    if (use_laser_ && meas_package.sensor_type_ == MeasurementPackage::LASER) {
      UpdateLidar(meas_package);
    }

    if (use_radar_ && meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      UpdateRadar(meas_package);
    }

    time_us_ = meas_package.timestamp_;
  }
}

double NormalizeAngle(double a) {
  while (a> M_PI) a-=2.*M_PI;
  while (a<-M_PI) a+=2.*M_PI;
  return a;
}


/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {

  // x_aug

  VectorXd x_aug = VectorXd(n_aug_);
  x_aug.head(5) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;

  // P_aug

  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
  P_aug.fill(0.0);
  P_aug.topLeftCorner(5,5) = P_;
  P_aug(5,5) = std_a_*std_a_;
  P_aug(6,6) = std_yawdd_*std_yawdd_;

  // Xsig_aug

  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
  MatrixXd L = P_aug.llt().matrixL();
  Xsig_aug.col(0)  = x_aug;
  for (int i = 0; i < n_aug_; i++)
  {
    Xsig_aug.col(i+1)        = x_aug + sqrt(lambda_+n_aug_) * L.col(i);
    Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_+n_aug_) * L.col(i);
  }

  // Xsig_pred

  MatrixXd Xsig_pred = MatrixXd(n_x_, 2 * n_aug_ + 1);
  for (int i = 0; i< 2*n_aug_+1; i++)
  {
    double p_x = Xsig_aug(0,i);
    double p_y = Xsig_aug(1,i);
    double v = Xsig_aug(2,i);
    double yaw = Xsig_aug(3,i);
    double yawd = Xsig_aug(4,i);
    double nu_a = Xsig_aug(5,i);
    double nu_yawdd = Xsig_aug(6,i);

    //predicted state values
    double px_p, py_p;
    if (fabs(yawd) > 0.001) {
      px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
      py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
    }
    else {
      px_p = p_x + v*delta_t*cos(yaw);
      py_p = p_y + v*delta_t*sin(yaw);
    }
    double v_p = v;
    double yaw_p = yaw + yawd*delta_t;
    double yawd_p = yawd;

    //add noise
    px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
    py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
    v_p = v_p + nu_a*delta_t;
    yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
    yawd_p = yawd_p + nu_yawdd*delta_t;

    Xsig_pred(0,i) = px_p;
    Xsig_pred(1,i) = py_p;
    Xsig_pred(2,i) = v_p;
    Xsig_pred(3,i) = yaw_p;
    Xsig_pred(4,i) = yawd_p;
  }

  // Predict x and P

  VectorXd weights = VectorXd(2*n_aug_+1);
  weights.fill(0.5/(n_aug_+lambda_));
  weights(0) = lambda_/(lambda_+n_aug_);

  x_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    x_ = x_ + weights(i) * Xsig_pred.col(i);
  }
  x_(3) = NormalizeAngle(x_(3));

  P_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    VectorXd x_diff = Xsig_pred.col(i) - x_;
    x_diff(3) = NormalizeAngle(x_diff(3));
    P_ = P_ + weights(i) * x_diff * x_diff.transpose() ;
  }
}


/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
}
