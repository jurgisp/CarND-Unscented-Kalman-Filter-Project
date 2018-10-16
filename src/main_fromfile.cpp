#include <iostream>
#include <math.h>
#include "ukf.h"
#include "tools.h"

using namespace std;


int main() {
  // Create a Kalman Filter instance
  UKF ukf;

  // used to compute the RMSE later
  Tools tools;
  vector<VectorXd> estimations;
  vector<VectorXd> ground_truth;

  std::ifstream infile(
          "/mnt/c/Users/jurgis/Documents/CarND-Extended-Kalman-Filter-Project/data/obj_pose-laser-radar-synthetic-input.txt");

  if (!infile.is_open()) {
    cout << "Can't open file" << endl;
    return -1;
  }

  std::string line;
  auto i = 0;
  while (std::getline(infile, line)) {
    if (i++ > 10)
      break;

    cout << line << endl;
    if (line == "")
      continue;
    std::istringstream iss(line);

    MeasurementPackage meas_package;
    long long timestamp;

    // reads first element from the current line
    string sensor_type;
    iss >> sensor_type;

    if (sensor_type.compare("L") == 0) {
      meas_package.sensor_type_ = MeasurementPackage::LASER;
      meas_package.raw_measurements_ = VectorXd(2);
      float px;
      float py;
      iss >> px;
      iss >> py;
      meas_package.raw_measurements_ << px, py;
      iss >> timestamp;
      meas_package.timestamp_ = timestamp;
    } else if (sensor_type.compare("R") == 0) {

      meas_package.sensor_type_ = MeasurementPackage::RADAR;
      meas_package.raw_measurements_ = VectorXd(3);
      float ro;
      float theta;
      float ro_dot;
      iss >> ro;
      iss >> theta;
      iss >> ro_dot;
      meas_package.raw_measurements_ << ro, theta, ro_dot;
      iss >> timestamp;
      meas_package.timestamp_ = timestamp;
    }
    float x_gt;
    float y_gt;
    float vx_gt;
    float vy_gt;
    iss >> x_gt;
    iss >> y_gt;
    iss >> vx_gt;
    iss >> vy_gt;
    VectorXd gt_values(4);
    gt_values(0) = x_gt;
    gt_values(1) = y_gt;
    gt_values(2) = vx_gt;
    gt_values(3) = vy_gt;
    ground_truth.push_back(gt_values);

    //Call ProcessMeasurment(meas_package) for Kalman filter
    ukf.ProcessMeasurement(meas_package);

    //Push the current estimated x,y positon from the Kalman filter's state vector

    VectorXd estimate(4);

    double p_x = ukf.x_(0);
    double p_y = ukf.x_(1);
    double v = ukf.x_(2);
    double yaw = ukf.x_(3);

    double v1 = cos(yaw) * v;
    double v2 = sin(yaw) * v;

    estimate(0) = p_x;
    estimate(1) = p_y;
    estimate(2) = v1;
    estimate(3) = v2;

    estimations.push_back(estimate);

    // Debug print
    MatrixXd estimate_vs_gt(4, 3);
    estimate_vs_gt << estimate, gt_values, (estimate-gt_values);
    cout << "estimate, truth, diff = " << endl << estimate_vs_gt.transpose() << endl;
  }

  VectorXd rmse = tools.CalculateRMSE(estimations, ground_truth);
  cout << "N = " << estimations.size() << endl;
  cout << "RMSE = " << rmse.transpose() << endl;
}