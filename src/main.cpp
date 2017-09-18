#include <fstream>
#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Dense"
#include "Eigen-3.3/Eigen/QR"
#include "json.hpp"
#include "spline.h"

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;

// for convenience
using json = nlohmann::json;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.find_first_of("}");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

double distance(double x1, double y1, double x2, double y2)
{
    return sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
}
int ClosestWaypoint(double x, double y, vector<double> maps_x, vector<double> maps_y)
{

    double closestLen = 100000; //large number
    int closestWaypoint = 0;

    for(int i = 0; i < maps_x.size(); i++)
    {
        double map_x = maps_x[i];
        double map_y = maps_y[i];
        double dist = distance(x,y,map_x,map_y);
        if(dist < closestLen)
        {
            closestLen = dist;
            closestWaypoint = i;
        }

    }

    return closestWaypoint;

}

int NextWaypoint(double x, double y, double theta, vector<double> maps_x, vector<double> maps_y)
{

    int closestWaypoint = ClosestWaypoint(x,y,maps_x,maps_y);

    double map_x = maps_x[closestWaypoint];
    double map_y = maps_y[closestWaypoint];

    double heading = atan2( (map_y-y),(map_x-x) );

    double angle = abs(theta-heading);

    if(angle > pi()/2)
    {
        closestWaypoint++;
    }

    return closestWaypoint;

}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
vector<double> getFrenet(double x, double y, double theta, vector<double> maps_x, vector<double> maps_y)
{
    int next_wp = NextWaypoint(x,y, theta, maps_x,maps_y);

    int prev_wp;
    prev_wp = next_wp-1;
    if(next_wp == 0)
    {
        prev_wp  = maps_x.size()-1;
    }

    double n_x = maps_x[next_wp]-maps_x[prev_wp];
    double n_y = maps_y[next_wp]-maps_y[prev_wp];
    double x_x = x - maps_x[prev_wp];
    double x_y = y - maps_y[prev_wp];

    // find the projection of x onto n
    double proj_norm = (x_x*n_x+x_y*n_y)/(n_x*n_x+n_y*n_y);
    double proj_x = proj_norm*n_x;
    double proj_y = proj_norm*n_y;

    double frenet_d = distance(x_x,x_y,proj_x,proj_y);

    //see if d value is positive or negative by comparing it to a center point

    double center_x = 1000-maps_x[prev_wp];
    double center_y = 2000-maps_y[prev_wp];
    double centerToPos = distance(center_x,center_y,x_x,x_y);
    double centerToRef = distance(center_x,center_y,proj_x,proj_y);

    if(centerToPos <= centerToRef)
    {
        frenet_d *= -1;
    }

    // calculate s value
    double frenet_s = 0;
    for(int i = 0; i < prev_wp; i++)
    {
        frenet_s += distance(maps_x[i],maps_y[i],maps_x[i+1],maps_y[i+1]);
    }

    frenet_s += distance(0,0,proj_x,proj_y);

    return {frenet_s,frenet_d};

}

// Transform from Frenet s,d coordinates to Cartesian x,y
vector<double> getXY(double s, double d, vector<double> maps_s, vector<double> maps_x, vector<double> maps_y)
{
    int prev_wp = -1;

    while(s > maps_s[prev_wp+1] && (prev_wp < (int)(maps_s.size()-1) ))
    {
        prev_wp++;
    }

    int wp2 = (prev_wp+1)%maps_x.size();

    double heading = atan2((maps_y[wp2]-maps_y[prev_wp]),(maps_x[wp2]-maps_x[prev_wp]));
    // the x,y,s along the segment
    double seg_s = (s-maps_s[prev_wp]);

    double seg_x = maps_x[prev_wp]+seg_s*cos(heading);
    double seg_y = maps_y[prev_wp]+seg_s*sin(heading);

    double perp_heading = heading-pi()/2;

    double x = seg_x + d*cos(perp_heading);
    double y = seg_y + d*sin(perp_heading);

    return {x,y};

}

//Calculate Jerk Minimizing Trajectory
vector<double> JMT(vector< double> start, vector <double> end, double T)
{

    MatrixXd A = MatrixXd(3,3);
    A << pow(T, 3), pow(T, 4), pow(T,5),
        3 * pow(T,2), 4 * pow(T,3), 5 * pow(T,4),
        6 * T, 12 * pow(T,2), 20 * pow(T,3);

    MatrixXd B = MatrixXd(3,1);
    B << end[0]-(start[0]+start[1]*T+.5*start[2]*T*T),
                end[1]-(start[1]+start[2]*T),
                end[2]-start[2];

    MatrixXd Ai = A.inverse();

    MatrixXd C = Ai * B;

    vector <double> result = {start[0], start[1], .5*start[2]};
    for(int i = 0; i < C.size(); i++)
    {
        result.push_back(C.data()[i]);
    }

    return result;

}

double solve_jmt(vector<double> jmt, double delta){
  return jmt[0] + jmt[1] * delta + jmt[2]* pow(delta,2) + jmt[3]*pow(delta,3) + jmt[4]*pow(delta,4) + jmt[5]*pow(delta,5);
}

// Function to calculate the speed of other cars on the road
double calculate_speed(double vx, double vy){
  return (sqrt(pow(vx,2) + pow(vy,2)))/0.447;
}

// Function to predict the future S coordinate of another car
double calculate_future_s(double s, double speed, int prev_size){
  return s + (double(prev_size)*0.02*speed*0.447);
}

// Function to predict the future D coordinate of another car
double calculate_future_d(double s, double d, double vx, double vy, int prev_size, vector<double> maps_s, vector<double> maps_x, vector<double> maps_y){
  vector<double> xy = getXY(s, d, maps_s, maps_x, maps_y);
  double future_x = xy[0] + (double(prev_size) * vx * 0.02);
  double future_y = xy[1] + (double(prev_size) * vy * 0.02);
  double theta = atan2(vy,vx);
  return getFrenet(future_x, future_y, theta, maps_x, maps_y)[1];
}

// Function to determine which lane a vehicle is in
int get_lane(double d){
  if (0 < d and d < 4){
    return 0;
  }
  else if (4 < d and d < 8){
    return 1;
  }
  return 2;
}

string string_lane(int car_lane){
  if (car_lane == 0){
    return "Left";
  } else if (car_lane == 1){
    return "Center";
  }
  return "Right";
}

// Function to check whether or not a car is in the center of its lane
bool in_center_of_lane(double d){
  if ((1 < d and d < 3) or (5 < d and d < 7) or (9 < d and d < 11)){
    return true;
  }
  return false;
}

// Determine the best lane based on the speed and location of other cars on the road
int get_best_lane(vector<double> scores){
  return std::distance(scores.begin(), std::max_element(scores.begin(), scores.end()));
}

struct lane {
  double front_speed;
  double rear_speed;
  double front_buffer;
  double rear_buffer;
  double score;
};

struct my_car {
  double target_speed;
  double current_speed;
  string state = "KL";
  int preferred_lane = 1;
  double ref_vel = 0.0;
};

struct other_car{
  double buffer;
  double speed;
};

my_car self;
other_car front;
other_car rear;
other_car left_front;
other_car right_front;
other_car left_rear;
other_car right_rear;

lane Left;
lane Center;
lane Right;

int main() {
  uWS::Hub h;

  // Load up map values for waypoint's x,y,s and d normalized normal vectors
  vector<double> map_waypoints_x;
  vector<double> map_waypoints_y;
  vector<double> map_waypoints_s;
  vector<double> map_waypoints_dx;
  vector<double> map_waypoints_dy;

  // Waypoint map to read from
  string map_file_ = "../data/highway_map_bosch1_final.csv";
  // The max s value before wrapping around the track back to 0
  double max_s = 6945.554;

  ifstream in_map_(map_file_.c_str(), ifstream::in);

  string line;
  while (getline(in_map_, line)) {
    istringstream iss(line);
    double x;
    double y;
    float s;
    float d_x;
    float d_y;
    iss >> x;
    iss >> y;
    iss >> s;
    iss >> d_x;
    iss >> d_y;
    map_waypoints_x.push_back(x);
    map_waypoints_y.push_back(y);
    map_waypoints_s.push_back(s);
    map_waypoints_dx.push_back(d_x);
    map_waypoints_dy.push_back(d_y);
  }

  h.onMessage([&map_waypoints_x,&map_waypoints_y,&map_waypoints_s,&map_waypoints_dx,&map_waypoints_dy](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                     uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    //auto sdata = string(data).substr(0, length);
    //cout << sdata << endl;
    if (length && length > 2 && data[0] == '4' && data[1] == '2') {

      auto s = hasData(data);

      if (s != "") {
        auto j = json::parse(s);

        string event = j[0].get<string>();

        if (event == "telemetry") {
          // j[1] is the data JSON object

            // Main car's localization Data
            double car_x = j[1]["x"];
            double car_y = j[1]["y"];
            double car_s = j[1]["s"];
            double car_d = j[1]["d"];
            double car_yaw = j[1]["yaw"];
            double car_speed = j[1]["speed"];
            double speed_limit = 49.95;
            double target_speed = speed_limit;
            double front_speed;
            double critical_buffer = 13;
            double preferred_buffer = 18;
            double initial_buffer = 28;
            front.buffer = 1000;
            front.speed = speed_limit;
            rear.buffer = 1000;
            rear.speed = speed_limit;
            left_front.buffer = 1000;
            left_front.speed = speed_limit;
            left_rear.buffer = 1000;
            right_front.buffer = 1000;
            right_front.speed = speed_limit;
            right_rear.buffer = 1000;
            Right.front_buffer = 1000;
            Right.front_speed = speed_limit;
            Right.rear_buffer = 1000;
            Center.front_buffer = 1000;
            Center.front_speed = speed_limit;
            Center.rear_buffer = 1000;
            Left.front_buffer = 1000;
            Left.front_speed = speed_limit;
            Left.rear_buffer = 1000;
            bool front_clear = false;
            bool left_clear = false;
            bool right_clear = false;
            int target_lane = 1;

            // Previous path data given to the Planner
            auto previous_path_x = j[1]["previous_path_x"];
            auto previous_path_y = j[1]["previous_path_y"];
            // Previous path's end s and d values
            double end_path_s = j[1]["end_path_s"];
            double end_path_d = j[1]["end_path_d"];
            int prev_size = previous_path_x.size();

            // Sensor Fusion Data, a list of all other cars on the same side of the road.
            auto sensor_fusion = j[1]["sensor_fusion"];

            int car_lane = get_lane(end_path_d);
            bool in_center = in_center_of_lane(end_path_d);

            // Process the sensor fusion data to determine the speed and position of the closest cars in each lane
            for (int i = 0; i < sensor_fusion.size(); i++){
              double sensed_s = sensor_fusion[i][5];
              double sensed_d = sensor_fusion[i][6];
              double speed = calculate_speed(sensor_fusion[i][3], sensor_fusion[i][4]);
              double future_s = calculate_future_s(sensed_s, speed, prev_size);
              double future_d = calculate_future_d(sensed_s, sensed_d, sensor_fusion[i][3], sensor_fusion[i][4], prev_size, map_waypoints_s, map_waypoints_x, map_waypoints_y);
              int sensed_lane = get_lane(future_d);
              if (sensed_lane == 0){
                if (future_s > end_path_s){
                  if (future_s - end_path_s < Left.front_buffer){
                    Left.front_buffer = future_s - end_path_s;
                    Left.front_speed = speed;
                  }
                } else {
                  if (fabs(future_s - end_path_s) < Left.rear_buffer){
                    Left.rear_buffer = fabs(future_s - end_path_s);
                    Left.rear_speed = speed;
                  }
                }
              }
              if (sensed_lane == 1) {
                if (future_s > end_path_s){
                  if (future_s - end_path_s < Center.front_buffer){
                    Center.front_buffer = future_s - end_path_s;
                    Center.front_speed = speed;
                  }
                } else {
                  if (fabs(future_s - end_path_s) < Center.rear_buffer){
                    Center.rear_buffer = fabs(future_s - end_path_s);
                    Center.rear_speed = speed;
                  }
                }
              }
              if (sensed_lane == 2) {
                if (future_s > end_path_s){
                  if (future_s - end_path_s < Right.front_buffer){
                    Right.front_buffer = future_s - end_path_s;
                    Right.front_speed = speed;
                  }
                } else {
                  if (fabs(future_s - end_path_s) < Right.rear_buffer){
                    Right.rear_buffer = fabs(future_s - end_path_s);
                    Right.rear_speed = speed;
                  }
              }
            }
            }

            //Evaluate the sensor fusion data in relation to current lane
            if (car_lane == 0){
              front.buffer = Left.front_buffer;
              front.speed = Left.front_speed;
              rear.buffer = Left.rear_buffer;
              rear.speed = Left.rear_speed;
              right_front.buffer = Center.front_buffer;
              right_front.speed = Center.front_speed;
              right_rear.buffer = Center.rear_buffer;
              right_rear.speed = Center.rear_speed;
            }
            if (car_lane == 1){
              front.buffer = Center.front_buffer;
              front.speed = Center.front_speed;
              rear.buffer = Center.rear_buffer;
              rear.speed = Center.rear_speed;
              right_front.buffer = Right.front_buffer;
              right_front.speed = Right.rear_speed;
              right_rear.buffer = Right.rear_buffer;
              right_rear.speed = Right.rear_speed;
              left_front.buffer = Left.front_buffer;
              left_front.speed = Left.rear_speed;
              left_rear.buffer = Left.rear_buffer;
              left_rear.speed = Left.rear_speed;
            }
            if (car_lane == 2){
              front.buffer = Right.front_buffer;
              front.speed = Right.front_speed;
              rear.buffer = Right.rear_buffer;
              rear.speed = Right.rear_speed;
              left_front.buffer = Center.front_buffer;
              left_front.speed = Center.rear_speed;
              left_rear.buffer = Center.rear_buffer;
              left_rear.speed = Center.rear_speed;
            }

          critical_buffer *= car_speed/front.speed;
          preferred_buffer *= car_speed/front.speed;
          initial_buffer *= car_speed/front.speed;

          if (front.buffer > preferred_buffer){
            front_clear = true;
          }
          if ( car_lane - 1 >= 0 and left_front.buffer > preferred_buffer){
            left_clear = true;
          }
          if ( car_lane + 1 <= 2 and right_front.buffer > preferred_buffer){
            right_clear = true;
          }

          Left.score = Left.front_buffer * pow(Left.front_speed, 3);
          Center.score = Center.front_buffer * pow(Center.front_speed, 3) * 1.2;
          Right.score = Right.front_buffer * pow(Right.front_speed, 3);
          int best_lane = get_best_lane({Left.score, Center.score, Right.score});
          cout << "Best Lane: " << string_lane(best_lane) << endl;
          if (car_lane == 0){
            if (best_lane > 0){
              if (Center.front_buffer > preferred_buffer and Center.rear_buffer > 7 + pow(Center.rear_speed,5)/pow(car_speed,5)
               and self.state == "kl" and Left.front_buffer < initial_buffer * 1.5){
                self.preferred_lane += 1;
                self.state = "lcr";
              }
            }
          }
          if (car_lane == 2){
            if (best_lane < 2){
              if (Center.front_buffer > preferred_buffer and Center.rear_buffer > 7 + pow(Center.rear_speed,5)/pow(car_speed,5)
               and self.state == "kl" and Right.front_buffer < initial_buffer * 1.5) {
                self.preferred_lane -= 1;
                self.state = "lcl";
              }
            }
          }
          if (car_lane == 1){
            if (best_lane == 0){
              if (Left.front_buffer > preferred_buffer and Left.rear_buffer > 7 + pow(Center.rear_speed,5)/pow(car_speed,5) and self.state == "kl" and Center.front_buffer < initial_buffer){
                self.preferred_lane -= 1;
                self.state = "lcl";
              }
              else if (Right.score > Center.score and Right.front_buffer > preferred_buffer and Right.rear_buffer > 7 + pow(Center.rear_speed,5)/pow(car_speed,5) and self.state == "kl" and Center.front_buffer < initial_buffer){
                self.preferred_lane += 1;
                self.state = "lcr";
              }
            }
            if (best_lane == 2){
              if (Right.front_buffer > preferred_buffer and Right.rear_buffer > 7 + pow(Center.rear_speed,5)/pow(car_speed,5) and self.state == "kl" and Center.front_buffer < initial_buffer){
                self.preferred_lane += 1;
                self.state = "lcr";
              }
              else if (Left.score > Center.score and Left.front_buffer > preferred_buffer and Left.rear_buffer > 7 + pow(Center.rear_speed,5)/pow(car_speed,5) and self.state == "kl" and Center.front_buffer < initial_buffer){
                self.preferred_lane -= 1;
                self.state = "lcl";
              }
            }
          }

          if (self.preferred_lane < 0) {
            self.preferred_lane = 0;
          }
          if (self.preferred_lane > 2){
            self.preferred_lane = 2;
          }

          if (car_lane == self.preferred_lane and in_center == true){
              self.state = "kl";
          }

          if (front.buffer < critical_buffer){
            target_speed = front.speed - 2;
          } else if (front.buffer < preferred_buffer and self.state == "kl"){
            target_speed = front.speed;
          }
          if (car_lane != best_lane and car_speed < 40 and self.state == "kl" and front.buffer < 25){
            if (car_lane == 0 and abs(Left.front_buffer - Center.front_buffer) - 22) {
              target_speed = Left.front_speed -2;
            }
            if (car_lane == 2 and abs(Right.front_buffer - Center.front_buffer) - 22) {
              target_speed = Right.front_speed -2;
            }
          }

          cout << "===========================" << endl;
          cout << "Left Lane: " << endl;
          cout << "---------------------------" << endl;
          cout << "Front Buffer: " << Left.front_buffer << endl;
          cout << "Front Speed: " << Left.front_speed << endl;
          cout << "Rear Buffer: " << Left.rear_buffer << endl;
          cout << "Rear Speed: " << Left.rear_speed << endl;
          cout << "===========================" << endl;
          cout << "Center Lane: " << endl;
          cout << "---------------------------" << endl;
          cout << "Front Buffer: " << Center.front_buffer << endl;
          cout << "Front Speed: " << Center.front_speed << endl;
          cout << "Rear Buffer: " << Center.rear_buffer << endl;
          cout << "Rear Speed: " << Center.rear_speed << endl;
          cout << "===========================" << endl;
          cout << "Right Lane: " << endl;
          cout << "---------------------------" << endl;
          cout << "Front Buffer: " << Right.front_buffer << endl;
          cout << "Front Speed: " << Right.front_speed << endl;
          cout << "Rear Buffer: " << Right.rear_buffer << endl;
          cout << "Rear Speed: " << Right.rear_speed << endl;
          cout << "===========================" << endl;
          cout << "My Car: " << endl;
          cout << "---------------------------" << endl;
          cout << "Current Lane: " << string_lane(car_lane) << endl;
          cout << "Target Lane: " << string_lane(self.preferred_lane) << endl;
          if (self.state == "kl"){
              cout << "Current State = Keep  Lane " << endl;
            } else if (self.state == "lcl") {
              cout << "Current State = Lane Change Left " << endl;
            } else if (self.state == "lcr") {
              cout << "Current State = Lane Change Right " << endl;
            } else {
              cout << "Current State = Slow Down! " << endl;
            }
          cout << "Current Speed: " << car_speed << endl;
          cout << "Target Speed: " << target_speed << endl;

          vector<double> ptsx, ptsy;
          double ref_x = car_x;
          double ref_y = car_y;
          double ref_yaw = deg2rad(car_yaw);

          if (prev_size < 2){
            double prev_car_x = car_x - cos(car_yaw);
            double prev_car_y = car_y - sin(car_yaw);

            ptsx.push_back(prev_car_x);
            ptsx.push_back(car_x);
            ptsy.push_back(prev_car_y);
            ptsy.push_back(car_y);
          }
          else {
            ref_x = previous_path_x[prev_size-1];
            ref_y = previous_path_y[prev_size-1];

            double ref_x_prev = previous_path_x[prev_size-2];
            double ref_y_prev = previous_path_y[prev_size-2];
            ref_yaw = atan2(ref_y - ref_y_prev, ref_x - ref_x_prev);

            ptsx.push_back(ref_x_prev);
            ptsx.push_back(ref_x);
            ptsy.push_back(ref_y_prev);
            ptsy.push_back(ref_y);
          }
          double car_acc = 0;
          vector<double> start_d = {car_d, 0, car_acc};
          vector<double> end_d = {(2 + (double (self.preferred_lane) * 4)), 0, car_acc};
          vector<double> jmt_d = JMT(start_d, end_d, 1.8);

          vector<double> next_wp0 = getXY(car_s+45, solve_jmt(jmt_d, 1), map_waypoints_s, map_waypoints_x, map_waypoints_y);
          vector<double> next_wp1 = getXY(car_s+90, solve_jmt(jmt_d, 2), map_waypoints_s, map_waypoints_x, map_waypoints_y);
          vector<double> next_wp2 = getXY(car_s+135, solve_jmt(jmt_d, 3), map_waypoints_s, map_waypoints_x, map_waypoints_y);

          ptsx.push_back(next_wp0[0]);
          ptsx.push_back(next_wp1[0]);
          ptsx.push_back(next_wp2[0]);

          ptsy.push_back(next_wp0[1]);
          ptsy.push_back(next_wp1[1]);
          ptsy.push_back(next_wp2[1]);

          for (int i = 0; i < ptsx.size(); i++){
            double shift_x = ptsx[i] - ref_x;
            double shift_y = ptsy[i] - ref_y;

            ptsx[i] = (shift_x * cos(0-ref_yaw) - shift_y * sin(0-ref_yaw));
            ptsy[i] = (shift_x * sin(0-ref_yaw) + shift_y * cos(0-ref_yaw));
          }

          tk::spline s;

          s.set_points(ptsx, ptsy);

          vector<double> next_x_vals;
          vector<double> next_y_vals;

          for (int i = 0; i < previous_path_x.size(); i++){
            next_x_vals.push_back(previous_path_x[i]);
            next_y_vals.push_back(previous_path_y[i]);
          }

          double target_x = 30.0;
          double target_y = s(target_x);
          double target_dist = sqrt(pow(target_x,2) + pow(target_y, 2));

          double x_add_on = 0;

          for (int i = 1; i <= 50 - previous_path_x.size(); i++){
            if (self.ref_vel > target_speed + 0.2){
              self.ref_vel -= 0.27;
            } else if (self.ref_vel < target_speed){
              if (self.ref_vel < 44 and front.buffer > initial_buffer and self.state == "kl"){
                self.ref_vel += .444;
              } else {
                self.ref_vel += min(.36, speed_limit - self.ref_vel);
              }
            }
            if (self.ref_vel > speed_limit){
              self.ref_vel = speed_limit - 0.015;
            }
            double N = (target_dist/(0.02*self.ref_vel/2.24));
            double x_point = x_add_on+target_x/N;
            double y_point = s(x_point);

            x_add_on = x_point;

            double x_ref = x_point;
            double y_ref = y_point;

            x_point = (x_ref * cos(ref_yaw) - y_ref * sin(ref_yaw));
            y_point = (x_ref * sin(ref_yaw) + y_ref * cos(ref_yaw));

            x_point += ref_x;
            y_point += ref_y;

            next_x_vals.push_back(x_point);
            next_y_vals.push_back(y_point);
          }

          json msgJson;
          // TODO: define a path made up of (x,y) points that the car will visit sequentially every .02 seconds
          msgJson["next_x"] = next_x_vals;
          msgJson["next_y"] = next_y_vals;

          auto msg = "42[\"control\","+ msgJson.dump()+"]";

          //this_thread::sleep_for(chrono::milliseconds(1000));
          ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);

        }
      } else {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }
  });

  // We don't need this since we're not using HTTP but if it's removed the
  // program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                     size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
    } else {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                         char *message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  if (h.listen(port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}
