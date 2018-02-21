# Udacity-Bosch Path Planning Challenge
Your challenge is to design a path planner that creates smooth, safe paths along a 3-lane highway with traffic. Successful path planners will be able to keep inside their lanes, avoid hitting other cars, and pass slower moving traffic - all by using localization, sensor fusion, and map data. Path planners will be ranked by the time it takes them to complete a lap around a test track.

#### The map of the highway is in data/highway_map_bosch1.csv
Each waypoint in the list contains  [x,y,s,dx,dy] values. x and y are the waypoint's map coordinate position, the s value is the distance along the road to get to that waypoint in meters, the dx and dy values define the unit normal vector pointing outward of the highway loop.

## Basic Build Instructions

1. Clone this repo.
2. Make a build directory: `mkdir build && cd build`
3. Compile: `cmake .. && make`
4. Run it: `./path_planning`.

Here is the data provided from the Simulator to the C++ Program

#### Main car's localization Data (No Noise)

["x"] The car's x position in map coordinates

["y"] The car's y position in map coordinates

["s"] The car's s position in frenet coordinates

["d"] The car's d position in frenet coordinates

["yaw"] The car's yaw angle in the map

["speed"] The car's speed in MPH

#### Previous path data given to the Planner

//Note: Return the previous list but with processed points removed, can be a nice tool to show how far along
the path has processed since last time. 

["previous_path_x"] The previous list of x points previously given to the simulator

["previous_path_y"] The previous list of y points previously given to the simulator

#### Previous path's end s and d values 

["end_path_s"] The previous list's last point's frenet s value

["end_path_d"] The previous list's last point's frenet d value

#### Sensor Fusion Data, a list of all other car's attributes on the same side of the road. (No Noise)

["sensor_fusion"] A 2d vector of cars and then that car's [car's unique ID, car's x position in map coordinates, car's y position in map coordinates, car's x velocity in m/s, car's y velocity in m/s, car's s position in frenet coordinates, car's d position in frenet coordinates. 

## Details

1. The car uses a perfect controller and will visit every (x,y) point it recieves in the list every .02 seconds. The units for the (x,y) points are in meters and the spacing of the points determines the speed of the car. The vector going from a point to the next point in the list dictates the angle of the car. Acceleration both in the tangential and normal directions is measured along with the jerk, the rate of change of total Acceleration. The (x,y) point paths that the planner recieves should not have a total acceleration that goes over 10 m/s^2, also the jerk should not go over 50 m/s^3. (NOTE: As this is BETA, these requirements might change. Also currently jerk is over a .02 second interval, it would probably be better to average total acceleration over 1 second and measure jerk from that.

2. There will be some latency between the simulator running and the path planner returning a path, with optimized code usually its not very long maybe just 1-3 time steps. During this delay the simulator will continue using points that it was last given, because of this its a good idea to store the last points you have used so you can have a smooth transition. previous_path_x, and previous_path_y can be helpful for this transition since they show the last points given to the simulator controller with the processed points already removed. You would either return a path that extends this previous path or make sure to create a new path that has a smooth transition with this last path.

## Tips

A really helpful resource for doing this project and creating smooth trajectories was using http://kluge.in-chemnitz.de/opensource/spline/, the spline function is in a single hearder file is really easy to use.

---

## Dependencies

* cmake >= 3.5
 * All OSes: [click here for installation instructions](https://cmake.org/install/)
* make >= 4.1
  * Linux: make is installed by default on most Linux distros
  * Mac: [install Xcode command line tools to get make](https://developer.apple.com/xcode/features/)
  * Windows: [Click here for installation instructions](http://gnuwin32.sourceforge.net/packages/make.htm)
* gcc/g++ >= 5.4
  * Linux: gcc / g++ is installed by default on most Linux distros
  * Mac: same deal as make - [install Xcode command line tools]((https://developer.apple.com/xcode/features/)
  * Windows: recommend using [MinGW](http://www.mingw.org/)
* [uWebSockets](https://github.com/uWebSockets/uWebSockets)
  * Run either `install-mac.sh` or `install-ubuntu.sh`.
  * If you install from source, checkout to commit `e94b6e1`, i.e.
    ```
    git clone https://github.com/uWebSockets/uWebSockets 
    cd uWebSockets
    git checkout e94b6e1
    ```

# Udacity Track Release
Go here to download the [simulator](https://github.com/udacity/Bosch-Challenge/releases/tag/v1.0) and original project track: https://github.com/udacity/CarND-Path-Planning-Project

# Bosch Track Release
The Bosch track is a highway that continuously goes down the same direction but with some slight curviness. 

The competitive challenge element added to this project is a brand new highway to test on, along with a timer element that measures how quickly the target 2.75 miles can be completed without incident. 

The guidelines for a valid submission are listed in the rubric below.

# Final Evaluation Track
Your submission will be run on a private evaluation track that will be similar to the Bosch track with the exception that the car placement won't be random.

# Project Rubric 


### Compilation
The code compiles correctly. Code must compile without errors with cmake and make.

### Valid Trajectories
**The car is able to drive at least 2.75 miles without incident**        

The top right screen of the simulator shows the current/best miles driven without incident. Incidents include exceeding acceleration/jerk/speed, collision, and driving outside of the lanes. Each incident case is also listed below in more detail.

----------------------------------------------------------------------------------------------------------------------------------------

**The car drives according to the speed limit**

The car doesn't drive faster than the speed limit of 50 miles per hour.

----------------------------------------------------------------------------------------------------------------------------------------

**Max acceleration and jerk are not exceeded**

The car does not exceed a total acceleration of 10 m/s^2 and a jerk of 10 m/s^3.

----------------------------------------------------------------------------------------------------------------------------------------

**Car does not have collisions**

The car must not come into contact with any of the other cars or objects on the road.

----------------------------------------------------------------------------------------------------------------------------------------

**The car stays in its lane except when changing lanes**

During lane changes, the car doesn't spend more than 3 seconds outside the lane lines. During normal operation the car stays inside the lane lines with the flow of traffic.

----------------------------------------------------------------------------------------------------------------------------------------

**The car completes its trip within 360 seconds**

Any submissions that don't finish within 360 seconds will be considered invalid.

----------------------------------------------------------------------------------------------------------------------------------------

**Submission file must meet naming standards**

Your file must be named main.cpp and include all the code required to successfully complete the Challenge. In addition, you must use the following map name in your code: highway_map_bosch1.csv

## Model Documentation:

The first step in my path planning model is to process the sensor fusion data to determine the speed and position of the closest cars in each lane (the closest car in front of my car and the closest car behind my car). Based on the current speed and trajectory of the other cars on the road, I am able to predict where they will be at n timesteps in the future and this helps me react to other cars beavior including changing lanes to and from my current lane.

After determining the speed and position of the nearby cars in each lane, I define a score for each lane based on the distance of the nearest car in front of my car, and that cars current speed. 
```
Left.score = Left.front_buffer * pow(Left.front_speed, 3);
Center.score = Center.front_buffer * pow(Center.front_speed, 3) * 1.2;
Right.score = Right.front_buffer * pow(Right.front_speed, 3);
```
As you can see, the speed is given more weight compared to the front buffer distance, and I have given a bias towards the center lane by adding 20% to the center lane score.

Any time my car is less than a specified distance away from the nearest car in front in it's current lane, then I apply a lane change criteria to determine whether or not I should change lanes, and if so, which lane I should change to. 

My lane changing criteria is like this:

* The lane with the highest score is considered to be the best lane.
* If my current lane is not the best lane, then I will attempt to change lanes towards the best lane.
* If it is not safe to change lanes based on the speed and position of cars in neighboring lanes, then I will maintain my current lane and slow down to the speed of the closest car in front of me.

To generate trajectories, I used a combination of spline fitting and Jerk Minimizing Trajectory (JMT). 
