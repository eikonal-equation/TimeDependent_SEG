/** ------ Libraries ---------------------------------------------------------*/
#include <chrono>
#include <iostream>
#include <vector>
#include <cmath>

/** ------ Project-specific header files -------------------------------------*/
#include "GlobalConfiguration.hpp"
#include "CAdversarialPlan.hpp"
#include "CTimeDependentPlanner.hpp"
#include "CMovingObserver.hpp"
#include "CStationaryTerrain.hpp"
#include "SimpleFunctions.hpp"
#include "ObserverPaths.hpp"

/** ------ Namespaces --------------------------------------------------------*/
using namespace std;

/*==============================================================================
  Example 1:
    - No obstacles
    - Two observer trajectories
    - Omnidirectional sensors
==============================================================================*/
CAdversarialPlan noObstacleTest() {
  /* Define observer paths */
  const vector<Function1D> observer_x_paths = {&xsimple1, &xsimple2};
  const vector<Function1D> observer_y_paths = {&ysimple1, &ysimple2};

  /* Create source and target */
  const vector<double> source = {0.1,0.1};
  const vector<double> target = {0.9, 0.9};

  /* Create terrain */
  const vector<vector<double>> obstacles = {};
  shared_ptr<CStationaryTerrain> terrain = make_shared<CStationaryTerrain>(obstacles);

  /* Set basic configurations for the grid */
  constexpr double phys_min = 0, phys_max = 1, max_time = 4;
  constexpr int n = 201;

  /* Calculate number of timesteps based on CFL condition: */
  double max_dt_allowed;
  if (NINE_POINT_STENCIL) {
    max_dt_allowed = (phys_max - phys_min)/(n-1);
  } else {
    max_dt_allowed = (phys_max - phys_min)/((n-1)*(SQRT2));
  }
  const int nt = floor(max_time / max_dt_allowed) + 1;

  /* Time factor = ratio of path tracer timesteps to HJB solver timesteps */
  constexpr int time_factor = 10;

  /* Create observers */
  const int n_observers = observer_x_paths.size();
  vector<shared_ptr<CObserver>> observer_pointers = {};
  for(int k = 0; k < n_observers; ++k) {
    observer_pointers.push_back(make_shared<CMovingObserver>(
        observer_x_paths[k], observer_y_paths[k], terrain, &inverseSquared, 
        n, nt, phys_min, phys_max, max_time));
  }

  /* Define fixed Lambda planner */
  shared_ptr<CFixedLambdaPlanner> fixedplan = make_shared<CTimeDependentPlanner>(
      observer_pointers, source, target, terrain, &inverseSquared, &constant2, 
      n, nt, n_observers, time_factor, phys_min, phys_max, max_time);
  
  return CAdversarialPlan(fixedplan);
}

/*==============================================================================
  Example 2:
    - One central obstacle
    - Two observer trajectories
    - Omnidirectional sensors
==============================================================================*/
CAdversarialPlan obstacleTest() {
  /* Define observer paths */
  const vector<Function1D> observer_x_paths = {&xsimple1, &xsimple2};
  const vector<Function1D> observer_y_paths = {&ysimple1, &ysimple2};

  /* Create source and target */
  const vector<double> source = {0.1,0.1};
  const vector<double> target = {0.9, 0.9};

  /* Create terrain */
  const vector<vector<double>> obstacles = {{0.5,0.7,0.4,0.6}};
  shared_ptr<CStationaryTerrain> terrain = make_shared<CStationaryTerrain>(obstacles);

  /* Set basic configurations for the grid */
  constexpr double phys_min = 0, phys_max = 1, max_time = 4;
  constexpr int n = 201;

  /* Calculate number of timesteps based on CFL condition: */
  double max_dt_allowed;
  if (NINE_POINT_STENCIL) {
    max_dt_allowed = (phys_max - phys_min)/(n-1);
  } else {
    max_dt_allowed = (phys_max - phys_min)/((n-1)*(SQRT2));
  }
  const int nt = floor(max_time / max_dt_allowed) + 1;

  /* Time factor = ratio of path tracer timesteps to HJB solver timesteps */
  constexpr int time_factor = 1;

  /* Create observers */
  const int n_observers = observer_x_paths.size();
  vector<shared_ptr<CObserver>> observer_pointers = {};
  for(int k = 0; k < n_observers; ++k) {
    observer_pointers.push_back(make_shared<CMovingObserver>(
        observer_x_paths[k], observer_y_paths[k], terrain, &inverseSquared, 
        n, nt, phys_min, phys_max, max_time));
  }

  /* Define fixed Lambda planner */
  shared_ptr<CFixedLambdaPlanner> fixedplan = make_shared<CTimeDependentPlanner>(
      observer_pointers, source, target, terrain, &inverseSquared, &constant2, 
      n, nt, n_observers, time_factor, phys_min, phys_max, max_time);
  
  return CAdversarialPlan(fixedplan);
}

/*==============================================================================
  Example 3:
    - Three rectangular obstacles
    - Two observer trajectories
    - Angle-restricted sensors
==============================================================================*/
CAdversarialPlan figureEightAnisoTest() {
  /* Define observer paths */
  const vector<Function1D> observer_x_paths = {&xtilted1, &xtilted2};
  const vector<Function1D> observer_y_paths = {&ytilted1, &ytilted2};
  const vector<Function1D> observer_x_derivs = {&xtiltedder1, &xtiltedder2};
  const vector<Function1D> observer_y_derivs = {&ytiltedder1, &ytiltedder2};

  /* Create source and target */
  const vector<double> source = {0.2,0.2};
  const vector<double> target = {0.9, 0.9};

  /* Create terrain */
  const vector<vector<double>> obstacles = {{0.1,0.3,0.7,0.9}, {0.4,0.6,0.4,0.6}, {0.7,0.9,0.1,0.3}};
  shared_ptr<CStationaryTerrain> terrain = make_shared<CStationaryTerrain>(obstacles);

  /* Set basic configurations for the grid */
  constexpr double phys_min = 0, phys_max = 1, max_time = 4;
  constexpr int n = 201;

  /* Calculate number of timesteps based on CFL condition: */
  double max_dt_allowed;
  if (NINE_POINT_STENCIL) {
    max_dt_allowed = (phys_max - phys_min)/(n-1);
  } else {
    max_dt_allowed = (phys_max - phys_min)/((n-1)*(SQRT2));
  }
  const int nt = floor(max_time / max_dt_allowed) + 1;

  /* Time factor = ratio of path tracer timesteps to HJB solver timesteps */
  constexpr int time_factor = 1;

  /* Create observers */
  const int n_observers = observer_x_paths.size();
  vector<shared_ptr<CObserver>> observer_pointers = {};
  for(int k = 0; k < n_observers; ++k) {
    observer_pointers.push_back(make_shared<CMovingObserver>(
        observer_x_paths[k], observer_y_paths[k], terrain, &inverseSquared, 
        n, nt, phys_min, phys_max, max_time, 
        observer_x_derivs[k], observer_y_derivs[k], 2*PI/3));
  }

  /* Define fixed Lambda planner */
  shared_ptr<CFixedLambdaPlanner> fixedplan = make_shared<CTimeDependentPlanner>(
      observer_pointers, source, target, terrain, &inverseSquared, &constant2, 
      n, nt, n_observers, time_factor, phys_min, phys_max, max_time);
  
  return CAdversarialPlan(fixedplan);
}

/*==============================================================================
  Example 4:
    - Maze-like domain
    - Four observer trajectories
    - Angle-restricted sensors
==============================================================================*/
CAdversarialPlan mazeTest() {
  /* Define observer paths */
  const vector<Function1D> observer_x_paths = {&xmazecircle1, &xmazecircle2, &xmazecircle3, &xmazecircle4};
  const vector<Function1D> observer_y_paths = {&ymazecircle1, &ymazecircle2, &ymazecircle3, &ymazecircle4};
  const vector<Function1D> observer_x_derivs = {&xmazecircleder1, &xmazecircleder2, &xmazecircleder3, &xmazecircleder4};
  const vector<Function1D> observer_y_derivs = {&ymazecircleder1, &ymazecircleder2, &ymazecircleder3, &ymazecircleder4};

  /* Create source and target */
  const vector<double> source = {0.95,0.15};
  const vector<double> target = {0.375, 0.5};

  /* Create terrain */
  const vector<vector<double>> obstacles = {{0.45,0.8,0.1,0.2}, {0.15,0.8,0.8,0.9}, 
                                            {0.7,0.8,0.2,0.75}, {0.3,0.55,0.35,0.45}, 
                                            {0.45,0.55,0.45,0.6}};
  shared_ptr<CStationaryTerrain> terrain = make_shared<CStationaryTerrain>(obstacles);

  /* Set basic configurations for the grid */
  constexpr double phys_min = 0, phys_max = 1, max_time = 4;
  constexpr int n = 201;

  /* Calculate number of timesteps based on CFL condition: */
  double max_dt_allowed;
  if (NINE_POINT_STENCIL) {
    max_dt_allowed = (phys_max - phys_min)/(n-1);
  } else {
    max_dt_allowed = (phys_max - phys_min)/((n-1)*(SQRT2));
  }
  const int nt = floor(max_time / max_dt_allowed) + 1;

  /* Time factor = ratio of path tracer timesteps to HJB solver timesteps */
  constexpr int time_factor = 1;

  /* Create observers */
  const int n_observers = observer_x_paths.size();
  vector<shared_ptr<CObserver>> observer_pointers = {};
  for(int k = 0; k < n_observers; ++k) {
    observer_pointers.push_back(make_shared<CMovingObserver>(
        observer_x_paths[k], observer_y_paths[k], terrain, &inverseSquared, 
        n, nt, phys_min, phys_max, max_time, 
        observer_x_derivs[k], observer_y_derivs[k], 2*PI/3));
  }

  /* Define fixed Lambda planner */
  shared_ptr<CFixedLambdaPlanner> fixedplan = make_shared<CTimeDependentPlanner>(
      observer_pointers, source, target, terrain, &inverseSquared, &constant2, 
      n, nt, n_observers, time_factor, phys_min, phys_max, max_time);
  
  return CAdversarialPlan(fixedplan);
}

int main(int argc, char* argv[]) {
  /* Start timer */
  auto t1 = chrono::high_resolution_clock::now();

  /* Select example based on command line argument */
  CAdversarialPlan adversarial_plan;
  string filename;
  if (argc < 2) {
    cout << "No argument, terminating..." << endl;
    return 0;
  } else if (string(argv[1]) == "1") {
    cout << "Running Example 1: No Obstacles" << endl;
    adversarial_plan = noObstacleTest();
    filename = "noobstacle";
  } else if (std::string(argv[1]) == "2") {
    std::cout << "Running Example 2: Central Obstacle" << std::endl;
    adversarial_plan = obstacleTest();
    filename = "obstacle";
  } else if (std::string(argv[1]) == "3") {
    std::cout << "Running Example 3: Angle-restricted observers" << std::endl;
    adversarial_plan = figureEightAnisoTest();
    filename = "figure8aniso";
  } else if (std::string(argv[1]) == "4") {
    std::cout << "Running Example 4: Maze-like domain" << std::endl;
    adversarial_plan = mazeTest();
    filename = "maze";
  } else {
    cout << "Invalid argument, terminating..." << endl;
    return 0;
  }

  /* Find Pareto Front (for Examples 1 and 2) */
  if (filename == "noobstacle" || filename == "obstacle") {
    adversarial_plan.getParetoFront(filename, 100);
  }

  /* Find Nash Equilibrium */
  if (filename == "figure8aniso") {
    /* This example is symmetric, so no need to optimize over lambda */
    adversarial_plan.findNash(filename, 0);
  } else {
    adversarial_plan.findNash(filename, 200);
  }

  /* Stop timer */
  auto t2 = std::chrono::high_resolution_clock::now();
  cout << "Test complete, took " << chrono::duration_cast<chrono::milliseconds>(t2-t1).count() << " milliseconds" << endl;
  return 0;
}