# MPC Bicycle Model Simulation

![MPC_Bicycle](https://github.com/rishabh-bavithran/MPC-Bicycle-Model/assets/145865695/08ea369a-752f-4846-b601-bc4958ae2137)

This project simulates an MPC (Model Predictive Controller) for a simplified bicycle model. The goal is to implement a control system that can effectively navigate a bicycle model based on predictive modeling and control techniques.

## Project Overview

The simulation assumes the following:
- Constant forward velocity and no-slip conditions.
- Small angle assumption for constant cornering stiffness and model linearization.
- Mathematical model derived for a turning bicycle.
- MPC controller applied with a horizon period of 5 seconds.

## Simulation Details

- **Total Simulation Time:** 10 seconds
- **MPC Calculation Interval:** 0.01 seconds
- **Horizion Period:** 5 seconds 
- **Cost Function:** Quadratic
- **Cost Function Parameters:** error along y-axis and change in steering angle
- **Predicted State Space Variables:** Obtained using Forward Euler Integration
- **Visualization:** Results are animated post-simulation
