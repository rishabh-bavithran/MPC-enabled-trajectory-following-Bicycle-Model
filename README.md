Copy code
# MPC Bicycle Model Simulation

![MPC_Bicycle_model](https://github.com/Kelvin4915/MPC_Bicycle_Model/assets/134540002/87e69b86-82bb-4c7e-b067-006fc6b6a1b6)

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
- **State Space Variables:** Interpolated with a timestep divided into 20 segments.
- **Visualization:** Results are animated post-simulation.
