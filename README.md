# F-16 Angle of Attack Estimation using 6-DoF Simulation and Extended Kalman Filter

**Avionic Systems Project** — MATLAB implementation of a nonlinear 6-DoF flight
dynamics simulation and Extended Kalman Filter (EKF) for estimating Angle of
Attack (AoA) on the F-16 Fighting Falcon without relying on a direct AoA vane measurement.

---

## Overview

The project consists of two parts:

1. **6-DoF Truth Model** — a full nonlinear simulation of F-16 flight dynamics
   using the NASA aerodynamic dataset and a 7-layer ISA atmosphere model.
2. **Extended Kalman Filter** — a 9-state EKF that fuses gyroscope, accelerometer,
   and airspeed measurements to estimate AoA in real time.

---

## File Structure

| File | Description |
|---|---|
| `F16_EOM.m` | 12-state 6-DoF equations of motion (force, moment, kinematic, navigation) |
| `F16_Atmosphere.m` | 7-layer ISA atmosphere model (temperature, pressure, density, speed of sound) |
| `F16_AeroTables.m` | Aerodynamic coefficient lookup tables with linear interpolation |
| `F16_Sensors.m` | Sensor models: accelerometers, gyroscopes, pitot tube, AoA vane, barometer |
| `F16_Trim.m` | Nonlinear trim solver for straight and level flight |
| `F16_EKF.m` | 9-state EKF — predict, update, and AoA extraction |
| `F16_EKF_Init.m` | EKF initialisation: state, P, Q, R |
| `Run_F16_Simulation.m` | Run the 6-DoF truth model and plot results |
| `Run_F16_EKF.m` | Run the full EKF simulation and compute RMS performance |
| `Test_AeroTables.m` | Sanity checks on aerodynamic coefficient tables |
| `Test_Atmosphere.m` | ISA model verification at key altitudes |
| `Test_Sensors.m` | Sensor noise and bias visualisation |

---

## How to Run

1. Open MATLAB and set the working directory to this folder.
2. Run the trim solver first to verify the flight condition:
   ```matlab
   run('F16_Trim.m')
   ```
3. Run the truth model simulation:
   ```matlab
   run('Run_F16_Simulation.m')
   ```
4. Run the EKF estimation:
   ```matlab
   run('Run_F16_EKF.m')
   ```

---

## Flight Condition

| Parameter | Value |
|---|---|
| Airspeed | 100 m/s |
| Altitude | 6000 m |
| Angle of attack | 8.2816° |
| Elevator | −0.3957° |
| Thrust | 6279.17 N |

---

## EKF Design Summary

- **States (9):** body velocities $(u, v, w)$, angular rates $(p, q, r)$, Euler angles $(\phi, \theta, \psi)$
- **Measurements (6):** gyroscopes $(p, q, r)$, accelerometers $(f_x, f_z)$, airspeed $(V)$
- **AoA estimate:** derived as $\hat{\alpha} = \arctan(w/u)$ — not a filter state
- **Jacobians:** computed numerically via finite difference
- **Covariance update:** Joseph form for numerical stability

---

## Requirements

- MATLAB R2019b or later
- Optimization Toolbox (for `fsolve` in `F16_Trim.m`)

---

## Reference

- Stevens, Lewis & Johnson, *Aircraft Control and Simulation*, 3rd ed., Wiley, 2015
- Nguyen et al., *NASA Technical Paper 1538*, 1979 (F-16 aerodynamic data)
