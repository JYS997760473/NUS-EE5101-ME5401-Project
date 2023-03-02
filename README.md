# NUS-EE5101/ME5401-Project

This repository is for EE5101/ME5401 Linear System Project in NUS.

## Abstract

Self-balancing two-wheeled vehicle is quite useful for people who have difficulty in riding a bicycle and it has been developing by applying modern control theory. In this project, we construct a linear state-space system to simulate the true state of this vehicle in real world, and use pole-placement control, LQR control , controller with observer, decoupling control and integral control separately to control and simulate this system using MATLAB and Simulink. Controller and observer’s performance on different inputs and disturbances are also included in this paper and final results prove that these control theories are suitable for this system.

## Requirments

The requirements are all in `project-2022.pdf`.

## Introduction

`ME5401.m` is the main function contain six main requirements in `project-2022.pdf`.

`place_pole.m` contains a pole placement method by using cannocial form.

`my_lqr.m` contains a pole placement method by LQR to find the positive definite solution of the Riccati equation.

`unity_place.m` contains a pole placement method by using Ackermann's formula unity rank method.

`model1-6.slx` are the Simulink models for six different requirements.

## Usage

### With MATLAB Desktop

Type the `run` button to run `ME5401.m` file in MATLAB Destkop.

### Without MATLAB Desktop (In Linux or MacOS terminal)

Run
```bash
./ME5401
```

## Demo and Results

Final results can be found in  `ME5401project.pdf`.