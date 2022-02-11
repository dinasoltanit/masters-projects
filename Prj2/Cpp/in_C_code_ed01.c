/*
Language: C
Author: Dina Soltani Tehrani.
Course: Advanced CFD1 - by Dr. Taeibi Rahni
Sharif University of Technology

Equation: ut = -c * ux
Method: Finite difference, FTCS, Explicit
*/
# include <cstdlib>
# include <iostream>
# include <fstream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>
//***** Constants ******
const int itermax = 100; // max number of iterations
const float dx = 0.5; // mesh size (m)
const int L = 300; // Length (m)
//const float dx1 = 0.5, dx2=0.1, dx3=0.05;

void main()
{
float dt = 0.75/itermax, err = 0, pi = 4*atan(1.0);
// dt: time-step size
// err: related to the accuracy criterion
nx = L/dx + 1;
float u_old[nx], u_new[nx]