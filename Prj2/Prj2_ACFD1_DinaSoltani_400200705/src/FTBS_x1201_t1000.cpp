# include <cstdlib>
# include <iostream>
# include <fstream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>

using namespace std;

int main ( );
int i4_modp ( int i, int j );
int i4_wrap ( int ival, int ilo, int ihi );
double *initial_condition ( int nx, double x[] );
double *r8vec_linspace_new ( int n, double a, double b );
void timestamp ( );

//****************************************************************************

int main ( )

//****************************************************************************
//
// Language: C++
//  Author: Dina Soltani Tehrani.
//  Course: Advanced CFD1 - by Dr. Taeibi Rahni
//  Sharif University of Technology
//  Licence: GNU LGPL
//
//  Equation: ut = -c * ux
//  Method: Finite difference, FTBS Method, Explicit
//  Pointer-based Method
//
//    The finite difference form is as bellow:
//
//      U(X,T+dt) - U(X,T)          ( U(X-dx,T+dt) + U(X+dx,T+dt) )
//      ------------------  = -c * --------------------------------- + F(X, T+dt)
//               dt                              2 * dx
//
{
  double a;
  double b;
  double c;
  double L;
  double duration;
  double  pi = 4*atan(1.0);
  string data_filename = "FTBS_x1201_t1000_p.txt";
  ofstream data_unit;
  double dt; // time-step
  double dx; // size-step
  int i;
  int j;
  int jm1;
  int jp1;
  int nx;
  int nt;
  int nt_step;
  int plotstep;
  double t;
  double *u; // old velocity
  double *unew; // new velocity
  double *x;

  cout << "\n";
  cout << "Finite Difference: FTBS for 1D Wave Equation:\n";
  cout << "  C++ language\n";
  cout << "\n";
  cout << "  Solve the constant-velocity advection equation in 1D,\n";
  cout << "    du/dt = - c du/dx\n";
  cout << "  over the interval:\n";
  cout << "    0.0 <= x <= 300.0\n";
  cout << "  with a given boundary conditions, and\n";
  cout << "  with a given initial condition\n";
  cout << "    u(0,x) = 80 * sin( (pi/60)* (x-40) ) for 40.0 <= x <= 100.0\n";
  cout << "           = 0 elsewhere.\n";
  cout << "\n";
  cout << "  We use a method known as FTBS:\n";
  cout << "   FT: Forward Time  : du/dt = (u(t+dt,x)-u(t,x))/dt\n";
  cout << "   BS: Backward Space: du/dx = (u(t,x)-u(t,x-dx))/dx\n";
  cout << "\n";

  timestamp ( );

  L = 300; // Length (m)
  nx = 1201; // Mid Space Resolution
  dx = L / ( double ) ( nx - 1 ); // size-step
  a = 0.0; // x-start
  b = 300.0; // x-end
  x = r8vec_linspace_new ( nx, a, b );
  nt = 1000; // Time Resolution
  duration = 0.6; // Duration (s)
  dt = duration / ( double ) ( nt ); // time-step
  c = 330.0;

  u = initial_condition ( nx, x ); // Initializing old velocity with initial condition.

//  ***** Opening data file, and writing solutions as they are computed. *****

  data_unit.open ( data_filename.c_str ( ) );

  // Writing the initial data.
  t = 0.0;
  data_unit << x[0]
            << "  " << u[0] << "\n";

  // Writing data as we move forward in space.
  for ( j = 0; j < nx; j++ )
  {
    data_unit << x[j]
              << "  " << u[j] << "\n";
  }
  data_unit << "zone";
  data_unit << "\n";



//*************************************************************************************************
  nt_step = 100;
  unew = new double[nx]; // Creating the new velocity vector

// ***** Core Calculation: begin *****
  for ( i = 0; i < nt; i++ )
  {
    //data_unit << "zone\n";
    // First do this for all x in xn:
    for ( j = 0; j < nx; j++ )
    {
      jm1 = i4_wrap ( j - 1, 0, nx - 1 ); // phi_n_j-1
      //jp1 = i4_wrap ( j + 1, 0, nx - 1 ); // phi_n_j+1
      unew[j] = u[j] - (c * dt / dx) * ( u[j] - u[jm1] );
      //if ( u[j] < 1e-04 ) { u[j] = 0; }
    }
    // Then, do this, again for all x in xn:
    for ( j = 0; j < nx; j++ )
    {
      u[j] = unew[j];
      //data_unit << x[j] << "  " << u[j] << "\n";
    }

    if ( i == nt_step - 1 )
    {
      t = ( double ) ( i ) * dt;
      for ( j = 0; j < nx; j++ )
      {
        //if ( u[j] < 1e-05 ) { u[j] = 0; }
        data_unit << x[j]
                  << "  " << u[j] << "\n";
      }
      data_unit << "zone";
      data_unit << "\n";
      nt_step = nt_step + 100;
    }

    //data_unit << "zone";
    //data_unit << "\n";
  }
// ***** Core Calculation: end *****


//
//  ***** Closing the data file as the computation is done. *****
//
  data_unit.close ( );

  cout << "\n";
  cout << "  Plot data written to the file \"" << data_filename << "\"\n";
//

//  Free memory.
//
  delete [] u;
  delete [] unew;
  delete [] x;
//
//  Terminate.
//
  cout << "\n";
  cout << "  Finite Difference: FTBS for 1D Wave Equation\n";
  cout << "  Pointer-based Method\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}

//
// ***** Functions *****

//****************************************************************************
int i4_modp ( int i, int j )
{
  int value;

  if ( j == 0 )
  {
    cerr << "\n";
    cerr << "I4_MODP - Fatal error!\n";
    cerr << "  I4_MODP ( I, J ) called with J = " << j << "\n";
    exit ( 1 );
  }

  value = i % j;

  if ( value < 0 )
  {
    value = value + abs ( j );
  }

  return value;
}

//****************************************************************************
int i4_wrap ( int ival, int ilo, int ihi )
{
  int jhi;
  int jlo;
  int value;
  int wide;

  if ( ilo <= ihi )
  {
    jlo = ilo;
    jhi = ihi;
  }
  else
  {
    jlo = ihi;
    jhi = ilo;
  }

  wide = jhi + 1 - jlo;

  if ( wide == 1 )
  {
    value = jlo;
  }
  else
  {
    value = jlo + i4_modp ( ival - jlo, wide );
  }

  return value;
}

//****************************************************************************
double *initial_condition ( int nx, double x[] )
{
  int i;
  double *u;
  double  pi = 4*atan(1.0);

  u = new double[nx];

  for ( i = 0; i < nx; i++ )
  {
    if  ( 40 <= x[i] && x[i] <= 100 )
    {
      u[i] = 80 * sin( (pi/60)*(x[i] - 40) );
    }
    else
    {
      u[i] = 0.0;
    }
  }
  return u;
}

//****************************************************************************
double *r8vec_linspace_new ( int n, double a_first, double a_last )
{
  double *a;
  int i;

  a = new double[n];

  if ( n == 1 )
  {
    a[0] = ( a_first + a_last ) / 2.0;
  }
  else
  {
    for ( i = 0; i < n; i++ )
    {
      a[i] = ( ( double ) ( n - 1 - i ) * a_first
             + ( double ) (         i ) * a_last )
             / ( double ) ( n - 1     );
    }
  }
  return a;
}

//****************************************************************************
void timestamp ( )
{
# define TIME_SIZE 40

  static char time_buffer[TIME_SIZE];
  const struct std::tm *tm_ptr;
  size_t len;
  std::time_t now;

  now = std::time ( NULL );
  tm_ptr = std::localtime ( &now );

  len = std::strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm_ptr );

  std::cout << time_buffer << "\n";

  return;
# undef TIME_SIZE
}
