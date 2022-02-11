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

//****************************************************************************80

int main ( )

//****************************************************************************80
/*
Language: C++
Author: Dina Soltani Tehrani.
Course: Advanced CFD1 - by Dr. Taeibi Rahni
Sharif University of Technology

Equation: ut = -c * ux
Method: Finite difference, FTCS, Explicit
*/
{
  double a;
  double b;
  double c;
  double L;
  double duration;
  string command_filename = "advection_commands.txt";
  ofstream command_unit;
  string data_filename = "advection_data.txt";
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

  timestamp ( );
  cout << "\n";
  cout << "Finite Difference: FTCS for 1D Wave Equation:\n";
  cout << "  C++ language\n";
  cout << "\n";
  cout << "  Solve the constant-velocity advection equation in 1D,\n";
  cout << "    du/dt = - c du/dx\n";
  cout << "  over the interval:\n";
  cout << "    0.0 <= x <= 300.0\n";
  cout << "  with periodic boundary conditions, and\n";
  cout << "  with a given initial condition\n";
  cout << "    u(0,x) = 80 * sin( (pi/60)* (x-40) ) for 40.0 <= x <= 100.0\n";
  cout << "           = 0 elsewhere.\n";
  cout << "\n";
  cout << "  We use a method known as FTCS:\n";
  cout << "   FT: Forward Time  : du/dt = (u(t+dt,x)-u(t,x))/dt\n";
  cout << "   CS: Centered Space: du/dx = (u(t,x+dx)-u(t,x-dx))/2/dx\n";
  cout << "\n";
  cout << "  Resolution:\n";
  cout << "   Space: 301, 601, 3001\n";
  cout << "   Time: 1000 iterations\n";

  L = 300; // Length (m)
  nx = 601; // Mid Space Resolution
  dx = L / ( double ) ( nx - 1 ); // size-step
  a = 0.0; // x-start
  b = 300.0; // x-end
  x = r8vec_linspace_new ( nx, a, b );
  nt = 1000; // Time Resolution
  duration = 1.0; // Duration (s)
  dt = duration / ( double ) ( nt ); // time-step
  c = 1.0;

  u = initial_condition ( nx, x ); // Initializing old velocity with initial condition.
//
//  ***** Openning data file, and writing solutions as they are computed. *****
//
  data_unit.open ( data_filename.c_str ( ) );

  t = 0.0; // Writing the initial data.
  data_unit << "  " << x[0]
            << "  " << t
            << "  " << u[0] << "\n";
  // Writing data as we move forward in space.
  for ( j = 0; j < nx; j++ )
  {
    data_unit << "  " << x[j]
              << "  " << t
              << "  " << u[j] << "\n";
  }
  data_unit << "\n";

  nt_step = 100;

  cout << "\n";
  cout << "  Number of nodes NX = " << nx << "\n";
  cout << "  Number of time steps NT = " << nt << "\n";
  cout << "  Constant velocity C = " << c << "\n";

  unew = new double[nx]; // Creating the new velocity vector

// ***** Core Calculation: begin *****
  for ( i = 0; i < nt; i++ )
  {
    // First do this for all x in xn:
    for ( j = 0; j < nx; j++ )
    {
      jm1 = i4_wrap ( j - 1, 0, nx - 1 ); // phi_n_j-1
      jp1 = i4_wrap ( j + 1, 0, nx - 1 ); // phi_n_j+1
      unew[j] = u[j] - c * dt / dx / 2.0 * ( u[jp1] - u[jm1] );
    }
    // Then, do this, again for all x in xn:
    for ( j = 0; j < nx; j++ )
    {
      u[j] = unew[j];
    }
    if ( i == nt_step - 1 )
    {
      t = ( double ) ( i ) * dt;
      for ( j = 0; j < nx; j++ )
      {
        data_unit << "  " << x[j]
                  << "  " << t
                  << "  " << u[j] << "\n";
      }
      data_unit << "\n";
      nt_step = nt_step + 100;
    }
  }
// ***** Core Calculation: end *****


//
//  ***** Closing the data file as the computation is done. *****
//
  data_unit.close ( );

    cout << "\n";
  cout << "  Plot data written to the file \"" << data_filename << "\"\n";
//
//  Write gnuplot command file.
//
  command_unit.open ( command_filename.c_str ( ) );

  command_unit << "set term png\n";
  command_unit << "set output 'advection.png'\n";
  command_unit << "set grid\n";
  command_unit << "set style data lines\n";
  command_unit << "unset key\n";
  command_unit << "set xlabel '<---X--->'\n";
  command_unit << "set ylabel '<---Time--->'\n";
  command_unit << "splot '" << data_filename << "' using 1:2:3 with lines\n";
  command_unit << "quit\n";

  command_unit.close ( );

  cout << "  Gnuplot command data written to the file \"" << command_filename << "\"\n";

//  Free memory.
//
  delete [] u;
  delete [] unew;
  delete [] x;
//
//  Terminate.
//
  cout << "\n";
  cout << "FFinite Difference: FTCS for 1D Wave Equation\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}

//
// ***** Functions *****

int i4_modp ( int i, int j )

//****************************************************************************
//
//  Purpose:
//
//    I4_MODP returns the nonnegative remainder of I4 division.
//
//  Parameters:
//
//    Input, int I, the number to be divided.
//
//    Input, int J, the number that divides I.
//
//    Output, int I4_MODP, the nonnegative remainder when I is
//    divided by J.
//
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

//****************************************************************************
//
//  Purpose:
//
//    I4_WRAP forces an I4 to lie between given limits by wrapping.
//
//  Parameters:
//
//    Input, int IVAL, an integer value.
//
//    Input, int ILO, IHI, the desired bounds for the integer value.
//
//    Output, int I4_WRAP, a "wrapped" version of IVAL.
//
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

//****************************************************************************
//
//  Purpose:
//
//    INITIAL_CONDITION sets the initial condition.
//
//  Parameters:
//
//    Input, int NX, the number of nodes.
//
//    Input, double X[NX], the coordinates of the nodes.
//
//    Output, double INITIAL_CONDITION[NX], the value of the initial condition.
//
{
  int i;
  double *u;

  u = new double[nx];

  for ( i = 0; i < nx; i++ )
  {
    if  ( 0.4 <= x[i] && x[i] <= 0.6 )
    {
      u[i] = pow ( 10.0 * x[i] - 4.0, 2 )
           * pow ( 6.0 - 10.0 * x[i], 2 );
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

//****************************************************************************
//
//  Purpose:
//
//    R8VEC_LINSPACE_NEW creates a vector of linearly spaced values.
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input, double A_FIRST, A_LAST, the first and last entries.
//
//    Output, double R8VEC_LINSPACE_NEW[N], a vector of linearly spaced data.
//
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

//****************************************************************************
//
//  Purpose:
//
//    TIMESTAMP prints the current YMDHMS date as a time stamp.
//
//  Example:
//
//    31 May 2001 09:45:54 AM
//
//  Parameters:
//
//    None
//
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
