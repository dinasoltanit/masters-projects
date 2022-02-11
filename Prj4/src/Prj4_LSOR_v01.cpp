# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <ctime>
# include <cmath>

//  Language: C++
//  Course: Advanced CFD1 - by Dr. Taeibi Rahni
//  Sharif University of Technology
//
//  Licence: GNU LGPL
//
//  Equation: ut = -c * ux
//  Method: Finite difference, FTCS, Explicit
//
//    The finite difference form is as bellow:
//
//      U(X,T+dt) - U(X,T)          ( U(X-dx,T+dt) + U(X+dx,T+dt) )
//      ------------------  = -c * --------------------------------- + F(X, T+dt)
//               dt                              2 * dx
//
//            +     c * dt / dx / 2    * U(X+dx,T+dt)
//            +     1                  * U(X,   T+dt)
//            -     c * dt / dx / 2    * U(X-dx,T+dt)
//      =                                U(X,   T)
//      +                        dt    * F(X,   T+dt)

using namespace std;
//****************************************************************************
// Functions: Begin
int main ( );
void dtable_data_write ( ofstream &output, int m, int n, double table[] );

void dtable_write ( string output_filename, int m, int n, double table[],
  bool header );

void f ( double a, double b, double t0, double t, int n, double x[],
  double value[] );

int r83_np_fa ( int n, double a[] );

double *r83_np_sl ( int n, double a_lu[], double b[], int job );

void timestamp ( );

void u0 ( double a, double b, double t0, int n, double x[], double value[] );

double ua ( double a, double b, double t0, double t );
double ub ( double a, double b, double t0, double t );
// Functions: End
//****************************************************************************

int main ( )
//****************************************************************************
{
  double *a;
  double *b;
  double *fvec;
  bool header;
  int i;
  int j;
  int job;
  double k;
  double *t;
  double dt;
  string t_file;
  double t_max;
  double t_min;
  int nt;
  double *u;
  string u_file;
  double w;
  double *x, *y;
  double dx, dy;
  string x_file;
  double x_max, y_max;
  double x_min, y_min;
  int nx, ny;
  string final_file;
  double c; // wave speed
  timestamp ( );

//****************************************************************************
//  Set X and Y values.
  x_min = 0.02;
  x_max = 0.04;
  y_min = 0.00;
  y_max = 0.02;
  dx = dy = 0.0025 // m
  nx = ((x_max - x_min) / dx) + 1;
  nx = ((y_max - y_min) / dy) + 1;
  x = new double[nx];
  y = new double[ny];

  for ( i = 0; i < nx; i++ )
  {
    x[i] = ( ( double ) ( nx - i - 1 ) * x_min
           + ( double ) (         i  ) * x_max )
           / ( double ) ( nx     - 1 );
  }

  for ( i = 0; i < ny; i++ )
  {
    y[i] = ( ( double ) ( ny - i - 1 ) * y_min
           + ( double ) (         i  ) * y_max )
           / ( double ) ( ny     - 1 );
  }

//****************************************************************************
//  Setting boundary conditions:

  u = new double[nx*nt];
  // returning the initial condition at the starting time:
  u0 ( x_min, x_max, t_min, nx, x, u );

//****************************************************************************
//  Defining the matrix A. This matrix is constant and does not change with time.
//  We can set it once, factor it once, and solve repeatedly.
//  Related to the left-hand-side of the system of equations:

  w = c* dt /dx /2; // w = stab as the stability checker.

  a = new double[3*nx];

  a[0+0*3] = 0.0;
  a[1+0*3] = 1.0;
  a[0+1*3] = 0.0;

  for ( i = 1; i < nx - 1; i++ )
  {
    a[2+(i-1)*3] =           - w;
    a[1+ i   *3] =           1.0;
    a[0+(i+1)*3] =           + w;
  }

  a[2+(nx-2)*3] =           -2*w;
  a[1+(nx-1)*3] =      1.0 + 2*w;
  a[2+(nx-1)*3] =            0.0;

//  Factor the matrix.

  // factoring the R83 system without pivoting:
  r83_np_fa ( nx, a );

//****************************************************************************
//  Related to the right-hand-side of the system of equations:

  b = new double[nx];
  fvec = new double[nx];

  for ( j = 1; j < nt; j++ )
  {

//  Set the right hand side B.
//
    b[0] = ua ( x_min, x_max, t_min, t[j] ); // returns the Dirichlet boundary condition at the left endpoint.

    // returning the right hand side of the 1d wave equation:
	// for this case, the excitation term is set to be zero;

	f ( x_min, x_max, t_min, t[j-1], nx, x, fvec );

    for ( i = 1; i < nx; i++ )
    {
      b[i] = u[i+(j-1)*nx] + dt * fvec[i];
    }

    //b[nx-1] = ub ( x_min, x_max, t_min, t[j] ); //returns the Dirichlet boundary condition at the right endpoint.

    delete [] fvec;

    job = 0;
	// solveing the R83 system factored by R83_NP_FA:
    fvec = r83_np_sl ( nx, a, b, job );

    for ( i = 0; i < nx; i++ )
    {
      u[i+j*nx] = fvec[i];
    }
  }
//****************************************************************************
// data handelling:

  x_file = "x.txt";
  header = false;
  dtable_write ( x_file, 1, nx, x, header ); // writes information to a DTABLE file.

  cout << "\n";
  cout << "  X data written to \"" << x_file << "\".\n";

  t_file = "t.txt";
  header = false;
  dtable_write ( t_file, 1, nt, t, header ); // writes information to a DTABLE file.

  cout << "  T data written to \"" << t_file << "\".\n";

  u_file = "u.txt";
  header = false;
  dtable_write ( u_file, nx, nt, u,  header ); // writes information to a DTABLE file.

  final_file = "EuImp_x2401_t100_ed2.txt";
  header = false;
  dtable_write ( final_file, nx, nt, u, header ); // writes information to a DTABLE file.

  cout << "  Final data written to \"" << final_file << "\".\n";

  delete [] a;
  delete [] b;
  delete [] fvec;
  delete [] t;
  delete [] u;
  delete [] x;
//
//  Terminate.
//
  cout << "\n";
  cout << "  Finite Difference: Euler Implicit for 1D Wave Equation\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************

void dtable_data_write ( ofstream &output, int nx, int nt, double table[] )

{
  int i;
  int j;
  //double *vect;
  double x_min = 0;
  double x_max = 300;
  double t_min = 0.0;
  double t_max = 0.6;
  double dx;
  //double dt;
  double vect[nt][nx];
  dx = ( x_max - x_min ) / ( double ) ( nx - 1 );
  //dt = ( t_max - t_min ) / ( double ) ( nt - 1 );

  for ( j = 0; j < nt; j++ )
  {
    output << "zone\n";
	for ( i = 0; i < nx; i++ )
    {
      vect[j][i] = table[i+j*nx];
	  output << x_min + i*dx << '\t' << vect[j][i] << '\n';
	  //output << setw(10) << table[i+j*dx] << "  ";
    }
    output << "\n";
  }

  return;
}
//****************************************************************************

void dtable_write ( string output_filename, int m, int n, double table[],
  bool header )

{
  ofstream output;

  output.open ( output_filename.c_str ( ) );

  if ( !output )
  {
    cerr << "\n";
    cerr << "DTABLE_WRITE - Fatal error!\n";
    cerr << "  Could not open the output file.\n";
    return;
  }

  if ( header )
  {
//  dtable_header_write ( output_filename, output, m, n );
  }

  dtable_data_write ( output, m, n, table );

  output.close ( );

  return;
}
//****************************************************************************

void f ( double a, double b, double t0, double t, int n, double x[],
  double value[] )

{
  int i;

  for ( i = 0; i < n; i++ )
  {
    value[i] = 0.0;
  }
  return;
}
//****************************************************************************

int r83_np_fa ( int n, double a[] )

{
  int i;

  for ( i = 1; i <= n-1; i++ )
  {
    if ( a[1+(i-1)*3] == 0.0 )
    {
      cout << "\n";
      cout << "R83_NP_FA - Fatal error!\n";
      cout << "  Zero pivot on step " << i << "\n";
      return i;
    }
//
//  Store the multiplier in L.
//
    a[2+(i-1)*3] = a[2+(i-1)*3] / a[1+(i-1)*3];
//
//  Modify the diagonal entry in the next column.
//
    a[1+i*3] = a[1+i*3] - a[2+(i-1)*3] * a[0+i*3];
  }

  if ( a[1+(n-1)*3] == 0.0 )
  {
    cout << "\n";
    cout << "R83_NP_FA - Fatal error!\n";
    cout << "  Zero pivot on step " << n << "\n";
    return n;
  }

  return 0;
}
//****************************************************************************

double *r83_np_sl ( int n, double a_lu[], double b[], int job )

{
  int i;
  double *x;

  x = new double[n];

  for ( i = 0; i < n; i++ )
  {
    x[i] = b[i];
  }

  if ( job == 0 )
  {
//
//  Solve L * Y = B.
//
    for ( i = 1; i < n; i++ )
    {
      x[i] = x[i] - a_lu[2+(i-1)*3] * x[i-1];
    }
//
//  Solve U * X = Y.
//
    for ( i = n; 1 <= i; i-- )
    {
      x[i-1] = x[i-1] / a_lu[1+(i-1)*3];
      if ( 1 < i )
      {
        x[i-2] = x[i-2] - a_lu[0+(i-1)*3] * x[i-1];
      }
    }
  }
  else
  {
//
//  Solve U' * Y = B
//
    for ( i = 1; i <= n; i++ )
    {
      x[i-1] = x[i-1] / a_lu[1+(i-1)*3];
      if ( i < n )
      {
        x[i] = x[i] - a_lu[0+i*3] * x[i-1];
      }
    }
//
//  Solve L' * X = Y.
//
    for ( i = n-1; 1 <= i; i-- )
    {
      x[i-1] = x[i-1] - a_lu[2+(i-1)*3] * x[i];
    }
  }

  return x;
}
//****************************************************************************

void timestamp ( )

{
# define TIME_SIZE 40

  static char time_buffer[TIME_SIZE];
  const struct tm *tm;
  time_t now;

  now = time ( NULL );
  tm = localtime ( &now );

  strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

  cout << time_buffer << "\n";

  return;
# undef TIME_SIZE
}
//****************************************************************************

void u0 ( double a, double b, double t0, int nx, double x[], int ny, double y[], double value[] )

{
  int i;
  double  pi = 4*atan(1.0);

  for ( i = 0; i < nx; i++ )
  {
    if  ( 40 <= x[i] && x[i] <= 100 )
    {
      value[i] = 80 * sin( (pi/60)*(x[i] - 40) );
    }
    else
    {
      value[i] = 0.0;
    }
  }

  return;
}
//****************************************************************************

double ua ( double a, double b, double t0, double t )

{
  double value;

  value = 0;

  return value;
}
//****************************************************************************

double ub ( double a, double b, double t0, double t )

{
  double value;

  value = 2.14618;

  return value;
}
