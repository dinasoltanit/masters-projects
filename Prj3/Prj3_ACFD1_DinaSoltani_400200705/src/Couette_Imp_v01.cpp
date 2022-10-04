# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <ctime>
# include <cmath>

using namespace std;
//****************************************************************************
// FUNCTIONS BEGIN
//****************************************************************************
int main ( );
void dtable_data_write ( ofstream &output, int ny, int nt, double table[] );
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

//****************************************************************************
// MAIN PROGRAM
//****************************************************************************
int main ( )

{
  double *a;
  double *b;
  double *fvec;
  bool header;
  int i,j;
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
  double *y;
  double dy;
  string final_file;
  string y_file;
  double y_max;
  double y_min;
  int ny;

  timestamp ( );

  k = 0.01307;
//
//  Set X values.
//
  y_min = 0.0;
  y_max = 0.05;
  ny = 81;
  dy = ( y_max - y_min ) / ( double ) ( ny - 1 );

  y = new double[ny];

  for ( i = 0; i < ny; i++ )
  {
    y[i] = ( ( double ) ( ny - i - 1 ) * y_min
           + ( double ) (         i     ) * y_max )
           / ( double ) ( ny     - 1 );
  }
//
//  Set T values.
//
  t_min = 0.0;
  t_max = 2.0;
  nt = 2001;
  dt = ( t_max - t_min ) / ( double ) ( nt - 1 );

  t = new double[nt];

  for ( j = 0; j < nt; j++ )
  {
    t[j] = ( ( double ) ( nt - j - 1 ) * t_min
           + ( double ) (         j     ) * t_max )
           / ( double ) ( nt     - 1 );
  }
//
//  Set the initial data, for time T_MIN.
//
  u = new double[ny*nt];

  u0 ( y_min, y_max, t_min, ny, y, u);
//
//  The matrix A does not change with time.  We can set it once,
//  factor it once, and solve repeatedly.
//
  w = k * dt / dy / dy;

  a = new double[3*ny];

  a[0+0*3] = 0.0;

  a[1+0*3] = 1.0;
  a[0+1*3] = 0.0;

  for ( i = 1; i < ny - 1; i++ )
  {
    a[2+(i-1)*3] =           - w;
    a[1+ i   *3] = 1.0 + 2.0 * w;
    a[0+(i+1)*3] =           - w;
  }

  a[2+(ny-2)*3] = 0.0;
  a[1+(ny-1)*3] = 1.0;

  a[2+(ny-1)*3] = 0.0;
//
//  Factor the matrix.
//
  r83_np_fa ( ny, a );

  b = new double[ny];
  fvec = new double[ny];

  for ( j = 1; j < nt; j++ )
  {
//
//  Set the right hand side B.
//
    b[0] = ua ( y_min, y_max, t_min, t[j] );

    f ( y_min, y_max, t_min, t[j-1], ny, y, fvec );

    for ( i = 1; i < ny - 1; i++ )
    {
      b[i] = u[i+(j-1)*ny] + dt * fvec[i];
    }

    b[ny-1] = ub ( y_min, y_max, t_min, t[j] );

    delete [] fvec;

    job = 0;
    fvec = r83_np_sl ( ny, a, b, job );

    for ( i = 0; i < ny; i++ )
    {
      u[i+j*ny] = fvec[i];
    }
  }


  final_file = "Couette_Imp_y81_t2001_visc100_bc3.txt";
  header = false;
  dtable_write ( final_file, ny, nt, u, header ); // writes information to a DTABLE file.


  delete [] a;
  delete [] b;
  delete [] fvec;
  delete [] t;
  delete [] u;
  delete [] y;
//
//  Terminate.
//
  cout << "\n";
  cout << "FD1D_HEAT_IMPLICIT\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";
  timestamp ( );

  return 0;
}
//****************************************************************************80

void dtable_data_write ( ofstream &output, int ny, int nt, double table[] )
{
  int i;
  int j;
  //double *vect;
  double y_min = 0.0;
  double y_max = 0.05;
  double t_min = 0.0;
  double t_max = 20.0;
  double dy;
  //double dt;
  double vect[nt][ny];
  dy = ( y_max - y_min ) / ( double ) ( ny - 1 );
  //dt = ( t_max - t_min ) / ( double ) ( nt - 1 );

  for ( j = 0; j < nt; j++ )
  {
    output << "zone\n";
	for ( i = 0; i < ny; i++ )
    {
      vect[j][i] = table[i+j*ny];
	  output << y_min + i*dy << '\t' << vect[j][i] << '\n';
	  //output << setw(10) << table[i+j*dx] << "  ";
    }
    output << "\n";
  }

  return;
}
//****************************************************************************80

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
//****************************************************************************80

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
//****************************************************************************80

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
//****************************************************************************80

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
//****************************************************************************80

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
//****************************************************************************80

void u0 ( double a, double b, double t0, int n, double x[], double value[] )
{
  int i;
  double U0 = 1.5;

  for ( i = 0; i < (n-1); i++ )
  {
    value[i] = 0.0;
  }
  value[n-1] = U0 * 0.0;
  return;
}
//****************************************************************************80

double ua ( double a, double b, double t0, double t )
{
  double value;
  double U0 = 1.5;

  value = 0.0;

  return value;
}
//****************************************************************************80

double ub ( double a, double b, double t0, double t )
{
  double value;
  double U0 = 1.5;
  double pi = 4.0*atan(1.0);


  value = U0 * sin(2 * pi * t);


  return value;
}
