# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <ctime>
# include <cmath>

using namespace std;
//****************************************************************************
// Functions: Begin
int main ( );
void dtable_data_write ( ofstream &output, int m, int n, double table[] );

void dtable_write ( string output_filename, int m, int n, double table[],
  bool header );

void timestamp ( );

//void u0 ( double a, double b, double t0, int n, double x[], double value[] );

// Functions: End
//****************************************************************************

int main ( )
//****************************************************************************
{
  string data_filename = "Laplace_PJ_v01.txt";
  ofstream data_unit;
  bool header;
  double omega;
  double error_treshold;
  int iterations;
  double sor_adjustment;
  double delta_v;
  int i,j,row,col;
  int job;
  double k;
  double x_min = 0.02;
  double x_max = 0.04;
  double y_min = 0.00;
  double y_max = 0.02;
  double dx = 0.00125;
  double dy = 0.00125;
  double *x, *y;
//  double *u;
//  double *gs;
  string u_file;
  double w;
  string x_file;
  string y_file;
  int nx, ny;
  string final_file;
  timestamp ( );

//****************************************************************************
//  Set X and Y values.

  nx = ((x_max - x_min) / dx) + 1;
  ny = ((y_max - y_min) / dy) + 1;
  cout << "nx =" << nx << '\n';
  cout << "ny =" << ny << '\n';
  double u[nx][ny];
  double gs[nx][ny];
  double s[nx][ny];


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

//  u  = new double[nx][ny];
//  gs = new double[nx][ny]; // gs == gauss seidel
//u0 ( x_min, x_max, y_min, y_max, nx, ny, x,y, u  );
//u0 ( x_min, x_max, y_min, y_max, nx, ny, x,y, gs );

  for ( i = 0; i < nx; i++ )
  {
    u[i][ny] = 5.0 * (0.0004 + (0.0004 + (x[i] *x[i]) ) / sqrt(1 + 0.0004 * (1/(x[i] *x[i])) ) / sqrt(0.0004 + (x[i] *x[i])));
    gs[i][ny] = 5.0 * (0.0004 + (0.0004 + (x[i] *x[i]) ) / sqrt(1 + 0.0004 * (1/(x[i] *x[i])) ) / sqrt(0.0004 + (x[i] *x[i])));
    s[i][ny] = 5.0 * (0.0004 + (0.0004 + (x[i] *x[i]) ) / sqrt(1 + 0.0004 * (1/(x[i] *x[i])) ) / sqrt(0.0004 + (x[i] *x[i])));
  }

  for ( i = 0; i < nx; i++ )
  {
    u[i][0] = 5.0 * (0.0004 + ((x[i] *x[i]))) / sqrt((x[i] *x[i]));
    gs[i][0] = 5.0 * (0.0004 + ((x[i] *x[i]))) / sqrt((x[i] *x[i]));
    s[i][0] = 5.0 * (0.0004 + ((x[i] *x[i]))) / sqrt((x[i] *x[i]));
  }

  for ( i = 0; i < ny; i++ )
  {
    u[nx][i] = 5.0 * (0.0004 + (0.0016 + (y[i] * y[i]))) / sqrt(1 + 625.0 * (y[i] * y[i])) / sqrt(0.0016 + (y[i] * y[i]));
    gs[nx][i] = 5.0 * (0.0004 + (0.0016 + (y[i] * y[i]))) / sqrt(1 + 625.0 * (y[i] * y[i])) / sqrt(0.0016 + (y[i] * y[i]));
    s[nx][i] = 5.0 * (0.0004 + (0.0016 + (y[i] * y[i]))) / sqrt(1 + 625.0 * (y[i] * y[i])) / sqrt(0.0016 + (y[i] * y[i]));
  }

  for ( i = 0; i < ny; i++ )
  {
    u[0][i] = 5.0 * (0.0004 + (0.0004 + (y[i] * y[i]))) / sqrt(1 + 2500.0 * (y[i] * y[i])) / sqrt(0.0004 + ((y[i] * y[i])));
    gs[0][i] = 5.0 * (0.0004 + (0.0004 + (y[i] * y[i]))) / sqrt(1 + 2500.0 * (y[i] * y[i])) / sqrt(0.0004 + ((y[i] * y[i])));
    s[0][i] = 5.0 * (0.0004 + (0.0004 + (y[i] * y[i]))) / sqrt(1 + 2500.0 * (y[i] * y[i])) / sqrt(0.0004 + ((y[i] * y[i])));
  }
//****************************************************************************
//  Solution based on LSOR with omega = 1.4:
  //omega = 1;
  error_treshold = 1E-4;
  iterations = 0;
  //sor_adjustment = 0.0;
  delta_v = 0.0;
    while(true)
        {
            for ( row = 1; row < nx-1; row++ )
            {
                for ( col = 1; col < ny-1; col++ )
                {
                    gs[row][col] = (u[row-1][col] + u[row+1][col] + u[row][col-1] + u[row][col+1])/4;
                    delta_v += abs(u[row][col] - gs[row][col]);
                }
            }

            iterations += 1;

            if (delta_v < error_treshold)
            {
                break;
            }
            else
            {
                delta_v = 0;  // Restart counting delta_v for the next iteration
                // copy u to gs

                for ( row = 1; row < nx; row++ )
                {
                    for ( col = 1; col < nx-1; col++ )
                    {
                        u[row][col] = gs[row][col];
                    }
                }
            }
        }

//****************************************************************************

//****************************************************************************
// data handelling:

  /*
  x_file = "x.txt";
  header = false;
  dtable_write ( x_file, 1, nx, x, header ); // writes information to a DTABLE file.

  cout << "\n";
  cout << "  X data written to \"" << x_file << "\".\n";

  y_file = "y.txt";
  header = false;
  dtable_write ( y_file, 1, ny, y, header ); // writes information to a DTABLE file.

  cout << "  Ydata written to \"" << y_file << "\".\n";

  //u_file = "u.txt";
  //header = false;
  //dtable_write ( u_file, nx, ny, u,  header ); // writes information to a DTABLE file.

  //final_file = "LSOR_w14_v01.txt";
  //header = false;
  //dtable_write ( final_file, nx, ny, u, header ); // writes information to a DTABLE file.
  */
  data_unit.open( data_filename.c_str ( ) );
  cout << "file opened.";

  if( !data_unit ) { // file couldn't be opened
      cerr << "Error: file could not be opened" << endl;
      exit(1);
   }
   for ( row = 0; row < nx; row++ )
    {
        cout << "in row, data-out for" << '\n';
        for ( col = 0; col < ny; col++ )
        {
           cout << "in col, data-out for" << '\n';
           data_unit << x[row] << '\t' << y[col] << '\t' << u[row][col] << '\n';
           cout << "printed: x = " << x[row] << '\t' << y[col] << '\t' << "u = " << u[row][col] << '\n';
        }
    }

  cout << "number of iterations = " << iterations << '\n';
  cout << "Calculation, done.\n";
  data_unit.close();

  delete [] x;
  delete [] y;
//
//  Terminate.
//
//  cout << "\n";
//  cout << "  Finite Difference: Euler Implicit for 1D Wave Equation\n";
//  cout << "  Normal end of execution.\n";
//  cout << "\n";
  timestamp ( );

  return 0;
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
