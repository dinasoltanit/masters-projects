# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <ctime>
# include <cmath>


//****************************************************************************
// This code is used both for the GS and the SOR Methods
//****************************************************************************


using namespace std;
//****************************************************************************
// Functions: Begin
int main ( );

void timestamp ( );

void u0 ( double x_min, double x_max, double y_min, double y_max, int nx, int ny, double x[], double y[], double **u );
void psi0 ( double x_min, double x_max, double y_min, double y_max, int nx, int ny, double x[], double y[], double **psi );

void Point_SOR(double **u, double **gs, double dx, double dy, double x_min, double y_min, double x_max, double y_max, double Beta, double w, double delta_v, int nx, int ny, double x[], double y[]);

void MyWrite(double ** phi, ofstream &fout, string output_filename, double dx, double dy, double x_min, double y_min, double x_max, double y_max, int nx, int ny);
// Functions: End
//****************************************************************************

int main ( )
//****************************************************************************
{
  string u_data_filename = "Laplace_phi_PGS.plt";
  ofstream u_data_unit;
  string psi_data_filename = "Laplace_psi_PGS.plt";
  ofstream psi_data_unit;
  bool header;
  double omega;
  double error_treshold;
  int iterations;
  double sor_adjustment;
  double delta_v;
  int i,j,row,col;
  int job;
  double k;
  double x_max = -0.02;
  double x_min = -0.04;
  double y_min = 0.00;
  double y_max = 0.02;
  double dx = 0.00125;
  double dy = 0.00125;
  double *x, *y;
  double **u;
  double **gs;
  double **u_psi;
  double **gs_psi;
  double w;
  int nx, ny;
  double Beta;
  timestamp ( );

//****************************************************************************
//  Set X and Y values.

  nx = ((x_max - x_min) / dx) + 1;
  ny = ((y_max - y_min) / dy) + 1;
  cout << "nx =" << nx << '\n';
  cout << "ny =" << ny << '\n';

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
//  New u**, gs**, u_psi**, gs_psi:

    u = new double *[nx];
    //cout << "here10\n";
    for (size_t i = 0; i < nx; i++)
    {
        u[i] = new double[ny];
    }

    gs = new double *[nx];
    //cout << "here10\n";
    for (size_t i = 0; i < nx; i++)
    {
        gs[i] = new double[ny];
    }

    u_psi = new double *[nx];
    //cout << "here10\n";
    for (size_t i = 0; i < nx; i++)
    {
        u_psi[i] = new double[ny];
    }

    gs_psi = new double *[nx];
    //cout << "here10\n";
    for (size_t i = 0; i < nx; i++)
    {
        gs_psi[i] = new double[ny];
    }
//****************************************************************************
//  Setting boundary conditions:

    u0 ( x_min, x_max, y_min, y_max, nx, ny, x,y, u  );

    u0 ( x_min, x_max, y_min, y_max, nx, ny, x,y, gs  );

    psi0 ( x_min, x_max, y_min, y_max, nx, ny, x,y, u_psi  );

    psi0 ( x_min, x_max, y_min, y_max, nx, ny, x,y, gs_psi  );

//****************************************************************************
//  Solution based on LSOR with omega = 1.2:
  omega = 1.0;
  error_treshold = 1E-7;
  iterations = 0;
  sor_adjustment = 0.0;
  delta_v = 0.0;
  Beta = 1.0;
  Point_SOR(u, gs, dx, dy, x_min, y_min, x_max, y_max, Beta, omega, delta_v, nx, ny, x, y);

  Point_SOR(u_psi, gs_psi, dx, dy, x_min, y_min, x_max, y_max, Beta, omega, delta_v, nx, ny, x, y);

//****************************************************************************

//****************************************************************************
// data handelling:

  MyWrite(u, u_data_unit, u_data_filename, dx, dy, x_min, y_min, x_max, y_max, nx, ny);

  MyWrite(u_psi, psi_data_unit, psi_data_filename, dx, dy, x_min, y_min, x_max, y_max, nx, ny);


  cout << "number of iterations = " << iterations << '\n';
  cout << "Calculation, done.\n";

  delete [] x;
  delete [] y;
  delete [] u;
  delete [] gs;
  delete [] u_psi;
  delete [] gs_psi;
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
//u0 ( x_min, x_max, y_min, y_max, nx, ny, x,y, u );
void u0 ( double x_min, double x_max, double y_min, double y_max, int nx, int ny, double x[], double y[], double **u )

{
  int i;
  double  pi = 4.0*atan(1.0);
//cout << "here12-1\n";
  for ( i = 0; i < nx; i++ )
  {
    // y = 0.02 m
    u[i][ny-1] = ( 0.004 + 5.0 * (x[i] *x[i]) ) / ( sqrt(1.0 + 0.0004/(x[i] *x[i]) ) ) / ( sqrt( 0.0004 + (x[i] *x[i]) ) );
    //u[i][ny-1] = 5.0 * ( 0.0004 + (0.0004 + (x[i] *x[i]) ) ) / sqrt(1.0 + 0.0004 * (1.0/(x[i] *x[i])) ) / sqrt(0.0004 + (x[i] *x[i]));
    //cout << u[i][ny-1];
    //gs[i][ny] = 5.0 * (0.0004 + (0.0004 + (x[i] *x[i]) ) / sqrt(1 + 0.0004 * (1/(x[i] *x[i])) ) / sqrt(0.0004 + (x[i] *x[i])));
  }
//cout << "here12-2\n";
  for ( i = 0; i < nx; i++ )
  {
    // y = 0.0 m
    u[i][0] = ( 0.002 + 5.0 * (x[i] *x[i]) )/( sqrt(x[i] *x[i]) );
    //u[i][0] = 5.0 * (0.0004 + ((x[i] *x[i]))) / sqrt((x[i] *x[i]));
    //gs[i][0] = 5.0 * (0.0004 + ((x[i] *x[i]))) / sqrt((x[i] *x[i]));
  }
//cout << "here12-3\n";
  for ( i = 0; i < ny; i++ )
  {
    // x = -0.04 m
    u[0][i] = ( 0.01 + 5.0 * (y[i] * y[i]) ) / ( sqrt( 0.0016 + (y[i] * y[i]) ) ) / ( sqrt( 1.0 + 625.0 * (y[i] * y[i]) ) );
    //u[nx-1][i] = 5.0 * (0.0004 + (0.0016 + (y[i] * y[i]))) / (sqrt(1 + 625.0 * (y[i] * y[i])) * sqrt(0.0016 + (y[i] * y[i])));
    //gs[nx][i] = 5.0 * (0.0004 + (0.0016 + (y[i] * y[i]))) / sqrt(1 + 625.0 * (y[i] * y[i])) / sqrt(0.0016 + (y[i] * y[i]));
  }
//cout << "here12-4\n";
  for ( i = 0; i < ny; i++ )
  {
    // x = -0.02 m
    u[nx-1][i] = ( 0.004 + 5.0 * (y[i] * y[i]) ) / ( sqrt(0.0004 + (y[i] * y[i])) ) / ( sqrt(1.0 + 2500.0 * (y[i] * y[i])) );
    //u[0][i] = 5.0 * (0.0004 + (0.0004 + (y[i] * y[i]))) / sqrt(1.0 + 2500.0 * (y[i] * y[i])) / sqrt(0.0004 + ((y[i] * y[i])));
    //cout << "u[0][" << i << "] = " << u[0][i] << '\n';
    //gs[0][i] = 5.0 * (0.0004 + (0.0004 + (y[i] * y[i]))) / sqrt(1 + 2500.0 * (y[i] * y[i])) / sqrt(0.0004 + ((y[i] * y[i])));
  }
//cout << "here12-5\n";
  //return;
}
//****************************************************************************
//void psi0 ( double x_min, double x_max, double y_min, double y_max, int nx, int ny, double x[], double y[], double **psi );
void psi0 ( double x_min, double x_max, double y_min, double y_max, int nx, int ny, double x[], double y[], double **psi )

{
  int i;
  double  pi = 4.0*atan(1.0);
//cout << "here12-1\n";
  for ( i = 0; i < nx; i++ )
  {
    // y = 0.02 m
    psi[i][ny-1] = ( 0.1 * x[i] ) / ( sqrt(1.0 + 0.0004/(x[i] *x[i]) ) ) / ( sqrt( 0.0004 + (x[i] *x[i]) ) );
  }
//cout << "here12-2\n";
  for ( i = 0; i < nx; i++ )
  {
    // y = 0.0 m
    psi[i][0] = 0.0;
  }
//cout << "here12-3\n";
  for ( i = 0; i < ny; i++ )
  {
    // x = -0.04 m
    psi[0][i] = ( ( -0.2 * y[i] ) * (0.0012 + (y[i] * y[i])) *( sqrt( 1.0 + 625.0 * (y[i] * y[i]) ) ) ) / ( 0.0016 + (y[i] * y[i]) ) / ( sqrt(0.0016 + (y[i] * y[i])) );
  }
//cout << "here12-4\n";
  for ( i = 0; i < ny; i++ )
  {
    // x = -0.02 m
    psi[nx-1][i] = ( ( -0.1 * (y[i] * y[i] * y[i])) * (sqrt(1.0 + 2500.0 * (y[i] * y[i]))) ) / ( 0.0004 + (y[i] * y[i]) ) / ( sqrt(0.0004 + (y[i] * y[i])) );
  }
//cout << "here12-5\n";
  //return;
}

void Point_SOR(double **u, double **gs, double dx, double dy, double x_min, double y_min, double x_max, double y_max, double Beta, double w, double delta_v, int nx, int ny, double x[], double y[])
{

  double sor_adjustment = 0.0;
  double omega = w;
  int iterations = 0;
  double error_treshold = 1E-7;
  while(true)
        {
            for ( int row = 1; row < nx-1; row++ )
            {
                for ( int col = 1; col < ny-1; col++ )
                {
                    gs[row][col] = (gs[row-1][col] + gs[row+1][col] + gs[row][col-1] + gs[row][col+1])/4;
                    //sor_adjustment = (1.0/omega) * (gs[row][col] - u[row][col]);
                    //u[row][col] = sor_adjustment + u[row][col];
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
                for ( int row = 1; row < nx-1; row++ )
                {
                    for ( int col = 1; col < ny-1; col++ )
                    {
                        u[row][col] = gs[row][col];
                    }
                }

            }
        }
        cout << "number of iterations = " << iterations << '\n';

}

void MyWrite(double ** phi, ofstream &fout, string output_filename, double dx, double dy, double x_min, double y_min, double x_max, double y_max, int nx, int ny)
{
    double x = x_min;
    double y = y_min;
    fout.open( output_filename.c_str ( ) );
    cout << "file opened.";

    if( !fout ) { // file couldn't be opened
        cerr << "Error: file could not be opened" << endl;
        exit(1);
     }
    fout << defaultfloat << "VARIABLE= \"I\", \"J\", \"Phi\"\n" << "ZONE I=\t" << nx << " , J=\t" << ny << " , F=point"<< endl;
    for (int i = 0; i < nx; i++)
    {
        y = y_min;
        for (int j = 0; j < ny; j++)
        {
            fout << scientific << setprecision(10) << x << "\t" << y << "\t" << phi[i][j] << endl;
            y += dy;
        }
        x += dx;
    }
  //cout << "number of iterations = " << iterations << '\n';

  fout.close();
  cout << "file closed.\n";
}

//*************************************temp*************************************************************

