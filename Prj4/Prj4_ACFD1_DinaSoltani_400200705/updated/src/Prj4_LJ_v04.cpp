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

void u0 ( double x_min, double x_max, double y_min, double y_max, int nx, int ny, double x[], double y[], double **u );
void psi0 ( double x_min, double x_max, double y_min, double y_max, int nx, int ny, double x[], double y[], double **psi );

void Line_SOR(double **u, double dx, double dy, double x_min, double y_min, double x_max, double y_max, double Beta, double w, int nx, int ny, double x[], double y[]);

void Initial_Guessed(double **phi, int node_number_x, int node_number_y);

void myError( double **u, double **phi_k, int nx, int ny, int &E);

void MyWrite(double ** phi, ofstream &fout, string output_filename, double dx, double dy, double x_min, double y_min, double x_max, double y_max, int nx, int ny);
// Functions: End
//****************************************************************************

int main ( )
//****************************************************************************
{
  string u_data_filename = "Laplace_phi_LJ.plt";
  ofstream u_data_unit;
  string psi_data_filename = "Laplace_psi_LJ.plt";
  ofstream psi_data_unit;
  bool header;
  double omega;
  double Beta;
  double error_treshold;
  int iterations;
  double sor_adjustment;
  double delta_v;
  int i,j,row,col;
  int job;
  double k;
  double **u;
  double **psi;
  //double **gs;
  string u_file;
  double w;
  double *x, *y;
  double dx, dy;
  string x_file;
  string y_file;
  double x_max, y_max;
  double x_min, y_min;
  int nx, ny;
  string final_file;
  timestamp ( );

//****************************************************************************
//  Set X and Y values.
  x_max = -0.02;
  x_min = -0.04;
  y_min = 0.00;
  y_max = 0.02;
  dx = dy = 0.00125; // m
  nx = ((x_max - x_min) / dx) + 1;
  cout << "nx = " << nx << '\n';
  ny = ((y_max - y_min) / dy) + 1;
  //cout << "here1\n";
  x = new double[nx];
  //cout << "here2\n";
  y = new double[ny];
  //cout << "here3\n";

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
//cout << "here4\n";
//****************************************************************************
//  Setting boundary conditions for phi:
    u = new double *[nx];
    //cout << "here10\n";
    for (size_t i = 0; i < nx; i++)
    {
        u[i] = new double[ny];
    }
//  u  = new double[nx][ny];
//  gs = new double[nx][ny]; // gs == gauss seidel
    u0 ( x_min, x_max, y_min, y_max, nx, ny, x,y, u  );
//  u0 ( x_min, x_max, y_min, y_max, nx, ny, x,y, gs );

//****************************************************************************
//  Setting boundary conditions for psi:
    psi = new double *[nx];
    //cout << "here10\n";
    for (size_t i = 0; i < nx; i++)
    {
        psi[i] = new double[ny];
    }

    psi0 ( x_min, x_max, y_min, y_max, nx, ny, x,y, psi  );

//****************************************************************************
//  Solution based on LSOR with omega = 1.4:
//  void Line_SOR(double **u, double dx, double dy, double x_min, double y_min, double x_max, double y_max, double Beta, double w, int nx, int ny)
  omega = 1;
  error_treshold = 1E-7;
  iterations = 0;
  Beta = 1.0;
  //cout << "here5\n";
  Line_SOR(u, dx, dy, x_min, y_min, x_max, y_max, Beta, omega, nx, ny, x, y);
  Line_SOR(psi, dx, dy, x_min, y_min, x_max, y_max, Beta, omega, nx, ny, x, y);
  //cout << "here6\n";
//****************************************************************************

  MyWrite(u, u_data_unit, u_data_filename, dx, dy, x_min, y_min, x_max, y_max, nx, ny);

  MyWrite(psi, psi_data_unit, psi_data_filename, dx, dy, x_min, y_min, x_max, y_max, nx, ny);

//****************************************************************************

//  cout << "number of iterations = " << iterations << '\n';
  cout << "Calculation, done.\n";
  delete [] u;
  delete [] psi;
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
//****************************************************************************
void Line_SOR(double **u, double dx, double dy, double x_min, double y_min, double x_max, double y_max, double Beta, double w, int nx, int ny, double x[], double y[])
{
    int E = 1;
    int k = 0;

    double * coefficient_JM1;
    double * coefficient_J;
    double * coefficient_JP1;
    double * constant;
    double * A;
    double * B;
    double ** phi_k;

    //cout << "here7\n";
    coefficient_JM1 = new double[ny - 2];
    //cout << "here8\n";
    coefficient_J = new double[ny - 2];
    coefficient_JP1 = new double[ny - 2];
    constant = new double[ny - 2];
    A = new double[ny - 2];
    B = new double[ny - 2];
    //cout << "here9\n";
    phi_k = new double *[nx];
    //cout << "here10\n";
    for (size_t i = 0; i < nx; i++)
    {
        phi_k[i] = new double[ny];
    }
//cout << "here11\n";
    u0 ( x_min, x_max, y_min, y_max, nx, ny, x,y, phi_k );

    //u0 ( x_min, x_max, y_min, y_max, nx, ny, x,y, u );
    //cout << "here12\n";
    //Boundry_Condition(phi_k, dx, dy, x_min, y_min, x_max, y_max, nx, ny);
    Initial_Guessed(phi_k, nx, ny);
    //cout << "here13\n";

    while (E)
    {
        //cout << "here14\n";
        E = 0;
        for (int i = 1; i < (nx - 1); i++)
        {
            for (int j = 1; j < (ny - 1); j++)
            {
                if (j == 1)
                {
                    //cout << "in if 1, j = 1; i =" << i << "\n";
                    coefficient_JM1[j - 1] = 0;
                    //cout << "in if, here1\n";
                    coefficient_J[j - 1] = 2 * (1 + Beta * Beta);
                    //cout << "in if, here2\n";
                    coefficient_JP1[j - 1] = 0;
                    //cout << "in if, here3\n";
                    //cout << w << '\n';
                    //cout << u[i - 1][j] << '\n';
                    //cout << phi_k[i + 1][j] << '\n';
                    //cout << u[i][j - 1] << '\n';
                    //cout << phi_k[i][j] << '\n';
                    constant[j - 1] = 1.0 * Beta * Beta * phi_k[i][j - 1] + 1.0 * Beta * Beta * phi_k[i][j + 1] + u[i - 1][j] + u[i + 1][j];
                    //cout << constant[j - 1] << '\n';
                    //cout << "in if, here4\n";
                }
                else if (j == (ny - 2))
                {
                    //cout << "in if 2, j = ny-2; i =" << i << "\n";
                    coefficient_JM1[j - 1] = 0;
                    coefficient_J[j - 1] = 2 * (1 + Beta * Beta);
                    coefficient_JP1[j - 1] = 0;
                    constant[j - 1] = 1.0 * u[i - 1][j] + 1.0 * u[i + 1][j] + 1.0 * Beta * Beta * phi_k[i][j + 1] + Beta * Beta * phi_k[i][j - 1];
                }
                else
                {
                    //cout << "in else, j = " << j << ", i =" << i << "\n";
                    coefficient_JM1[j - 1] = 0;
                    coefficient_J[j - 1] = 2 * (1 + Beta * Beta);
                    coefficient_JP1[j - 1] = 0;
                    constant[j - 1] = 1.0 * u[i - 1][j] + 1.0 * u[i + 1][j] + Beta * Beta * phi_k[i][j + 1] + Beta * Beta * phi_k[i][j - 1];
                }
            }
            A[0] = - coefficient_JP1[0] / coefficient_J[0];
            B[0] = constant[0] / coefficient_J[0];
            for (int j = 1; j < (ny - 2); j++)
            {
                double q = coefficient_J[j] + coefficient_JM1[j] * A[j - 1];
                A[j] = -1 * coefficient_JP1[j] / q;
                B[j] = (constant[j] - coefficient_JM1[j] * B[j - 1]) / q;
            }
            u[i][ny - 2] = B[ny - 3];
            for (int j = (ny - 3); j >= 1; j--)
            {
                u[i][j] = A[j - 1] * u[i][j + 1] + B[j - 1];
            }
        }
        //cout << "here15\n";
        myError(u, phi_k, nx, ny, E);
        //cout << "here16\n";
        k++;
    }

    for (int i = 0; i < nx; i++)
    {
        delete[] phi_k[i];
    }
    delete[] phi_k;

    delete[] coefficient_JM1;
    delete[] coefficient_J;
    delete[] coefficient_JP1;
    delete[] constant;
    delete[] A;
    delete[] B;
}

void Initial_Guessed(double **phi, int nx, int ny)
{
    for (int i = 1; i < (nx - 1); i++)
    {
        for (int j = 1; j < (ny - 1); j++)
        {
            phi[i][j] = 0;
        }
    }
}

void myError( double **u, double **phi_k, int nx, int ny, int &E)
{
    //cout << "here14-1\n";
    double error_treshold = 1E-7;
    //cout << "here14-2\n";
    for (int i=0; i<nx; i++)
    {
        for(int j=0; j<ny; j++)
        {
            if( abs( u[i][j] - phi_k[i][j] ) > error_treshold )
            {
                for (int row=0; row<nx; row++)
                {
                    for (int col=0; col<ny; col++)
                    {
                        phi_k[row][col] = u[row][col];
                    }
                }
                E = 1;
                break;
            }
        }
        if(E) { break; }
    }
    //cout << "here14-3\n";

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
