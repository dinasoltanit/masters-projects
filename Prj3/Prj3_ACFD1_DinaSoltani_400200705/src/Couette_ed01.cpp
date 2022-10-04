# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <fstream>
# include <ctime>
# include <cmath>


//****************************************************************************
// Couette FLow
//****************************************************************************


using namespace std;
//****************************************************************************
// Functions: Begin
int main ( );

void timestamp ( );

void u0 ( double t_min, double t_max, double y_min, double y_max, int nt, int ny, double t[], double y[], double **u );

void Couette_Flow (double **u, double **u_old, double dt, double dy, double t_min, double y_min, double t_max, double y_max, double Beta, int nt, int ny, double t[], double y[]);

void MyWrite(double ** phi, ofstream &fout, string output_filename, double dt, double dy, double t_min, double y_min, double t_max, double y_max, int nt, int ny, double t[], double y[]);
// Functions: End
//****************************************************************************

int main ( )
//****************************************************************************
{
  string u_data_filename = "Couette_Exp_t001_y000625.txt";
  ofstream u_data_unit;
  bool header;
  double error_treshold;
  int iterations;
  double delta_v;
  int i,j,row,col;
  double k;
  double t_min = 0.0;
  double t_max = 2.0;
  double y_min = 0.00;
  double y_max = 0.05;
  double dt = 0.001;
  double dy = 0.000625;
  double *t, *y;
  double **u;
  double **u_old;
  double w;
  double omega;
  double nu;
  int nt, ny;
  double Beta;
  timestamp ( );

//****************************************************************************
//  Set X and Y values.

  nt = ((t_max - t_min) / dt) + 1;
  ny = ((y_max - y_min) / dy) + 1;
  cout << "nt =" << nt << '\n';
  cout << "ny =" << ny << '\n';

  t = new double[nt];
  y = new double[ny];
cout << "here1\n";
  for ( i = 0; i < nt; i++ )
  {
    t[i] = ( ( double ) ( nt - i - 1 ) * t_min
           + ( double ) (         i  ) * t_max )
           / ( double ) ( nt     - 1 );
  }
cout << "here2\n";
  for ( i = 0; i < ny; i++ )
  {
    y[i] = ( ( double ) ( ny - i - 1 ) * y_min
           + ( double ) (         i  ) * y_max )
           / ( double ) ( ny     - 1 );
  }
//****************************************************************************
//  New u**, u_old**:

    u = new double *[nt];
    cout << "here3\n";
    for (size_t i = 0; i < nt; i++)
    {
        u[i] = new double[ny];
    }

    u_old = new double *[nt];
    cout << "here4\n";
    for (size_t i = 0; i < nt; i++)
    {
        u_old[i] = new double[ny];
    }

//****************************************************************************
//  Setting boundary conditions:
cout << "here5\n";
    u0 ( t_min, t_max, y_min, y_max, nt, ny, t,y, u  );
cout << "here6\n";
    u0 ( t_min, t_max, y_min, y_max, nt, ny, t,y, u_old  );

//****************************************************************************
//  Solution based on LSOR with omega = 1.2:
  nu = 0.0001307;
  cout << "here7\n";
  Beta = nu * dt / (dy * dy);
  cout << "Beta = " << Beta << '\n';
  cout << "here8\n";
  Couette_Flow(u, u_old, dt, dy, t_min, y_min, t_max, y_max, Beta, nt, ny, t, y);
cout << "here9\n";
//****************************************************************************

//****************************************************************************
// data handelling:
cout << "here10\n";
  MyWrite(u, u_data_unit, u_data_filename, dt, dy, t_min, y_min, t_max, y_max, nt, ny, t, y);

  cout << "number of iterations = " << iterations << '\n';
  cout << "Calculation, done.\n";

  delete [] t;
  delete [] y;
  delete [] u;
  delete [] u_old;
//
//  Terminate.
  cout << "\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";

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
void u0 ( double t_min, double t_max, double y_min, double y_max, int nt, int ny, double t[], double y[], double **u )

{
  int i; //
  double U0 = 1.5;  // m/s
  double H  = 0.05; //m
  double pi = 4.0*atan(1.0);
cout << "\n";
//cout << "y = 0.05 m" << "    ";
  for ( i = 0; i < nt; i++ )
  {
    // at y = 0.05 m for all time steps
    u[i][ny-1] = U0;
    //cout << u[i][ny-1];
  }
//cout << "\n";

//cout << "t = 0.0 s" << "    ";
for ( i = 0; i < (ny-1); i++ )
  {
    // at t = 0.00 s for all the space steps
    //u[0][i] = (U0/H) * y[i];
    u[0][i] = 0.0;
    //cout << u[0][i] << '\n';
  }
  //u[0][0] = 0.0;
  u[0][ny-1] = U0;
  //cout << u[0][ny-1];
  //cout << "\n";

  //cout << "y = 0.00 m" << "    ";
  for ( i = 0; i < nt; i++ )
  {
    // at y = 0.0 m for all time steps
    u[i][0] = 0.0;
    //cout << u[i][0];

  }
//cout << "\n";

  /*
  for ( int row = 0; row < nt; row++ )
  {

      for ( int col = 0; col < ny; col++ )
      {

          cout << "row=" << row << " col=" << col << '\n';
          cout << "t=" << t[row] << "   "
               << "u=" << u[row][col] << "  " << '\n';
      }
      cout << '\n';
  }
  */

cout << "here6-4\n";
  //return;
}
//****************************************************************************
void Couette_Flow (double **u, double **u_old, double dt, double dy, double t_min, double y_min, double t_max, double y_max, double Beta, int nt, int ny, double t[], double y[])
{

  cout << "in couette_flow function \n";
//  while(true)
//        {
            for ( int row = 1; row < nt-1; row++ )
            //for ( int row = 1; row < 2; row++ )
            {
                cout << "new time-step\n";
                for ( int col = 1; col < ny-1; col++ )
                //for ( int col = 1; col < 5; col++ )
                {
                    // Core Calculation
                    cout << "u_old = " << u_old[row][col] << '\n';
                    cout << "u_old -1 = " << u_old[row][col-1] << '\n';
                    cout << "u_old +1 = " << u_old[row][col+1] << '\n';
                    u[row][col] = u_old[row][col] + Beta * ( u_old[row][col+1] - 2.0*(u_old[row][col]) + u_old[row][col-1] );
                    cout << "u = " << u[row][col] << '\n';
                    //u_old[row][col] = u[row][col];
                    u[row][col-1] = 1.5;
                    //delta_v += abs(u[row][col] - u_old[row][col]);
                }

                cout << "here8-2\n";

                for ( int col = 1; col < ny-1; col++ )
                {
                     //Core Calculation
                    u[row][col] = u_old[row][col];
                }
                cout << "here8-3\n";
            }

//            iterations += 1;

//            if (delta_v < error_treshold)
//            {
//                break;
//            }
//            else
//            {
//                delta_v = 0;  // Restart counting delta_v for the next iteration
//            }
//        }
        //cout << "number of iterations = " << iterations << '\n';

}

void MyWrite(double ** phi, ofstream &fout, string output_filename, double dt, double dy, double t_min, double y_min, double t_max, double y_max, int nt, int ny, double t[], double y[])
{
//    double t = t_min;
//    double y = y_min;
    fout.open( output_filename.c_str ( ) );
    cout << "file opened.";

    if( !fout ) { // file couldn't be opened
        cerr << "Error: file could not be opened" << endl;
        exit(1);
     }
    //fout << defaultfloat << "VARIABLE= \"I\", \"J\", \"Phi\"\n" << "ZONE I=\t" << nt << " , J=\t" << ny << " , F=point"<< endl;

    for (int i = 0; i < nt; i++)
    {
        fout << "zone\n";
        //y = y_min;
        for (int j = 0; j < ny; j++)
        {
            //fout << scientific << setprecision(10) << t << "\t" << y << "\t" << phi[i][j] << endl;
            fout << scientific << setprecision(10) << y[j] << "\t" << phi[i][j] << endl;
            //y += dy;
        }
        //t += dt;
    }
  //cout << "number of iterations = " << iterations << '\n';

  fout.close();
  cout << "file closed.\n";
}

//*************************************temp*************************************************************

