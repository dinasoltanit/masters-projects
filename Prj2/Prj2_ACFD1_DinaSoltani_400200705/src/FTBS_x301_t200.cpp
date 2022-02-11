# include <cstdlib>
# include <iostream>
# include <fstream>
# include <iomanip>
# include <cmath>
# include <ctime>
# include <cstring>
# include <stdio.h>      /* printf */
# include <math.h>       /* sin */
const double pi = 3.14159265359;
using namespace std;
using std::ofstream;
using std::cerr;
using std::endl;

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

int main ( );

//****************************************************************************

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
  cout << "   We use a method known as FTBS:\n";
  cout << "   FT: Forward Time  : du/dt = (u(t+dt,x)-u(t,x))/dt\n";
  cout << "   BS: Backward Space: du/dx = (u(t,x)-u(t,x-dx))/dx\n";
  cout << "\n";

  timestamp ( );

{
  double a;
  double b;
  double c;
  double L;
  double Endtime;
  //string command_filename = "advection_commands.txt";
  //ofstream command_unit;
  //string data_filename = "advection_data.txt";
  //ofstream data_unit;
  string data_filename = "x301_t200.txt";
  ofstream data_unit;
  //ofstream output_t10_x6001;

  double dt; // time-step
  double dx; // size-step
  double nu = 1.48 * 1E-5; // viscosity / density (SI)
  double stab; // stability checker --> Corant Number
  int i;
  int j;
  //int jm1;
  //int jp1;
  int nt = 200; // time steps
  int nx = 301; // mesh numbers
  //int nt_step;
  //int plotstep;
  //double t;
  double u[nt][nx]; // old velocity
  double unew[nt][nx]; // new velocity
  double x[nx];


  L = 300; // Length (m)
  a = 0.0; // x-start
  b = 300.0; // x-end
  c = 330; // Propagation Velocity
  Endtime = 0.6; // Duration (s)

//****************************************************************************
  //nx = 6001; // nx: Number of nodes
  dx = L / ( double ) ( nx - 1 ); // mesh size
//****************************************************************************
  // a loop to initialize the velocity vector.
  for ( i = 0; i < nx; i++ )
  {
      x[i] = x[i] + dx*(i);
      cout << x[i] << " ";
      cout << '\n';
  }
  cout << "x-direction mesh vector, set.\n";
//****************************************************************************
  // setting boundary conditions:
  for ( i = 0; i < ( ((nx-1)/b) *40 ); i++ )
  {
      u[0][i] = 0;
      unew[0][i] = 0;
      cout << " u[0][" << i <<"] = " << u[0][i] << '\n';
  }
  cout << "B.C. part 1, set.\n";

  for ( i = ( ((nx-1)/b) *40 ); i < ( ((nx-1)/b) *100 ); i++ )
  {
      cout << "x = " << x[i] << '\n';
      u[0][i] =  80 * sin( (pi)*(x[i] - 40) / 60.0 );
      unew[0][i] =  80 * sin( (pi)*(x[i] - 40) / 60.0 );
  }
  cout << "B.C. part 2, set.\n";

  for ( i = (((nx-1)/b) *100 ); i < (nx ); i++ )
  {
      u[0][i] = 0;
      unew[0][i] = 0;
  }
  cout << "B.C. part 3, set.\n";
//****************************************************************************
  //nt = 10; // Time Resolution
  dt = Endtime / ( double ) ( nt ); // time-step
  cout << "Time resolution, set. \n";
//****************************************************************************
  // setting initial condition:
  cout << "setting initial condition ... \n";
  for ( i = 0; i < nt; i++ )
  {
      //cout << i << " ";
      u[i][0]=0;
      unew[i][0] = 0;
      //cout << i << "-->" << u[i][0] << " ";
  }
  cout << "Initial condition, set. \n";
//****************************************************************************
  stab = nu * dt / dx / dx;
  cout << " courant number = " << stab << '\n';
  if ( stab > 0.5 )
  {
    cout << " --> The stability criteria do not match <-- ";
    return 0;
  }
//****************************************************************************
// ***** Core Calculation: start *****
  //data_unit.open("C:\\Users\\Diamonds\\Documents\\CFD1\\x6001_t20.txt");
  data_unit.open( data_filename.c_str ( ) );
  cout << "file opened.";

  if( !data_unit ) { // file couldn't be opened
      cerr << "Error: file could not be opened" << endl;
      exit(1);
   }
  //data_unit.open("x6001_t20.txt", ios::app);
  for ( i = 0; i < nt; i++ )
  {
    //cout << "in calc loop \n";
    data_unit << "zone\n";
    //cout << "zone printed. \n";
    //data_unit << fixed << setprecision(2) << x[0] << '\t' << fixed << setprecision(2) << u[i][0] << '\n';
    //data_unit << '\t' << x[0] << '\t' << i << '\t' << u[i][0] << '\n';
    data_unit << x[0] << '\t' << u[i][0] << '\n';
    //fprintf(fptr, "%f\t%f\n", x[0], u[i][0]);
    // First do this for all x in xn:
    for ( j = 1; j < nx; j++ )
    {
        unew[i+1][j] = u[i][j] - (c * dt / dx) * ( u[i][j] - u[i][j-1] );
        //cout << "u_new = " << unew[i+1][j] << '\n';
        u[i+1][j] = unew[i+1][j];
        //fprintf(fptr, "%f\t%f\n", x[j], u[i][j]);
        //data_unit << fixed << setprecision(2) << x[j] << '\t' << fixed << setprecision(2) << u[i][j] << '\n';
        //data_unit << '\t' << x[j] << '\t' << i << '\t' << u[i][j] << '\n';
        data_unit << x[j] << '\t' << u[i][j] << '\n';
        //cout << "printed: x = " << x[j] << '\t' << "u = " << u[i][j] << '\n';
    }
    data_unit << "\n";
  }
  cout << "Calculation, done.\n";
  cout << " courant number = " << stab << '\n';
  data_unit.close();

// ***** Core Calculation: end   *****

//  Free memory.
//
//
//  Terminate.
//
  cout << "\n";
  cout << "  Finite Difference: FTBS for 1D Wave Equation\n";
  cout << "  Array-based Algorithm\n";
  cout << "  Normal end of execution.\n";
  cout << "\n";

  return 0;
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