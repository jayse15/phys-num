#include <iostream>
#include <fstream>
#include <iomanip>
#include "ConfigFile.h" // Il contient les methodes pour lire inputs et ecrire outputs
#include <valarray>
#include <cmath>
#include <numeric>
using namespace std;

class Exercice3
{

private:
  double t, dt, tFin;
  double rS, rJ, a, Om;
  double GM=6.674e-11;
  double mS, mJ, xS, xJ;
  int N_excit, nsteps;
  int sampling;
  int last;
  int  nsel_physics;
  int norder;
  bool adapt;
  double tol= 0e0;
  valarray<double> x0 = std::valarray<double>(0.e0, 4); // Correctly initialized
  valarray<double> x  = std::valarray<double>(0.e0, 4); // Correctly initialized
  ofstream *outputFile;

  void printOut(bool write)
  {
    if((!write && last>=sampling) || (write && last!=1))
    {
      double Energy = compute_energy(x[0],x[1],x[2],x[3]);
      *outputFile << t << " "<< x[0] << " " << x[1] << " "<< x[2] << " " << x[3] << " " \
      << Energy<< " "<< nsteps<< endl; // write output on file
      last = 1;
    }
    else
    {
      last++;
    }
  }

  std::valarray<double> get_f(const std::valarray<double>& x) {
    std::valarray<double> f(0.0, 4);

    if ((nsel_physics!=1) and (nsel_physics!=2)){
        cerr << "No dynamics corresponds to this index" << endl;
        return f;
    }

    double grav_term_j;

    if (nsel_physics==1){
      grav_term_j=0;
    } else{
      grav_term_j = -GM*mJ/pow(dist(x, xJ), 3);
    }
    double grav_term_s = -GM*mS/pow(dist(x, xS), 3);
    double ax = (x[2]-xS)*grav_term_s + (x[2]-xJ)*grav_term_j + 2*Om*f[1] + pow(Om, 2)*f[2]; // accéleration selon x
    double ay = x[3]*(grav_term_s + grav_term_j) - 2*Om*f[0] + pow(Om, 2)*f[3]; // accéleration selon y

    f[2]      = x[0];
    f[3]      = x[1];
    f[0]      = ax;
    f[1]      = ay;

    return f;
  }

// Function to compute distance
double dist(const std::valarray<double>& x, double X){
  return sqrt(pow(x[2]-X, 2) + pow(x[3], 2));
}

// Function to compute mechanical energy per mass
double compute_energy(const std::valarray<double>& x) {
  double jup;

  if (nsel_physics==1){
    jup=0;
  } else{
    jup = GM*mJ/dist(x, xJ);
  }

    return 0.5*(pow(x[0], 2) + pow(x[1], 2)) - GM*mS/dist(x, xS) - jup
                             - 0.5*pow(Om, 2)*(pow(x[2], 2) + pow(x[3], 2));
}
void initial_condition(void){
  if(nsel_physics==1){
    xS = 0;
    Om=0;
  }
  else{
    xS = -a*mJ/(mS+mJ);
    xJ = a*mS/(mS+mJ);
    Om = sqrt(GM*mS/(pow(a, 2)*xJ));
  }
  x0[2] = xS+2*a;
  x0[3] = 0;
  x0[1] -= Om*x0[2];
}

std::valarray<double> RK4_do_onestep(const std::valarray<double>& yold, double dt) {
    std::valarray<double> k1, k2, k3, k4;

    k1 = dt*get_f(yold);
    k2 = dt*get_f(yold+k1/2);
    k3 = dt*get_f(yold+k2/2);
    k4 = dt*get_f(yold+k3);

    return yold + 1/6*(k1 + 2*k2 + 2*k3 + k4);
}


public:
  Exercice3(int argc, char* argv[])
  {
    const double pi=3.1415926535897932384626433832795028841971e0;
    string inputPath("configuration.in"); // Fichier d'input par defaut
    if(argc>1) // Fichier d'input specifie par l'utilisateur ("./Exercice3 config_perso.in")
      inputPath = argv[1];

    ConfigFile configFile(inputPath); // Les parametres sont lus et stockes dans une "map" de strings.

    for(int i(2); i<argc; ++i) // Input complementaires ("./Exercice3 config_perso.in input_scan=[valeur]")
      configFile.process(argv[i]);

    tFin         = configFile.get<double>("tFin");            // t final (overwritten if N_excit >0)
    dt       = configFile.get<double>("dt");            // time step (overwritten if nsteps_per >0)
    mJ         = configFile.get<double>("mJ");              // mass of Jupiter
    mS         = configFile.get<double>("mS");              // mass of the Sun
    rS         = configFile.get<double>("rS");              // Radius of the Sun
    rJ          = configFile.get<double>("rJ");              // Radius of Jupiter
    a            = configFile.get<double>("a");               // demi grand-axe (solei-terre hyp MCU)
    x0[0]        = configFile.get<double>("vx0");             // vitesse initiale x
    x0[1]        = configFile.get<double>("vy0");             // vitesse initiale y
    nsel_physics = configFile.get<int>("nsel_physics");       //1) one body problem around mass#2 or 2) one body in rotating reference frame of {1,2}
    adapt        = configFile.get<bool>("adapt");             //if 1=true -> adaptive dt, if 0=false -> fixed dt
    tol          = configFile.get<double>("tol");             //tolerance of the adaptive scheme
    sampling     = configFile.get<int>("sampling");     // number of time steps between two writings on file
    nsteps       = configFile.get<int>("nsteps");        // number of time step per period
    norder       = configFile.get<int>("norder");        // number of time step per period


    // Ouverture du fichier de sortie
    outputFile = new ofstream(configFile.get<string>("output").c_str());
    outputFile->precision(15);
    dt=tFin/nsteps;
  }

  ~Exercice3()
  {
    outputFile->close();
    delete outputFile;
  };

  void run()
  {
    t = 0.;
    initial_condition();
    x=x0;
    last = 0;
    printOut(true);
    std::valarray<double> y1;
    std::valarray<double> y2;
    double d;

    while (t<tFin){
      if (adapt==false){
        x = RK4_do_onestep(x, dt);
      }
      else{
        y1 = RK4_do_onestep(x, dt);
        y2 = RK4_do_onestep(x, dt/2);
        d = abs(y1-y2).max();

        if (d <= tol){
          dt *= pow(tol/d, 1/(norder+1));
        } else {
          while (d>tol){
            dt *= 0.999*pow(tol/d, 1/(norder+1));
            y1 = RK4_do_onestep(x, dt);
            y2 = RK4_do_onestep(x, dt/2);
            d = abs(y1-y2).max();
          }
        }
        x = y2;
      }
    }



    printOut(true); // ecrire le dernier pas de temps
  };

};

int main(int argc, char* argv[])
{
  Exercice3 engine(argc, argv);
  engine.run();

  return 0;
}
