#include <iostream>
#include <fstream>
#include <iomanip>
#include "ConfigFile.h" // Il contient les methodes pour lire inputs et ecrire outputs
#include <valarray>
#include <cmath> // Se usi fun
using namespace std;

class Exercice2
{

private:
  double t, dt, tFin;
  double m, g, L;
  double Omega, kappa;
  double theta, thetadot;
  double B0, B1;
  double mu;
  double om_0, om_1, nu, Ig;

  int N_excit, nsteps_per, Nperiod;
  int sampling;
  int last;
  ofstream *outputFile;

  void printOut(bool force)
  {
    if((!force && last>=sampling) || (force && last!=1))
    {
      double emec = 0.5*Ig*pow(thetadot, 2) - mu*B0*cos(theta);
      double demec = Ig*thetadot*acceleration(theta, thetadot, t).sum() + mu*B0*sin(theta)*thetadot;
      double pnc  = -mu*B1*sin(Omega*t)*sin(theta)*thetadot - kappa*pow(thetadot, 2);

      *outputFile << t << " " << theta << " " << thetadot << " " << emec << " " << pnc << " " << demec << endl;
      last = 1;
    }
    else
    {
      last++;
    }
  }

    // for a[0] = function of (theta,t), a[1] = function of (w)
  valarray<double> acceleration(double phi, double w, double t_)
  {
    valarray<double> acc = valarray<double>(2);

    acc[0] = -mu/Ig *(B0 + B1*sin(Omega*t_))*sin(phi); // angular acceleration depending on theta and t only
    acc[1] = -kappa/Ig *w; // angular acceleration depending on w only

    return acc;
  }


    void step()
  {
    valarray<double> a0 = acceleration(theta, thetadot, t);
    theta += thetadot*dt + 0.5*(a0[0]+a0[1])*pow(dt, 2);

    double thetadot_half = thetadot + 0.5*(a0[0] + a0[1])*dt;
    valarray<double> a1 = acceleration(theta, thetadot_half, t+dt); // We can do this because a1[1] doesn't depend on t.

    thetadot += 0.5*(a0[0] + a1[0])*dt + a1[1]*dt;

  }


public:
  Exercice2(int argc, char* argv[])
  {
    const double pi=3.1415926535897932384626433832795028841971e0;
    string inputPath("configuration.in"); // Fichier d'input par defaut
    if(argc>1) // Fichier d'input specifie par l'utilisateur ("./Exercice2 config_perso.in")
      inputPath = argv[1];

    ConfigFile configFile(inputPath); // Les parametres sont lus et stockes dans une "map" de strings.

    for(int i(2); i<argc; ++i) // Input complementaires ("./Exercice2 config_perso.in input_scan=[valeur]")
      configFile.process(argv[i]);

    Omega    = configFile.get<double>("Omega");     // frequency of oscillating magnetic field
    kappa    = configFile.get<double>("kappa");     // coefficient for friction
    m        = configFile.get<double>("m");         // mass
    L        = configFile.get<double>("L");         // length
    B1       = configFile.get<double>("B1");     //  B1 part of magnetic fields (oscillating part amplitude)
    B0       = configFile.get<double>("B0");     //  B0 part of magnetic fields (static part amplitude)
    mu       = configFile.get<double>("mu");     //  magnetic moment
    theta    = configFile.get<double>("theta0");    // initial condition in theta
    thetadot = configFile.get<double>("thetadot0"); // initial condition in thetadot
    sampling = configFile.get<int>("sampling");     // number of time steps between two writings on file
    N_excit  = configFile.get<int>("N_excit");      // number of periods of excitation
    Nperiod  = configFile.get<int>("Nperiod");      // number of periods of oscillation of the eigenmode
    nsteps_per= configFile.get<int>("nsteps");      // number of time step per period

    // Ouverture du fichier de sortie
    outputFile = new ofstream(configFile.get<string>("output").c_str());
    outputFile->precision(15);

    Ig = m*pow(L,2)/12; // Moment d'inertie
    om_0 = sqrt(mu*B0/Ig); // Mode propre
    Omega = 2*om_0; // Mode excité

    if(N_excit>0){
      // simulate N_excit periods of excitation
      tFin = N_excit*2*pi/Omega;
      dt   = 2*pi/(Omega*nsteps_per);
    }
    else{
      // simulate Nperiod periods of the eigenmode
      tFin = Nperiod*2*pi/om_0;
      dt   = 2*pi/(om_0*nsteps_per);
    }


    cout << "final time is "<<"  "<< tFin << endl;

  }


  ~Exercice2()
  {
    outputFile->close();
    delete outputFile;
  };

    void run()
  {
    t = 0.;
    last = 0;
    printOut(true);

    while( t < tFin-0.5*dt )
    {
      step();
      t += dt;
      printOut(false);
    }
    printOut(true);
  };

};

int main(int argc, char* argv[])
{
  Exercice2 engine(argc, argv);
  engine.run();

  return 0;
}
