#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <string.h>
#include <numeric>
#include "ConfigFile.tpp"
#include <algorithm>

using namespace std;


double energy(const std::vector<double>& fnow, double dx) {
  double ener = 0.0; // TODO: compute quantity E 
    return 0.0;
}


void boundary_condition(vector<double> &fnext, vector<double> &fnow, double const& A, double om, \
		double const& t,double const& dt, \
		vector<double> &beta2, string &bc_l, string &bc_r, int &N)
{
      if (bc_l == "fixe"){
        fnext[0] = 0.0; 
	// NB: on peut aussi utiliser la condition "excitation" et poser A=0
      }else if(bc_l == "libre"){
        fnext[0] = 0.0; // TODO : Modifier pour imposer la condition au bord gauche libre
      }else if (bc_l =="sortie"){
        fnext[0] = 0.0; // TODO : Modifier pour imposer la condition au bord gauche "sortie de l'onde"
      }else if (bc_l == "excitation"){
        fnext[0] = 0.0; // TODO : Modifier pour imposer la condition au bord gauche sinusoidale
      }else{
        cerr << "Merci de choisir une condition aux bord gauche valide" << endl;
      }
	      
      if (bc_r == "fixe"){
        fnext[N-1] = 0.0; 
	// NB: on peut aussi utiliser la condition "excitation" et poser A=0	
      }else if(bc_r == "libre"){
        fnext[N-1] = 0.0; // TODO : Modifier pour imposer la condition au bord droit libre
      }else if (bc_r =="sortie"){
        fnext[N-1] = 0.0; // TODO : Modifier pour imposer la condition au bord droit "sortie de l'onde"
      }else if (bc_l == "excitation"){ 
        fnext[N-1] = 0.0; // TODO : Modifier pour imposer la condition au bord droit sinusoidale
      }else{
        cerr << "Merci de choisir une condition aux bord droit valide" << endl;
      }
}

double finit(double x, double n_init, double L, double f_hat, double x1, double x2, string initialization)
{
  double finit_(0.);
  const double PI = 3.1415926535897932384626433832795028841971e0;

if(initialization=="mode"){
  // TODO: initialiser la fonction f(x,t=0) selon un mode propre
  finit_ = 0.0;
}
else{
  // TODO: initialiser la fonction f(x,t=0) selon la donnée du problème
  finit_ =0.0;
}
  return finit_;
}

//
// Surcharge de l'operateur pour ecrire les elements d'un tableau
//
template <class T> ostream& operator<< (ostream& o, vector<T> const& v)
{
  unsigned int len(v.size());
  for(unsigned int i(0); i < (len - 1); ++i)
    o << v[i] << " ";
  if(len > 0)
    o << v[len-1];
  return o;
}

//
// Main
//
int main(int argc, char* argv[])
{
  const double PI = 3.1415926535897932384626433832795028841971e0;
  const double g  = 9.81;
  double dx;
  double dt;
  double t;
  double Nsteps;
  int stride(0);

  string inputPath("configuration.in"); // Fichier d'input par defaut
  if(argc>1) // Fichier d'input specifie par l'utilisateur ("./Exercice7 config_perso.in")
    inputPath = argv[1];

  ConfigFile configFile(inputPath); // Les parametres sont lus et stockes dans une "map" de strings.

  for(int i(2); i<argc; ++i) // Input complementaires ("./Exercice7 config_perso.in input_scan=[valeur]")
    configFile.process(argv[i]);

  // Parametres de simulation :
  double tfin    = configFile.get<double>("tfin");
  int nx         = configFile.get<int>("nx"); // nb intervalles
  double CFL     = configFile.get<double>("CFL");
  double nsteps  = configFile.get<double>("nsteps");
  double A       = configFile.get<double>("A");
  double f_hat   = configFile.get<double>("f_hat");
  double n_init  = configFile.get<double>("n_init");
  double hL      = configFile.get<double>("hL");
  double hR      = configFile.get<double>("hR");
  double h00     = configFile.get<double>("h00"); // profondeur, cas uniforme
  double x1      = configFile.get<double>("x1");
  double x2      = configFile.get<double>("x2");
  double xa      = configFile.get<double>("xa");
  double xb      = configFile.get<double>("xb");
  double L       = configFile.get<double>("L");
  double om      = configFile.get<double>("om");
  int n_stride(configFile.get<int>("n_stride"));

  int N = nx+1;                                // nb pts de maillage

// Conditions aux bords:
  string bc_l           = configFile.get<string>("cb_gauche");
  string bc_r           = configFile.get<string>("cb_droite");

// Type de forme initiale de la vague: selon donnée Eq.(4) ou mode propre
// ('mode' pour mode propre, autrement Eq.(4))
  string initialization = configFile.get<string>("initialization"); 

// Onde partant vers la gauche ou vers la droite ou statique
// (par exemple 'left', 'right', 'static')
  string initial_state = configFile.get<string>("initial_state");

// Selecteur pour le cas h0 uniforme:
  bool v_uniform        = configFile.get<bool>("v_uniform");

// Selecteur pour choisir le pas de temps:
// true --> dt=tfin/nsteps; t final est exactement tfin
// false --> dt tel que beta_CFL=1; attention, t final n'est pas exactement tfin
  bool impose_nsteps    = configFile.get<bool>("impose_nsteps");

  vector<double> h0(N) ;
  vector<double> vel2(N) ;
  vector<double> x(N) ;
  vector<double> fpast(N), fnow(N), fnext(N), beta2(N);

  dx = L / (N-1);
  bool ecrire_f = configFile.get<bool>("ecrire_f"); // Exporter f(x,t) ou non
 // Eq.(1) ou Eq.(2) [ou Eq.(6) (faculattif)]: Eq1, Eq2 ou Eq6
  string equation_type = configFile.get<string>("equation_type");
  

  for(int i(0); i<N; ++i){ 
     x[i] = i * dx ;
     h0[i] = 0.0;
     if(v_uniform){
        h0[i]  = h00;
     } 
     else {
       // TODO: programmer la fonction h(x) selon la donnée
       h0[i]  = 999.999; // MODIFIER
     }
     vel2[i]  = g* h0[i];
  }
  // maiximal value of u^2 (to be used to set dt)
  auto max_vel2 = std::max_element(vel2.begin(), vel2.end());
  // TODO: set dt for given CFL
  dt = 1.0; // MODIFY
  // TODO: define dt and CFL with given nsteps
  if(impose_nsteps){
    dt  = 1.0; // MODIFY
    CFL = 1.0; // MODIFY
  }

  // Fichiers de sortie :
  string output = configFile.get<string>("output");

  ofstream fichier_x((output + "_x").c_str());
  fichier_x.precision(15);

  ofstream fichier_v((output + "_v").c_str());
  fichier_v.precision(15);

  ofstream fichier_f((output + "_f").c_str());
  fichier_f.precision(15);

  ofstream fichier_en((output + "_en").c_str());
  fichier_en.precision(15);

  // Initialisation des tableaux du schema numerique :

  //TODO initialize f and beta
  for(int i(0); i<N; ++i)
  {
    fpast[i] = 0.;
    fnow[i]  = 0.;
    beta2[i] = 1.0; // TODO: Modifier pour calculer beta^2 aux points de maillage

    fnow[i]  = finit(x[i], n_init,  L, f_hat, x1, x2, initialization);

    if(initial_state =="static"){
      fpast[i] = 0.0; // TODO: system is at rest for t<=0, 
    }
    else if(initial_state =="right"){ 
      fpast[i] = 0.0; // TODO: propagation to the right
    }
    else if(initial_state =="left"){
      fpast[i] = 0.0; // TODO: propagation to the left
    }
  }


  cout<<"beta2[0] is "<<beta2[0]<<endl;
  cout<<"dt is "<< dt <<endl;


  // Boucle temporelle :
  for(t=0.; t<tfin-.5*dt; t+=dt)
  {
    // Ecriture :
    if(stride%n_stride == 0)
    {
      if(ecrire_f) fichier_f << t << " " << fnow << endl;
      fichier_en << t << " " << energy(fnow,dx) << endl;
     }
    ++stride;

    // Evolution :
    for(int i(1); i<N-1; ++i)
    {
      fnext[i] = 0.0; // TODO : Schémas pour les 3 cas, Equation A ou B ou C
    }

    // Impose boundary conditions
    boundary_condition(fnext, fnow, A, om, t, dt, beta2, bc_l, bc_r, N);

    // Mise a jour et préparer le pas suivant:
    fpast = fnow;
    fnow  = fnext;
  }

  if(ecrire_f) fichier_f << t << " " << fnow << endl;
  fichier_x << x << endl;
  fichier_v << vel2 << endl;
  fichier_en << t << " " << energy(fnow,dx) << endl;

  fichier_f.close();
  fichier_x.close();
  fichier_v.close();
  fichier_en.close();

  return 0;
}
