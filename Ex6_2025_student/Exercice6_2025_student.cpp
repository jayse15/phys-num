#include "ConfigFile.tpp"
#include <chrono>
#include <cmath>
#include <complex> // Pour les nombres complexes
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

using namespace std;
typedef vector<complex<double>> vec_cmplx;

// Fonction resolvant le systeme d'equations A * solution = rhs
// où A est une matrice tridiagonale
template<class T>
void triangular_solve(vector<T> const& diag,  vector<T> const& lower, vector<T> const& upper,
                 vector<T> const& rhs, vector<T>& solution)
{
    vector<T> new_diag = diag;
    vector<T> new_rhs = rhs;

    // forward elimination
    for (int i(1); i < diag.size(); ++i) {
        T pivot = lower[i - 1] / new_diag[i - 1];
        new_diag[i] -= pivot * upper[i - 1];
        new_rhs[i] -= pivot * new_rhs[i - 1];
    }

    solution.resize(diag.size());

    // solve last equation
    solution[diag.size() - 1] = new_rhs[diag.size() - 1] / new_diag[diag.size() - 1];

    // backward substitution
    for (int i = diag.size() - 2; i >= 0; --i) {
        solution[i] = (new_rhs[i] - upper[i] * solution[i + 1]) / new_diag[i];
    }
}

// Potentiel V(x) :
double V(double om, double x_a, double x_b, double x_L, double x_R, double V0, double x)
{
    const double PI = 3.1415926535897932384626433832795028841971e0;

    if (x>=x_L){
      if (x<x_a){
        return 0.5*pow(om, 2)*pow((x-x_a)/(1.0-x_a/x_L), 2);
      }
      if (x<x_b){
        return V0*pow(sin(PI*(x-x_a)/(x_b-x_a)), 2);
      }
      if (x<=x_R){
        return 0.5*pow(om, 2)*pow((x-x_b)/(1.0-x_b/x_R), 2);
      }
    }
    return 0;
}

// Declaration des diagnostiques de la particule d'apres sa fonction d'onde psi :
//  - prob: calcule la probabilite de trouver la particule dans un intervalle [x_i, x_j]
//  - E:    calcule son energie,
//  - xmoy: calcule sa position moyenne,
//  - x2moy:calcule sa position au carre moyenne,
//  - pmoy: calcule sa quantite de mouvement moyenne,
//  - p2moy:calcule sa quantite de mouvement au carre moyenne.


// Calcule la probabilite de trouver la particule dans un intervalle [x_i, x_j]
double prob(int x_i, int x_j, vec_cmplx const& psi, double hx)
{
    double sum = 0;
    for (int x=x_i; x<x_j; x++){
      sum += norm(psi[x]) + norm(psi[x+1]);
    }
    return 0.5*hx*sum;
}

// Calcule l'energie
double E(vec_cmplx const& dH, vec_cmplx const& aH, vec_cmplx const& cH, vec_cmplx const& psi, double hx){
    complex<double> sum = 0;
    int N = psi.size();
    for (int i = 1; i < N-1; i++) {
        sum += conj(psi[i]) * dH[i] * psi[i];
        sum += std::conj(psi[i]) * aH[i - 1] * psi[i - 1];
        sum += std::conj(psi[i]) * cH[i] * psi[i + 1];
    }

    return hx*sum.real(); // No 0.5 factor because psi is 0 at edges
}

// Calcule x moyen
double xmoy(vector<double> x, vec_cmplx const& psi, double hx)
{
    complex<double> sum = 0;
    for (int i = 1; i < psi.size()-1; i++) {
        sum += conj(psi[i]) * x[i] * psi[i];
      }
    return hx*sum.real();
}

// Calculer x^2 moyen
double x2moy(vector<double> x, vec_cmplx const& psi, double hx)
{
    complex<double> sum = 0;
    for (int i = 1; i < psi.size()-1; i++) {
        sum += conj(psi[i]) * x[i]*x[i] * psi[i];
      }
    return hx*sum.real();
}

// Calcule p moyen
double pmoy(vec_cmplx const& psi, double hx)
{
    int N = psi.size();
    std::complex<double> sum = 0.0;

    for (int i = 1; i < N - 1; ++i) {
        std::complex<double> dpsi  = (psi[i + 1] - psi[i - 1]) / (2.0 * hx);
        sum += std::conj(psi[i]) * (-1i * dpsi);
    }

    return hx * sum.real();
}

// Calcule p^2 moyen
double p2moy(vec_cmplx const& psi, double hx)
{
    int N = psi.size();
    std::complex<double> sum = 0.0;

    for (int i = 1; i < N - 1; ++i) {
        std::complex<double> dpsi2  = (psi[i + 1] - 2.0*psi[i] + psi[i - 1]) / (hx * hx);
        sum += std::conj(psi[i]) * (-dpsi2);
    }

    return hx * sum.real();
}

// Calcule la normalization
vec_cmplx normalize(vec_cmplx const& psi, double const& dx)
{
    vec_cmplx psi_norm(psi.size());
    double integral = 0.0;
    for (int i=0; i<psi.size()-1; i++){
        integral+= dx*(norm(psi[i]) + norm(psi[i+1]))/2.0;
    }
    for (int i=0; i<psi.size(); i++){
        psi_norm[i] = psi[i]/sqrt(integral);
    }
    return psi_norm;
}




int
main(int argc, char** argv)
{
    complex<double> complex_i = complex<double>(0, 1); // Nombre imaginaire i
    const double PI = 3.1415926535897932384626433832795028841971e0;

    string inputPath("configuration.in.example"); // Fichier d'input par defaut
    if (argc > 1) // Fichier d'input specifie par l'utilisateur ("./Exercice8 config_perso.in")
        inputPath = argv[1];

    ConfigFile configFile(
      inputPath); // Les parametres sont lus et stockes dans une "map" de strings.

    for (int i(2); i < argc;
         ++i) // Input complementaires ("./Exercice8 config_perso.in input_scan=[valeur]")
        configFile.process(argv[i]);

    // Set verbosity level. Set to 0 to reduce printouts in console.
    const int verbose = configFile.get<int>("verbose");
    configFile.setVerbosity(verbose);

    // Parametres physiques :
    double hbar = 1.;
    double m = 1.;
    double tfin = configFile.get<double>("tfin");
    double xL = configFile.get<double>("xL");
    double xR = configFile.get<double>("xR");
    double xa = configFile.get<double>("xa");
    double xb = configFile.get<double>("xb");
    double V0 = configFile.get<double>("V0");
    double om0 = configFile.get<double>("om0");
    double n  = configFile.get<int>("n"); // Read mode number as integer, convert to double
    double L = (xR - xL);
    double x0 = configFile.get<double>("x0");
    double sigma0 = configFile.get<double>("sigma_norm") * L;

    int Nsteps = configFile.get<int>("Nsteps");
    int Nintervals = configFile.get<int>("Nintervals");

    // Initialise le paquet d'onde, equation (4.116) du cours

    double k0 = n*2*PI/L;

    int Npoints = Nintervals + 1;
    double dx = L / Nintervals;
    double dt = tfin / Nsteps;

    const auto simulationStart = std::chrono::steady_clock::now();

    // Maillage :
    vector<double> x(Npoints);
    for (int i(0); i < Npoints; ++i)
        x[i] = xL + i * dx;

    // Initialisation de la fonction d'onde :
    vec_cmplx psi(Npoints);

    // initialization time and position to check Probability
    double t = 0;
    unsigned int Nx0 = floor((0 - xL)/(xR-xL)*Npoints);

    // initialize psi
    for (int i(0); i < Npoints; ++i)
    	psi[i] = exp(complex_i*k0*x[i]) * exp(-pow((x[i]-x0), 2)/(2*sigma0*sigma0));

    // Modifications des valeurs aux bords :
    psi[0] = complex<double>(0., 0.);
    psi[Npoints - 1] = complex<double>(0., 0.);

    // Normalisation :
    psi = normalize(psi, dx);

    // Matrices (d: diagonale, a: sous-diagonale, c: sur-diagonale) :
    vec_cmplx dH(Npoints), aH(Nintervals), cH(Nintervals); // matrice Hamiltonienne
    vec_cmplx dA(Npoints), aA(Nintervals),
      cA(Nintervals); // matrice du membre de gauche de l'equation (4.100)
    vec_cmplx dB(Npoints), aB(Nintervals),
      cB(Nintervals); // matrice du membre de droite de l'equation (4.100)

    complex<double> a =
      complex_i * hbar * dt / (4.*m*dx*dx); // Coefficient complexe a

    complex<double> b =
      complex_i * dt / (2.*hbar); // Coefficient complexe b

    // Ces matrices sont stockées sous forme tridiagonale, d:diagonale, c et a: diagonales
    // supérieures et inférieures
    for (int i(0); i < Npoints; ++i) // Boucle sur les points de maillage
    {
        dH[i] = hbar*hbar/(m*(dx*dx)) + V(om0, xa, xb, xL, xR, V0, x[i]);
        dA[i] = 1.0 + 2.*a + b*V(om0, xa, xb, xL, xR, V0, x[i]); 
        dB[i] = 1.0 - 2.*a - b*V(om0, xa, xb, xL, xR, V0, x[i]);
    }
    for (int i(0); i < Nintervals; ++i) // Boucle sur les intervalles
    {
        aH[i] = -hbar*hbar/(2.*m*(dx*dx));
        aA[i] = -a;
        aB[i] = a;
        cH[i] = aH[i];
        cA[i] = aA[i];
        cB[i] = aB[i];
    }

    // Conditions aux limites: psi nulle aux deux bords
    dA[0] = dA[dA.size()-1] = 1;
    dB[0] = dB[dB.size()-1] = 1;

    cA[0] = cA[cA.size()-1] = aA[0] = aA[aA.size()-1] = 0;
    cB[0] = cB[cB.size()-1] = aB[0] = aB[aB.size()-1] = 0;




    // Fichiers de sortie :
    string output = configFile.get<string>("output");

    ofstream fichier_potentiel((output + "_pot.out").c_str());
    fichier_potentiel.precision(15);
    for (int i(0); i < Npoints; ++i)
        fichier_potentiel << x[i] << " " << V(om0, xa, xb, xL, xR, V0, x[i]) << endl;
    fichier_potentiel.close();

    ofstream fichier_psi((output + "_psi2.out").c_str());
    fichier_psi.precision(6);

    ofstream fichier_observables((output + "_obs.out").c_str());
    fichier_observables.precision(15);

    // t0 writing
    for (int i(0); i < Npoints; ++i){
        fichier_psi << abs(psi[i])  << " " << real(psi[i]) << " "  << imag(psi[i]) << " ";
        }
    fichier_psi << endl;

    // Ecriture des observables :
    fichier_observables << t << " " << prob(0, Nx0, psi, dx) << " " << prob(Nx0, Nintervals, psi, dx)
                << " " << E(dH, aH, cH, psi, dx) << " " << xmoy (x, psi, dx) << " "
                << x2moy(x, psi, dx) << " " << pmoy (psi, dx) << " " << p2moy(psi, dx) << endl;

    // Boucle temporelle :
    while (t < tfin) {

        // Multiplication psi_tmp = B * psi :
        vec_cmplx psi_tmp(Npoints, 0.);
        for (int i(0); i < Npoints; ++i)
            psi_tmp[i] = dB[i] * psi[i];
        for (int i(0); i < Nintervals; ++i) {
            psi_tmp[i] += cB[i] * psi[i + 1];
            psi_tmp[i + 1] += aB[i] * psi[i];
        }

        // Resolution de A * psi = psi_tmp :
        triangular_solve(dA, aA, cA, psi_tmp, psi);
        t += dt;

        // Writing
        for (int i(0); i < Npoints; ++i){
            fichier_psi << abs(psi[i])  << " " << real(psi[i]) << " "  << imag(psi[i]) << " ";
            }
        fichier_psi << endl;

        // Ecriture des observables :
        fichier_observables << t << " " << prob(0, Nx0, psi, dx) << " " << prob(Nx0, Nintervals, psi, dx)
                << " " << E(dH, aH, cH, psi, dx) << " " << xmoy (x, psi, dx) << " "
                << x2moy(x, psi, dx) << " " << pmoy (psi, dx) << " " << p2moy(psi, dx) << endl;

    } // Fin de la boucle temporelle





    fichier_observables.close();
    fichier_psi.close();

    const auto simulationEnd = std::chrono::steady_clock::now();
    const std::chrono::duration<double> elapsedSeconds = simulationEnd - simulationStart;
    std::cout << "Simulation finished in " << setprecision(3) << elapsedSeconds.count()
              << " seconds" << std::endl;
}
