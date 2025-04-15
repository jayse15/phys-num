#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "ConfigFile.tpp"


using namespace std;

const double PI=3.1415926535897932384626433832795028841971e0;
const double e0=8.85418781e-12;

// Résolution d'un système d'équations linéaires par élimination de
// Gauss-Jordan:
template<class T>
vector<T>
solve(const vector<T>& diag,
      const vector<T>& lower,
      const vector<T>& upper,
      const vector<T>& rhs)
{
    vector<T> solution(diag.size());
    vector<T> new_diag(diag);
    vector<T> new_rhs(rhs);

    for (int i = 1; i < diag.size(); ++i) {
        double pivot = lower[i - 1] / new_diag[i - 1];
        new_diag[i] -= pivot * upper[i - 1];
        new_rhs[i] -= pivot * new_rhs[i - 1];
    }

    solution[diag.size() - 1] =
      new_rhs[diag.size() - 1] / new_diag[diag.size() - 1];

    for (int i(diag.size() - 2); i >= 0; --i)
        solution[i] = (new_rhs[i] - upper[i] * solution[i + 1]) / new_diag[i];

    return solution;
}


double epsilon(double r, double r1, double R, double eA, double eB)
{
    if ((0<=r) and (r<r1)){
      return eA;
    } else if ((r1<=r) and (r<=R)){
      return eB;
    }
    return 0;
}

double rho_epsilon(double rh0, double r, double r1, double R, bool unif)
{
    if (unif){
      return 1;
    }
    if ((0<=r) and (r<r1)){
      return rh0*sin(PI*r/r1);
    }
    return 0;
}

int main(int argc, char* argv[])
{

    // Read the default input
    string inputPath = "configuration.in";
    // Optionally override configuration file.
    if (argc > 1)
        inputPath = argv[1];

    ConfigFile configFile(inputPath);
    // Override settings
    for (int i = 2; i < argc; i++)
        configFile.process(argv[i]);

    // Set verbosity level. Set to 0 to reduce printouts in console.
    const int verbose = configFile.get<int>("verbose");
    configFile.setVerbosity(verbose);

    // Read geometrical inputs
    const double R  = configFile.get<double>("R");
    const double r1 = configFile.get<double>("r1");

    const double rho0 = configFile.get<double>("rho0");

    // For the analytical comparison
    const bool uniform_rho_case = configFile.get<bool>("uniform_rho_case");

    // Dielectric relative permittivity
    const double epsilon_a = configFile.get<double>("epsilon_a");
    const double epsilon_b = configFile.get<double>("epsilon_b");

    // Boundary conditions
    const double VR = configFile.get<double>("VR");

    // Discretization
    const int N1 = configFile.get<int>("N1");
    const int N2 = configFile.get<int>("N2");

    // Fichiers de sortie:
    string fichier = configFile.get<string>("output");
    string fichier_phi = fichier+"_phi.out";
    string fichier_E   = fichier+"_E.out";
    string fichier_D   = fichier+"_D.out";
    string fichier_rho = fichier+"_rho.out";

    // Create our finite elements
    const int pointCount = N1 + N2 + 1; // Number of finite elements
    const double h1 = (r1 - 0) / N1;
    const double h2 = (R  - r1) / N2;

    // Position of elements
    vector<double> r(pointCount);

    for (int i=0; i<r.size(); i++){
      if (i<N1){
        r[i] = i*h1;
      } else {
        r[i] = N1*h1 + (i-N1)*h2;
      }
    }

    // Arrays initialization
    vector<double> h(pointCount-1); 	// Distance between grid points
    vector<double> midPoint(pointCount-1);  // Midpoint of each grid element

    for (int i=0; i<h.size(); i++){
      if (i<N1){
        h[i] = h1;
      } else {
        h[i] = h2;
      }
    }

    for (int i=0; i<midPoint.size(); i++){
      if (i<N1){
        midPoint[i] = r[i] + h1/2;
      } else {
        midPoint[i] = r[i] + h2/2;
      }
    }

    // Construct the matrix and right-hand side
    vector<double> diagonal(pointCount, 0.0);  // Diagonal
    vector<double> lower(pointCount - 1, 0.0); // Lower diagonal
    vector<double> upper(pointCount - 1, 0.0); // Upper diagonal
    vector<double> rhs(pointCount, 0.0);       // Right-hand-side

    // Loop over the intervals: add the contributions to matrix and rhs
    for (int k = 0; k < pointCount-1; ++k) {
      double integral = epsilon(midPoint[k], r1, R, epsilon_a, epsilon_b)*midPoint[k]/h[k];
      diagonal[k] += integral;
      lower[k] -= integral;
      upper[k] -= integral;
      diagonal[k+1] += integral;

      rhs[k] += rho_epsilon(rho0, midPoint[k], r1, R, uniform_rho_case)*0.5*midPoint[k]*h[k];
      rhs[k+1] += rho_epsilon(rho0, midPoint[k], r1, R, uniform_rho_case)*0.5*midPoint[k]*h[k];
      }


    // boundary condition at r=R (modify the lines below)
    lower[lower.size() - 1] = 0.0;
    diagonal[diagonal.size() - 1] = 1.0;
    rhs[rhs.size() - 1] = VR;


    // Solve the system of equations
    vector<double> phi = solve(diagonal, lower, upper, rhs);

    // Calculate electric field E and displacement vector D
    vector<double> E(pointCount - 1, 0);
    vector<double> D(pointCount - 1, 0);
    vector<double> rho(pointCount - 2, 0);
    for (int i = 0; i < pointCount - 1; ++i) {
        E[i] = (phi[i]-phi[i+1])/h[i];
        D[i] = E[i]*epsilon(midPoint[i], r1, R, epsilon_a, epsilon_b)*e0;
    }
    // Calculate rho with rho = Div D
    for (int i = 0; i < pointCount - 2; ++i) {
      rho[i] = (midPoint[i+1]*D[i+1] - midPoint[i]*D[i])/(midPoint[i+1]-midPoint[i]) /(r[i+1]*e0);
    }

    // Export data
    {
        // Electric potential phi
        ofstream ofs(fichier_phi);
        ofs.precision(15);

        if (r.size() != phi.size())
            throw std::runtime_error("error when writing potential: r and "
                                     "phi does not have size");

        for (int i = 0; i < phi.size(); ++i) {
            ofs << r[i] << " " << phi[i] << endl;
        }
    }

    {
        // Electric field E
        ofstream ofs(fichier_E);
        ofs.precision(15);

        if (r.size() != (E.size() + 1))
            throw std::runtime_error("error when writing electric field: size of "
                                     "E should be 1 less than r");

        for (int i = 0; i < E.size(); ++i) {
            ofs << midPoint[i] << " " << E[i] << endl;
        }
    }
    {
        // Displacement field D
        ofstream ofs(fichier_D);
        ofs.precision(15);

        if (E.size() != D.size())
            throw std::runtime_error("error when writing displacement field: size of "
                                     "D should be equal to E");

        for (int i = 0; i < D.size(); ++i) {
            ofs << midPoint[i] << " " << D[i] << endl;
        }
    }

    {
      // Charge density rho/e0 calculated with D
      ofstream ofs(fichier_rho);
      ofs.precision(15);

      if (rho.size() != D.size()-1)
          throw std::runtime_error("error when writing charge density : size of "
                                   "rho should be 1 less than D");

      for (int i = 0; i < rho.size(); ++i) {
          ofs << r[i+1] << " " << rho[i] << endl;
      }
    }

    return 0;
}
