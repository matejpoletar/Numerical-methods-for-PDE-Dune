/*
* Implementacija Newtonove metode za nelinerane jednadzbe
*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <stdlib.h>
#include <iostream>
#include <iomanip>

#include <dune/istl/bvector.hh>
#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/common/parallel/mpihelper.hh>

// Uvedimo globalno tipove s kojima ćemo raditi.
const int dim = 2;
using Vector = Dune::FieldVector<double, dim>;
using Matrix = Dune::FieldMatrix<double, dim, dim>;

/* Rjesavamo nelinerani sustav
 *    sin(x)-y=0
 *    exp(y)-x=0
*/

// dimenzija sustava F(x) = 0.
// Funkcija koja u točki x računa rezidual res (funkciju F(x))
void test_residual(Vector const & x, Vector & res){
    res[0] = std::sin (x[0]) - x[1];
    res[1] = std::exp (x[1]) - x[0];
}

// Funkcija koja u točki x računa jakobijan jac (funkciju grad F(x))
void test_jacobian(Vector const & x, Matrix & jac){
    jac[0][0] = std::cos (x[0]);
    jac[0][1] = -1;
    jac[1][0] = -1;
    jac[1][1] = std::exp(x[1]);
}

/*
//2. testni primjer dim =3
void test_residual(Vector const & x, Vector & res){
    res[0] = std::sin (x[0]) - x[1];
    res[1] = std::exp (x[1]) - x[2];
    res[2] = std::exp (x[2]) - 2;
}

// Funkcija koja u točki x računa jakobijan jac (funkciju grad F(x))
void test_jacobian(Vector const & x, Matrix & jac){
    jac[0][0] = std::cos (x[0]);
    jac[0][1] = -1;
    jac[0][2] = 0;
    jac[1][0] = 0;
    jac[1][1] = std::exp(x[1]);
    jac[1][2] = -1;
    jac[2][0] = 0;
    jac[2][1] = 0;
    jac[2][2] = std::exp(x[2]);
}
*/


//Parametrizirana klasa Newton koja implementira Newtonovu metodu.
//Parametri su tip matrice i tip vektora koji moraju biti kompatibilni.
template <typename Matrix, typename Vector>
class Newton{
 public:
    // Tipovi pokazivača na funkcije.
    using Residual = void (*)(Vector const &, Vector &);
    using Jacobian = void (*)(Vector const &, Matrix &);
    // Konstruktori.
    Newton() : residual(nullptr), jacobian(nullptr), EPS(1E-8), no_iter(-1), converged(false) {}
    void set_residual(Residual res_) { residual = res_; }
    void set_jacobian(Jacobian jac_) { jacobian = jac_; }
    Vector & get_solution() { return sol; }

    void solve(Vector const & x0){
        Vector x = x0;

        double res_Uk; //norma reziduala u trenutnoj iteraciji u(k)
		Vector iter; //iduca iteracija, tj u(k+1)

        residual(x, res);
        jacobian(x, jac);

        while(no_iter < 50){
            no_iter++;
			res_Uk = res.two_norm();			
            jac.solve(dir, res);
            iter = x - dir;
            residual(iter, res);
			
            int k;
            for (k = 0; k < 10; ++k){
				if (res.two_norm() > 0.95*res_Uk){
                    //raspolavljamo, ali ne vise od 10 puta
                    dir /=2;
                    iter = x - dir;
					residual(iter, res);
				}
				else
                    break; //kad je uvjet zadovoljen prihvatimo tu iteraciju i izađemo iz petlje
            }

            if (k == 10)
                break;

            if (res.two_norm() < EPS) {
                converged = true;
				break;
				}

            x = iter;
            jacobian(x,jac);
        }
        sol = iter;
    }

    bool is_converged() const { return converged; }
    int no_iteration() const { return no_iter; }

private:
    Residual residual; // Pokazivač na funkciju koja računa rezidual
    Jacobian jacobian; // Pokazivač na funkciju koja računa jakobijan
    double EPS;        // Tolerancija
    Matrix jac;        // Jakobijan
    Vector sol;        // rješenje u
    Vector dir;        // smjer traženja
    Vector res;        // rezidual F(u)
    int no_iter;       // broj iteracija
    bool converged;    // zastavica konvergirao/nije konvergirao
};

int main(int argc, char **argv)
{
    Dune::MPIHelper::instance(argc, argv);

    Newton<Matrix, Vector> method;
    method.set_jacobian(test_jacobian);
    method.set_residual(test_residual);

    Vector x0;
    x0[0]=x0[1] = 1.5;  // Neka početna iteracija
    method.solve(x0);
    if(method.is_converged()){
        std::cout << " Converged.\n";
        std::cout << " No iter = " <<  method.no_iteration() << "\n";
        std::cout << " Solution = " << method.get_solution() << std::endl;
    }
    else{
        std::cout << " Not converged.\n";
        std::cout << " Last iteration = " << method.get_solution() << std::endl;
    }

    return 0;
}


