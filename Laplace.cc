
/*
*Rjesavanje PDJ -Laplace(u)=f na 2D domeni [0,1]x[0,1]
*/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <cassert>
#include <iostream>
#include <iomanip>

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>

#include <dune/istl/bvector.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/operators.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/preconditioners.hh>

#include <dune/common/parallel/mpihelper.hh>



int const dim = 2;

//Egzaktno rjesenje
template <int dim>
double fun_exact(Dune::FieldVector<double, dim> const & x){
    return std::exp (-x[0]*x[0]-x[1]*x[1]);
}


//Funkcija desne strane, izračunata iz egzaktnog rješenja.
template <int dim>
double fun_f(Dune::FieldVector<double, dim> const & x){
    return 4*(1-x[0]*x[0]-x[1]*x[1])* std::exp (-x[0]*x[0]-x[1]*x[1]);
}

//Dirichletov rubni uvjet
template <int dim>
double fun_g(Dune::FieldVector<double, dim> const & x){
    return fun_exact(x);
}

template <typename F>
void calculateExact(F & exact, double hx, double hy){
   Dune::FieldVector<double, dim> x;

   for (unsigned int i = 0; i < exact.N(); ++i){		//petlja ide po blokovima
       x[1] = (i+1)*hy; //druga koordinata tj y je konstanta za svaki blok, pa ju jednom racunamo za svaki blok
       for (unsigned int j = 0; j < exact[0].dim(); ++j){ //petlja ide po elementima u bloku
           x[0]= (j+1)*hx;
           exact[i][j] = fun_exact(x); //j-ta komponeneta u i-tom bloku
       }
   }


}

//Blokovi za izgradnju matrice A - dva su tipa blokova
template <typename Block>
void Block1 (Block & B1, double a, double b){
    int n = B1.N();
	B1[0][0] = a;
	B1[0][1] = b;
    for (int i = 1; i < n-1; ++i){
        B1[i][i-1] = b;
		B1[i][i] = a;
        B1[i][i+1] = b;
    }
	B1[n-1][n-2] = b;
	B1[n-1][n-1] = a;
}

template <typename Block>
void Block2 (Block & B2, double c){
    for (int i = 0; i < B2.N(); ++i)
        B2[i][i] = c;
}

//Asembliranje matrice sustava A
template <typename F, typename Block>
void assembleA(F & A,Block &B1, Block &B2){

    int n = A.N();
    for (unsigned int i = 0; i < n; ++i){
        if (i == 0 || i == n-1)
            A.setrowsize(i,2);
        else
            A.setrowsize(i,3);
    }
    A.endrowsizes();

    for (unsigned int i = 0; i < n; ++i){
        A.addindex(i,i);
        if (i != n-1) A.addindex(i,i+1);
        if (i != 0) A.addindex(i,i-1);
    }
    A.endindices();

    for (unsigned int i = 0; i < n; ++i){
       A[i][i] = B1;
       if (i != n-1) A[i][i+1] = B2;
       if (i != 0) A[i][i-1] = B2;
   }

}

//Asembliranje vektora desne strane F
template <typename G>
void assembleF(G & F, double hx, double hy, double b, double c){

	int n = F.N(); //broj blokova
	int m = F[0].dim(); //dimenzija blokova
	
    Dune::FieldVector<double, dim> x;
	
    for (int i = 0; i < n; ++i){		//petlja ide po blokovima
		x[1] = (i+1)*hy; //druga koordinata tj y je konstanta za svaki blok, pa ju jednom racunamo za svaki blok
        for (int j = 0; j < m; ++j){ //petlja ide po elementima u bloku
			x[0]= (j+1)*hx; 
			F[i][j] = fun_f(x); //j-ta komponeneta u i-tom bloku
		}
        //Rubni uvjeti-prva i zadnja komponenta u svakom bloku:
		x[0]=0;
		F[i][0] -= b* fun_g(x);
		x[0]=1;
		F[i][m-1] -= b* fun_g(x);

	}
    //Treba jos korigirati prvi i zadnji blok (sve komponente u oba bloka):
	for (unsigned int j = 0; j < m; ++j ){
		x[0] = (j+1)*hx;
		x[1] = 0; 
		F[0][j] -= c*fun_g(x); //prvi blok
		x[1] = 1;
		F[n-1][j] -= c*fun_g(x); //zadnji blok
	}
	
}

//L2 norma greske
template <typename F>
double calculateL2Error (F &x, F &exact){
    int n = x.N();
    int m = x[0].dim();
    assert (n==exact.N());
    assert (m==exact[0].dim());

    double norma = 0;

    for (int i = 0; i < n; ++i)
        for (int j = 0; j < m; ++j)
            norma += (std::abs(x[i][j]-exact[i][j])) * (std::abs(x[i][j]-exact[i][j]));

    return std::sqrt(norma);
    }


int main(int argc, char * argv[])
{
    Dune::MPIHelper::instance(argc, argv);

    const int Nx = 10; //broj elemenata u x smjeru
    const int Ny = 10; //broj elemenata u y smjeru

    using VBlock = Dune::FieldVector<double, Nx-1>;
    using BlockVector = Dune::BlockVector<VBlock>;
    using MBlock = Dune::FieldMatrix<double, Nx-1, Nx-1>;
    using BlockMatrix = Dune::BCRSMatrix<MBlock>;

    double a = 2 *( Nx * Nx + Ny * Ny );
    double b = - Nx * Nx;
    double c = - Ny * Ny;

    double hx = 1.0 / Nx;
    double hy = 1.0 / Ny;

    // x = vektor rješenja, exact = vektor egzaktnog rješenja
    BlockVector x(Ny-1), exact(Ny-1);

    // F = vektor desne strane
    
    BlockVector F(Ny-1); //Ny-1 broj blokova
    assembleF(F, hx, hy, b, c);

    // A je matrica sustava
    BlockMatrix A(Ny-1,Ny-1,BlockMatrix::random);
    //asembliranje matrice A
    MBlock B1, B2;
    Block1(B1,a,b);
    Block2(B2, c);
    assembleA(A,B1,B2);

/*ispis matrice za provjeri (samo ne-nul blokovi)
    for (int i = 0; i < A.N(); ++i){
        for (int j = 0;j < A.N(); ++j)
            if (j == i+1 || j==i || j==i-1)std::cout << A[i][j] << std::endl;
        std::cout << std::endl;
    }

    std::cout << "F = \n" << F <<std::endl;
    */

//  Priprema linearnog solvera.
    Dune::MatrixAdapter<BlockMatrix,BlockVector,BlockVector> op(A);
    // inicijaliziraj ILU(n) prekondicioner, s n=1 i faktorom relaksacije 0.92
    Dune::SeqILUn<BlockMatrix,BlockVector,BlockVector> ilu1(A, 1, 0.92);
    //rješavač, red.greške=1e-8, verboznost=5000
    Dune::CGSolver<BlockVector> cgs(op, ilu1, 1e-8, 5000, 0);
    // Podaci o rješenju: broj učinjenih iteracija itd.
    Dune::InverseOperatorResult r;
    // inicijalizacija vektora rješenja
    x.resize(F.N(), false); // false = bez kopiranja vrijednosti iz F.
    x = 0.0;
    // riješi sustav
    BlockVector FF = F; // F može biti prebrisan pri rješavanju
    cgs.apply(x, F, r);

    // Izračunamo egzaktno rješenje
    calculateExact(exact, hx, hy);

    // informacije o konvergenciji
    if(r.converged){
      std::cout<< "Solver converged.\n";
      std::cout << "Solution=\n" << x << std::endl;
      std::cout << "No of iterations = " << r.iterations
            << ", reduction = " << r.reduction << std::endl;
      // Računanje greške u L2 normi.
      std::cout << "Exact solution=\n" << exact << std::endl;
      std::cout << "L2 norma greške = " << calculateL2Error(x, exact) << std::endl;
    }
    else
        std::cout<< "Solver did not converge.\n";

    return 0;
}
