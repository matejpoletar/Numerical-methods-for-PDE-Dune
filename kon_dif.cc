/*
* Implementacija metode konačnih diferencija za Laplaceovu jednadzbu
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

//globalne varijable
double L = 3.14;
double g0 = 0;
double g1 = 0;

//Rjesavamo dif. jednadzbu -u''=sin(x)
double f(double x){
    return std::sin(x);
}

//Egzaktno rjesenje je u(x)=sin(x)
double solution (double x){
    return std::sin(x);
}

template <typename Vec>
void calculateF(Vec &F){
    //racunamo vektor desne strane F
    int n = F.dim();
    double h = L / (n + 1);
    for (unsigned int i = 0; i < n; ++i)
        F[i] = h * h * f((i+1) * h);
    F[0] += g0;
    F[n-1] += g1;
}

template <typename Vec>
void calculateExact(Vec &exact){
    //racunamo vektor rjesenja exact
    int n = exact.dim();
    double h = L / (n + 1);
    for (unsigned int i = 0; i < n; ++i)
        exact[i] = solution((i+1)*h);
}

template <typename Mat>
void calculateMatA (Mat &A){
    //racunamo matricu A
    int n = A.N();
    std::cout << n << std::endl;
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
        A[i][i] = 2.0;
        if (i != n-1) A[i][i+1] = -1.0;
        if (i != 0) A[i][i-1] = -1.0;
    }
}

int main(int argc, char** argv){

    using VBlock = Dune::FieldVector<double, 1>;
    using BlockVector = Dune::BlockVector<VBlock>;
    using MBlock = Dune::FieldMatrix<double, 1, 1>;
    using BlockMatrix = Dune::BCRSMatrix<MBlock>;

    const int n = atoi(argv[1]);

    BlockVector x, f(n), exact(n);

    //Vektor desne strane f
    calculateF(f);

    //Matrica A
    BlockMatrix A;
    A.setSize(n, n);
    A.setBuildMode(BlockMatrix::random);

    calculateMatA(A);

    /*for (unsigned int i = 0; i < n; ++i){
        for (unsigned int j = 0; j < n; ++j)
            std::cout << A[i][j] << " ";
        std::cout << std::endl;
    }*/

    Dune::MatrixAdapter<BlockMatrix,BlockVector,BlockVector> op(A);
    Dune::SeqILUn<BlockMatrix,BlockVector,BlockVector> ilu1(A, 1, 1);
    Dune::CGSolver<BlockVector> bcgs(op, ilu1, 1e-15, 5000, 0);
      // Podaci o rješenju: broj učinjenih iteracija itd.
    Dune::InverseOperatorResult r;

      // inicijalizacija vektora rješenja
      x.resize(f.N(), false);
      x = 0.0;

      BlockVector B = f;
      bcgs.apply(x, f, r);

      if(r.converged){
        std::cout<< "Solver converged.\n";
        std::cout << "No of iterations = " << r.iterations
              << ", reduction = " << r.reduction << std::endl;
      }
      else
          std::cout<< "Solver did not converge.\n";

    std::cout << "Solution vector =\n" <<  g0 << "\n" << x << g1 << "\n" << std::endl;

    calculateExact(exact);
    std::cout << "Exact solution =\n" << g0 << "\n" << exact << g1 << "\n" << std::endl;

    BlockVector error = exact;
    for(unsigned int i=0; i< x.N(); ++i)
        error[i] -= x[i];

    //std::cout << "Error =\n" <<  error << std::endl;
    std::cout << "Greška= " << error.two_norm() << std::endl;

    return 0;


}
