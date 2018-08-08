/*
* Metoda konačnih diferencija za Laplaceovu jednadžbu
*/

#ifdef HAVE_CONFIG_H
# include "config.h"     
#endif

#include <dune/common/parallel/mpihelper.hh>

#include <dune/grid/yaspgrid.hh>         
#include <dune/grid/io/file/vtk.hh>      

#include <dune/istl/bvector.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/operators.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/solvers.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include <iostream>
#include <cmath>

constexpr int dim = 3;  // dimenzija mreže (mora biti konstantna)

/**
 * Koeficijenti jednadžbe
 */
template <int dim>
double fun_k(Dune::FieldVector<double, dim> const & x){
    return 1.0;
}
/**
  * Egzaktno rješenje - za testiranje. Napraviti vlastiti primjer!
 */
template <int dim>
double fun_exact(Dune::FieldVector<double, dim> const & x){
   //primjer 1:
        return  std::exp (-x[0]*x[0]-x[1]*x[1]-x[2]*x[2]);
   //primjer 2:
        //return std::sin(x[0] + x[1]) + std::cos(x[0]+x[2]);
}

/**
 * Funkcija desne strane, izračunata iz egzaktnog rješenja. 
 */
template <int dim>
double fun_f(Dune::FieldVector<double, dim> const & x){
    //primjer 1:
        return (6 - 4*(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]))*std::exp (-x[0]*x[0]-x[1]*x[1]-x[2]*x[2]);
    //primjer 2:
        //return 2*(std::sin(x[0] + x[1]) + std::cos(x[0]+x[2]));
}

/** 
 * Dirichletov rubni uvjet. Zadan je naprosto egzaktnim rješenjem. 
 */
template <int dim>
double fun_g(Dune::FieldVector<double, dim> const & x){
    return fun_exact(x);
}

/**
 * sparsityPattern[i] = skup svih indeksa susjednih elemenata 
 * Broj susjednih elemenata je maksimalno 2*dim, a manji može 
 * biti kad je element na rubu domene. Nepostojeće indekse
 * stavljamo na -1.
 */
std::vector<std::array<int,2*dim>> sparsityPattern;

// "Skalarni" blok-vektor i blok-matrica.
using Vector = Dune::BlockVector<Dune::FieldVector<double,1>>;
using Matrix = Dune::BCRSMatrix<Dune::FieldMatrix<double, 1, 1>>;
// x = vektor rješenja, F = vektor desne strane,
// exact = vektor egzaktnog rješenja
Vector x,F,exact;
// A je matrica sustava
Matrix A;

/**
 * Funkcija koja na danoj mreži formira globani vektor "exact"
 * koji sadrži vrijednosti egzaktnog rješenja u centrima elemenata.
 */
template <typename GV>
void calculateExact(GV const & gv){
    int n = sparsityPattern.size();
    exact.resize(n);
    auto &set = gv.indexSet();
    int index_elementa;
    for (auto const & element : elements(gv)){
        index_elementa = set.index(element);
        exact[index_elementa] = fun_exact (element.geometry().center());
    }
}

/**
 * Funkcija koja računa L2 normu greške. Točno rješenje "exact" i
 * aproksimativno rješenje "x" su globalne varijable.
 */
template <typename GV>
double calculateL2Error(GV const & gv){
    // INVARIJANTA: Rješenje x je izračunato
    auto &set = gv.indexSet();
    double norma = 0.0;
    for (auto const & element : elements(gv)){
        int i = set.index(element);
        norma += (x[i]-exact[i])*(x[i]-exact[i]);

    }
    return std::sqrt(norma);
  // ...
}

/**
 * Funkcija koja na danoj mreži formira globanu varijablu "sparsityPattern".
 * "sparsityPattern" za svaki element daje indekse njegovih susjeda. Elementi
 * se indeksiraju pomoću indexSet-a.
 */
template <typename GV>
void calculateSparsityPattern(GV const & gv)
{
	int count = 0;
	auto &set = gv.indexSet();

    sparsityPattern.resize(set.size(0));

    for (auto const & element : elements(gv)){

       int i = 0;
       // petlja po svim stranicama elementa
        for (auto const & side : intersections(gv, element))
        {
            if (side.boundary())  // Jesmo li na granici domene?
                sparsityPattern[count][i] = -1;
            else
                sparsityPattern[count][i] = set.index(side.outside());

            i++;

        } // kraj petlje po svim stranicama

        count++;
    } // kraj petlje po svim elementima
}

/**
  * Ispis "sparsityPattern" varijable radi kontrole.
 */
void printSparsityPattern()
{
    for (int i = 0; i < sparsityPattern.size(); ++i){
        std::cout << "Susjedni el. od elementa " << i << " su:" << std::endl;
        for (int j = 0; j < 2*dim; ++j )
            std::cout << sparsityPattern[i][j] << " " ;
        std::cout << std::endl;
    }
}


/** Dimenzioniranje matrice i vektora desne strane.
 *  Računanje profila matrice. Metoda koristi prethodno
 *  izračunati "sparsityPattern".
 */
void setMatrixProfile()
{
    int n = sparsityPattern.size();
    A.setSize(n,n);
    A.setBuildMode(Matrix::random);

    for (int i = 0; i < n; ++i){
        int count = 1;
        for (int j = 0; j < 2*dim; ++j){
            if (sparsityPattern[i][j] >= 0) count++;
        }
        A.setrowsize(i,count); //broj nenegativnih el. + 1
    }
    A.endrowsizes();

    for (int i = 0; i < n; ++i){
        A.addindex(i,i);
        for (int j = 0; j < 2*dim; ++j){
            if (sparsityPattern[i][j] >= 0)
                A.addindex(i, sparsityPattern[i][j]);
        }
    }
    A.endindices();

    F.resize(n);
}

/**
 * Asembliranje matrice i vektora desne strane. Dimenzioniranje
 * vektora i profil matrice su obavljeni u setMatrixProfile().
 */
template <typename GV>
void assemble(GV const & gv){

    auto &set = gv.indexSet();
    int index_elementa;
    double d;

    for (auto const & element : elements(gv)){
        index_elementa = set.index(element);
        F[index_elementa] = fun_f (element.geometry().center());
        // petlja po svim stranicama elementa
         for (auto const & side : intersections(gv, element))
         {
             if (side.boundary()){  // Jesmo li na granici domene?
                   d = element.geometry().volume()/side.geometry().volume();
                   A[index_elementa][index_elementa] += 2*fun_k(side.geometry().center())/(d*d);
                   F[index_elementa] += 2*fun_k(side.geometry().center())*fun_g(side.geometry().center())/(d*d);
             }
             else{
                 auto center1 = element.geometry().center();
                 auto center2 = side.outside().geometry().center();
                 d = 0.0;
                 for (int i = 0; i < dim ; ++i)
                     d += (center1[i]-center2[i])*(center1[i]-center2[i]);
                 d = std::sqrt(d);
                 d = side.geometry().volume()*fun_k(side.geometry().center())/(element.geometry().volume()*d) ;
                 A[index_elementa][index_elementa] += d;
                 A[index_elementa][set.index(side.outside())] -= d;
            }

         }
    }

}


int main(int argc, char * argv[])
{
    Dune::MPIHelper::instance(argc, argv);

    int Nel = 10;
    if(argc > 1) Nel = std::stoi(argv[1]);  // Broj elemenata se može zadati kao
                                            // argument komandne linije.

    using GridType = Dune::YaspGrid<dim>;
    using GridView = GridType::LeafGridView;


    // Konstrukcija YaspGid-a.
    Dune::FieldVector<double, dim> L(1.0);             // Duljina stranice
    std::array<int,dim>            s={Nel,Nel,Nel};        // broj ćelija po stranici
    std::bitset<dim>               periodic(false);    // periodičnost u svim smjerovima
    int overlap = 0;                                   // preklapanje za paralelni grid
    GridType grid(L, s, periodic, overlap); // serijska mreža
    GridView gv = grid.leafGridView();
      
    calculateSparsityPattern(gv);  // pomoćna rutina za računanje profila matrice
    printSparsityPattern();      // za debugging

    setMatrixProfile();            // profil matrice
    assemble(gv);
    // asembliranje matrice i vektora desne strane
//   Priprema linearnog solvera.
    Dune::MatrixAdapter<Matrix,Vector,Vector> op(A);
    // inicijaliziraj ILU(n) prekondicioner, s n=1 i faktorom relaksacije 0.92
    Dune::SeqILUn<Matrix,Vector,Vector> ilu1(A, 1, 0.92);
    // iterativni rješavač BiCGSTAB, redukcija greške = 10^(-15),
    //                               max broj iteracija = 5000, verboznost = 0
    Dune::BiCGSTABSolver<Vector> bcgs(op, ilu1, 1e-15, 5000, 0);
    // Podaci o rješenju: broj učinjenih iteracija itd.
    Dune::InverseOperatorResult r;
    // inicijalizacija vektora rješenja
    x.resize(F.N(), false); // false = bez kopiranja vrijednosti iz F.
    x = 0.0;
    // riješi sustav
    Vector Fsave = F; // F može biti prebrisan pri rješavanju
    bcgs.apply(x, F, r);
    // informacije o konvergenciji
    if(r.converged){
      std::cout<< "Solver converged.\n";
      std::cout << "No of iterations = " << r.iterations
            << ", reduction = " << r.reduction << std::endl;
    }
    else            std::cout<< "Solver did not converge.\n";

    // Iscrtavanje rješenja.
    Dune::VTKWriter<GridView> writer(gv);

    writer.addCellData(x,"approximation");  // aproksimativno rješenje
    // Izračunamo egzaktno rješenje
    calculateExact(gv);
    writer.addCellData(exact, "exact");     // egzaktno rješenje

    writer.write("fvlaplace",Dune::VTK::OutputType::ascii);

    // Računanje greške u L2 normi.
    std::cout << "L2 norma greške = " << calculateL2Error(gv) << std::endl;

    return 0;
}
