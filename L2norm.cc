/* 
* Računanje L2 norme funkcije
*/

#include "config.h"
#include <dune/grid/yaspgrid.hh>
#include "dune/common/parallel/mpihelper.hh" 
#include <exception>
#include <iostream>
#include <iomanip>

#include <dune/geometry/quadraturerules.hh>

//#include "functors.hh"
//#include "integrateentity.hh"

//primjer 1
double funA(const Dune::FieldVector<double,2>& x)
{   double f = std::abs(x[0]*x[0]*x[1]*x[1]+sin(x[0])*sin(x[1]));
    return f;
}

//primjer 2
/*
double funA(const Dune::FieldVector<double,2>& x)
{   return std::exp(-3.234*(x*x)); //f(x)>0,za svaki x pa ne treba abs()

}*/

// Integriraj po elementu sa zadanom točnošću
template<class Entity, class F>
double integrateEntity (const Entity &element, F const & fun, int p)
{
  // dimension of the entity
  const int dim = Entity::dimension;
  // Geometry objekt ima sve informacije o geometriji objekta.
  const auto geometry = element.geometry();
  // Uzmimo sada tip geometrije
  const auto gt = geometry.type();
  // Zatražimo kvadraturnu formulu (obavezno koristiti referencu)
  const auto& rule = Dune::QuadratureRules<double,dim>::rule(gt,p);
  //  provjeri da smo dobili traženi red.
  if (rule.order()<p) std::exit(-1);

  // compute approximate integral
  double result=0;
  for (auto const & qpoint : rule)
    {
      double fval = fun(geometry.global(qpoint.position()));
      double weight = qpoint.weight();
          // | det (grad g) | daje Geometry objekt
      double detjac = geometry.integrationElement(qpoint.position());
      result += fval* fval * weight * detjac; //integriramo kvadrat funkcije
    }

  return result;
}

template<class GV, class F>
double uniformintegration(GV& gridView, F const & fun,  int p) {
     double value = 0.0;
     for (auto const & element : elements(gridView))
            value += integrateEntity(element, fun,  p);
     return value;
}

template<class GV, class F>
double L2norm(GV& gridView, F const & fun){
    double value = 0.0;
    for (auto const & element : elements(gridView))
           value += integrateEntity(element, fun,  3);
    return std::sqrt (value);
}

int main(int argc, char **argv)
{
    Dune::MPIHelper::instance(argc, argv);

    const int dim = 2;
    typedef Dune::YaspGrid<dim> GridType;
    Dune::FieldVector<GridType::ctype, dim> L(1.0);
    std::array<int, dim> N {5,5};
    std::bitset<dim>    periodic(false);
    int overlap = 0;
    GridType grid(L, N, periodic, overlap);
    auto gridView = grid.leafGridView();

    // Integriramo 6 puta i svaki puta profinjujemo mrežu
    for (int k = 0; k < 6; k++) {
         double norma = L2norm(gridView, funA);
         std::cout << "br elemenata=" << std::setw(8) << std::right << gridView.size(0)
                   << " L2-norma=" << std::scientific << std::setprecision(12)
                   << norma<< std::endl;

         grid.globalRefine(1);
    }
    return 0;
}


