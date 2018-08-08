/* 
* Računanje H1 norme funkcije
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

using Point = Dune::FieldVector<double, 2>;


double funA(const Point& x) //vraca kvadrat funkcije
{     double f = x[0]*x[0]*x[1]*std::exp(-(x[0]*x[0]+x[1]*x[1]));
      return f * f;
}

double funGradA(const Point& x) //vraca kvadrat apsolutne vrijednosti gradijenta
{
    Point y(0.0);
    y[0] = 2*x[0]*x[1]*(1-x[0]*x[0])*std::exp(-(x[0]*x[0]+x[1]*x[1]));
    y[1] = (1 - 2*x[1]*x[1])*x[0]*x[0]*std::exp(-(x[0]*x[0]+x[1]*x[1]));
    return y[0]*y[0]+y[1]*y[1];
}

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
      result += fval * weight * detjac;
    }

  return result;
}

template<class GV, class F, class G>
double H1norm(GV& gridView, F const & fun, G const & funGrad, int p){
    double value1 = 0.0, value2 = 0.0;
    for (auto const & element : elements(gridView)){
           value1 += integrateEntity(element, fun,  p); //integral funkcije
           value2 += integrateEntity(element, funGrad,  p);//integral gradijenta
    }
    return std::sqrt (value1 + value2);
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

    for (int k = 1; k < 7; k++) {
         double norma = H1norm(gridView, funA, funGradA, k);
         std::cout << "Red int. formule = " << k << "   H1-norma = " << std::scientific <<
                      std::setprecision(12) << norma<< std::endl;
    }
    return 0;
}


