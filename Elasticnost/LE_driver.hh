#ifndef  __DRIVER_HH_IS_INCLUDED__
#define  __DRIVER_HH_IS_INCLUDED__

#include <dune/pdelab/finiteelementmap/qkfem.hh>
#include <dune/pdelab/constraints/conforming.hh>
#include <dune/pdelab/constraints/common/constraints.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/tags.hh>
#include <dune/pdelab/gridfunctionspace/subspace.hh>
#include <dune/pdelab/gridfunctionspace/vtk.hh>
#include <dune/pdelab/backend/istl.hh>
#include <dune/pdelab/gridoperator/gridoperator.hh>
#include <dune/pdelab/stationary/linearproblem.hh>

#include<dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include <dune/pdelab/function/minus.hh>

#include <string>
#include <stdexcept>

#include "LE_bctype.hh"
#include "LE_operator.hh"
#include "LE_exact.hh"

template <typename GV>
void driver(GV & gv, double E, double nu, double g_vert, double rho, std::string  name)
{
  const int dim = GV::Grid::dimension;
  const int Qk_order = 1;
  using Coord = typename GV::Grid::ctype;
  // skalarni prostor konačnih elemenata
  using FEM0 = Dune::PDELab::QkLocalFiniteElementMap<GV, Coord, double, Qk_order>;
  //using CON = Dune::PDELab::NoConstraints;
  using CON = Dune::PDELab::ConformingDirichletConstraints;
  using VEB0 = Dune::PDELab::istl::VectorBackend<>;
  using GFS0 = Dune::PDELab::GridFunctionSpace<GV,FEM0,CON,VEB0>;
  // Blok vektori su dimenzije dim (poznato pri kompilaciji)
  using VEB = Dune::PDELab::istl::VectorBackend<Dune::PDELab::istl::Blocking::fixed,dim>;
  // Vektorski grid function space. Blok varijabli odgovara entitetu! (kod nas vrhu elementa)
  using GFS = Dune::PDELab::PowerGridFunctionSpace<GFS0, dim, VEB, Dune::PDELab::EntityBlockedOrderingTag>;
  using CC = typename GFS:: template ConstraintsContainer<double>::Type;
  // vektorski rubni uvjeti
  using U_BCTypeParam = Dune::PDELab::PowerConstraintsParameters< BCTypeParam<GV>, dim >;
  using LOP = ElasticityLocalOperator<BCTypeParam<GV>, FEM0>;
  using MBE = Dune::PDELab::istl::BCRSMatrixBackend<>;
  // konstrukcija vektorskog grid  operatora
  using GO = Dune::PDELab::GridOperator<GFS, GFS, LOP, MBE, double, double, double, CC, CC>;
  // vektor komponenti
  using U = typename GO::Traits::Domain;
  // Interpoliramo rubni uvjet
  using BCE0 = BCExtension<GV, double>;
  // Konstruiraj vektorsku funkciju rubnog uvjeta
  using BCE = Dune::PDELab::PowerGridFunction<BCE0, dim>;
  // Linear solver -- sustav je simetričan, koristimo CG
  using LS = Dune::PDELab::ISTLBackend_SEQ_CG_ILU0;
  using SLP = Dune::PDELab::StationaryLinearProblemSolver<GO, LS, U>;
  // Uzmimo tipove za potprostore: 0 -> potprostor pve komponente, 1 -> druge itd.
  using U0SUB = Dune::PDELab::GridFunctionSubSpace<GFS, Dune::TypeTree::TreePath<0> >;
  using U1SUB = Dune::PDELab::GridFunctionSubSpace<GFS, Dune::TypeTree::TreePath<1> >;
   using U2SUB = Dune::PDELab::GridFunctionSubSpace<GFS, Dune::TypeTree::TreePath<2> >;
  // Napravi mrežnu funciju od svake komponenete rješenja
  using U0_DGF =  Dune::PDELab::DiscreteGridFunction<U0SUB, U> ;
  using U1_DGF =  Dune::PDELab::DiscreteGridFunction<U1SUB, U> ;
 using U2_DGF =  Dune::PDELab::DiscreteGridFunction<U2SUB, U> ;
  FEM0 fem0(gv);
  GFS0 gfs0(gv,fem0);
  GFS  gfs(gfs0);
  // rubni uvjet za komponentu
  BCTypeParam<GV> bc0(gv);
  U_BCTypeParam bc( bc0 );
  // odredi Dirichletovu granicu
  CC cc;
  Dune::PDELab::constraints(bc, gfs, cc);
  // Parametri za lokalni operator
  double mu = E/( 2*(1+nu) );
  double lambda = E*nu/( (1+nu)*(1-2*nu) );
  LOP lop(bc0, mu, lambda, g_vert, rho, true);  // use_exact = true
  MBE mbe(25);
  GO go(gfs, cc, gfs, cc, lop, mbe);
  U u(gfs, 0.0);
  BCE0 bce0(gv);
  BCE bce(bce0);
  Dune::PDELab::interpolate(bce, gfs, u);

  // ILI ako razne komponente rješenja imaju različite Dirichletove vrijednosti 
  // using BCE0 = BCExtension0<GV, double>;
  // using BCE1 = BCExtension0<GV, double>;
  // BCE0 bce0(gv);
  // BCE0 bce1(gv);
  // using BCE = Dune::PDELab::CompositeGridFunction<BCE0, BCE0>;
  // BCE bce(bce0, bce1);
  // Dune::PDELab::interpolate(bce, gfs, u);


  LS ls(5000, true);
  SLP slp(go,ls,u, 1e-8);
  slp.apply();
  std::cout << "Problem solved.\n";
  //  Dune::VTKWriter<GV> vtkwriter(gv,Dune::VTK::conforming);
  Dune::SubsamplingVTKWriter<GV> vtkwriter(gv, 1);
  U0SUB u0sub(gfs);
  U1SUB u1sub(gfs);
  U2SUB u2sub(gfs);

  U0_DGF u0_dgf(u0sub, u);
  U1_DGF u1_dgf(u1sub, u);
  U2_DGF u2_dgf(u2sub, u);

  typedef Dune::PDELab::DiscreteGridFunction<U0SUB, U>  ResultType0;
  typedef Dune::PDELab::DiscreteGridFunction<U1SUB, U>  ResultType1;
  typedef Dune::PDELab::DiscreteGridFunction<U2SUB, U>  ResultType2;

  typedef Exact<GV> ExactType;
  ExactType exact(gv);

  // Ispiši mrežne funkcije
  vtkwriter.addVertexData( std::make_unique<Dune::PDELab::VTKGridFunctionAdapter<U0_DGF>>(u0_dgf,"u_x"));
  vtkwriter.addVertexData( std::make_unique<Dune::PDELab::VTKGridFunctionAdapter<U1_DGF>>(u1_dgf,"u_y"));
  vtkwriter.addVertexData( std::make_unique<Dune::PDELab::VTKGridFunctionAdapter<U2_DGF>>(u2_dgf,"u_z"));
  vtkwriter.addVertexData(std::make_unique<Dune::PDELab::VTKGridFunctionAdapter<ExactType>>(exact,"exact"));
  vtkwriter.write(name, Dune::VTK::ascii);

}
#endif   
