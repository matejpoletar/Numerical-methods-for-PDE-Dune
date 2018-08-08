//===============================================================
// Driver za mješovite konačne elemente RT tipa.
//===============================================================
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>

#include <dune/istl/bvector.hh>
#include <dune/istl/io.hh>
#include <dune/istl/operators.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/superlu.hh>

#include <dune/pdelab/finiteelementmap/p0fem.hh>
#include <dune/pdelab/finiteelementmap/raviartthomasfem.hh>

#include <dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include <dune/pdelab/gridfunctionspace/interpolate.hh>
#include <dune/pdelab/gridfunctionspace/subspace.hh>

#include <dune/pdelab/constraints/common/constraints.hh>
#include <dune/pdelab/constraints/common/constraintsparameters.hh>
#include <dune/pdelab/constraints/raviartthomas0.hh>

#include <dune/pdelab/gridoperator/gridoperator.hh>

#include <dune/pdelab/backend/istl.hh>

#include <dune/pdelab/common/function.hh>
#include <dune/pdelab/common/vtkexport.hh>

#include "operator.hh"
#include "bctype.hh"
#include "exact.hh"


template <typename GV>
void driver(const GV &gv, int problem, std::string filename)
{
  using DF = typename GV::Grid::ctype;
  using R = double;
  const int dim = GV::dimension;
  // P0 konačni elementi za tlak
  using P0FEM = Dune::PDELab::P0LocalFiniteElementMap<DF, R, dim>;
  P0FEM p0fem(Dune::GeometryType(Dune::GeometryType::cube, dim));

  // RT0 konačni elementi za brzinu
  using RT0FEM = Dune::PDELab::RaviartThomasLocalFiniteElementMap<
                                       GV, DF, R, 0, Dune::GeometryType::cube>;
  RT0FEM rt0fem(gv);

  // make a grid function space
  using VBE = Dune::PDELab::istl::VectorBackend<>;
  // prostor za tlak Y
  using P0GFS = Dune::PDELab::GridFunctionSpace<GV, P0FEM, Dune::PDELab::NoConstraints, VBE>;
  P0GFS p0gfs(gv, p0fem);
  // prostor za brzinu W
  using RT0GFS = Dune::PDELab::GridFunctionSpace<GV, RT0FEM, Dune::PDELab::RT0Constraints, VBE>;
  RT0GFS rt0gfs(gv, rt0fem);
  // Produktni prostor W x Y
  using MGFS = Dune::PDELab::CompositeGridFunctionSpace<VBE, Dune::PDELab::LexicographicOrderingTag, RT0GFS, P0GFS>;
  MGFS mgfs(rt0gfs, p0gfs);

  // Rubni uvjeti. Na tlak nema rubnih uvjeta, pa su ustvari rubni uvjeti dani u BCTypeVelocity.
  // Dirichletov rubni uvjet se uzima varijacijski, a Neumannov se građuje u prostor. O tome brine
  // Dune::PDELab::RT0Constraints klas.
  BCTypePressure conP;
  BCTypeVelocity conV;
  using BCT = Dune::PDELab::CompositeConstraintsParameters<BCTypeVelocity,BCTypePressure>;
  BCT bct(conV, conP);

  // constraints
  typedef typename MGFS::template ConstraintsContainer<R>::Type T;
  T t;                                     // container for transformation
  Dune::PDELab::constraints(bct, mgfs, t, /* verbose = */ true); // fill container

  // Funkcije za interpolaciju rubnog uvjeta
  using VType = VelocityExtension<GV>;
  VType vext(gv);
  typedef Dune::PDELab::PiolaBackwardAdapter<VType> RVType;
  RVType rvext(vext);
  using GType = PressureExtension<GV>;
  GType pext(gv);
  typedef Dune::PDELab::CompositeGridFunction<RVType, GType> UType;
  UType u(rvext, pext);

  // Vektor koeficijenata
  typedef typename Dune::PDELab::Backend::Vector<MGFS, R> X;
  X x(mgfs, 0.0);

  // Interpolacija rubnog uvjeta
  Dune::PDELab::interpolate(u, mgfs, x);
  Dune::PDELab::set_nonconstrained_dofs(t, 0.0, x); // clear interior

  // Lokalni operator
  typedef DiffusionMixed<BCTypeVelocity> LOP;
  LOP lop(conV, 4, 2, problem); // Default 2,1, p=0/1 (verifikacija/problem)

  typedef Dune::PDELab::istl::BCRSMatrixBackend<> MBE;
  MBE mbe(9); // Maximal number of nonzeros per row can be cross-checked by
              // patternStatistics().
  typedef Dune::PDELab::GridOperator<MGFS, MGFS, LOP, MBE, R, R, R, T, T> GO;
  GO go(mgfs, t, mgfs, t, lop, mbe);

  // represent operator as a matrix
  typedef typename GO::Jacobian M;
  M m(go);
  std::cout << m.patternStatistics() << std::endl;
  m = 0.0;
  go.jacobian(x, m);
  //  Dune::printmatrix(std::cout,m.base(),"global stiffness matrix","row",9,1);

  // Sustav je indefiniran pa ga rješavamo direktnom metodom.
  typedef Dune::PDELab::Backend::Native<M> ISTLM;
#if HAVE_SUPERLU
  //  Dune::SuperLU<ISTLM> solver(Dune::PDELab::istl::raw(m), true);
  Dune::SuperLU<ISTLM> solver(Dune::PDELab::Backend::native(m), true);
  Dune::InverseOperatorResult stat;

  X r(mgfs, 0.0);
  go.residual(x, r);
  X z(mgfs, 0.0);
  solver.apply(z, r, stat);
  x -= z;
#else
#error No superLU support, please install and configure it.
#endif

  // Potprostori
  typedef Dune::PDELab::GridFunctionSubSpace<MGFS, Dune::TypeTree::TreePath<0>> VSUB;
  VSUB vsub(mgfs); // Prostor brzine
  typedef Dune::PDELab::GridFunctionSubSpace<MGFS, Dune::TypeTree::TreePath<1>> PSUB;
  PSUB psub(mgfs); // Prostor tlaka

  // make discrete function object
  typedef Dune::PDELab::DiscreteGridFunctionPiola<VSUB, X> RT0DGF;
  RT0DGF rt0dgf(vsub, x);
  typedef Dune::PDELab::DiscreteGridFunction<PSUB, X> P0DGF;
  P0DGF p0dgf(psub, x);


  typedef ExactVelocity<GV> ExactTypeVelocity;
  ExactTypeVelocity exactVel(gv);

  typedef ExactPressure<GV> ExactTypePressure;
  ExactTypePressure exactP(gv);

  // output grid function with VTKWriter
  // Dune::VTKWriter<GV> vtkwriter(gv,Dune::VTK::conforming);
  Dune::SubsamplingVTKWriter<GV> vtkwriter(gv, 1); // plot result
  vtkwriter.addCellData(
      std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<P0DGF>>(
          p0dgf, "pressure"));
  vtkwriter.addVertexData(
      std::make_shared<Dune::PDELab::VTKGridFunctionAdapter<RT0DGF>>(
          rt0dgf, "velocity"));
  if (problem == 0){
  vtkwriter.addVertexData(std::make_unique<Dune::PDELab::VTKGridFunctionAdapter<ExactTypeVelocity>>(exactVel,"exact_velocity"));
   vtkwriter.addCellData(std::make_unique<Dune::PDELab::VTKGridFunctionAdapter<ExactTypePressure>>(exactP,"exact_pressure"));
  }
   vtkwriter.write(filename, Dune::VTK::ascii);

    typedef typename Dune::PDELab::DiscreteGridFunctionPiola<VSUB, X>::Traits::RangeType RF;
   //Izračunavanje integrala
   double value = 0.0;
     double result = 0.0;
     int p = 2;
     for (auto const & element : elements(gv)){

         // Geometry objekt ima sve informacije o geometriji objekta.
         const auto geometry = element.geometry();
         // Uzmimo sada tip geometrije
         const auto gt = geometry.type();
         // Zatražimo kvadraturnu formulu (obavezno koristiti referencu)
         const auto& rule = Dune::QuadratureRules<double,dim>::rule(gt,p);
         result = 0.0;
         for (auto const & qpoint : rule)
         {
         RF fval;

         rt0dgf.evaluate(element,qpoint.position(), fval);
             double weight = qpoint.weight();
                 // | det (grad g) | daje Geometry objekt
             double detjac = geometry.integrationElement(qpoint.position());
             result += fval[0] * weight * detjac;
         }
         value += result;
      }

      std::cout << "Vrijednost traženog integrala iznosi " << value << std::endl;

}
