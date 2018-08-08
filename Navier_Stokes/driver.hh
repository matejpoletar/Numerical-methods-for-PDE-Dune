/* 
 * File:   driver.hh
 * 
 */

#ifndef DRIVER_HH
#define	DRIVER_HH

#include <dune/common/fvector.hh>
#include <dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include <dune/grid/io/file/vtk/vtksequencewriter.hh>

#include <dune/istl/bvector.hh>
#include <dune/istl/operators.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/superlu.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/io.hh>

#include <dune/pdelab/finiteelementmap/pkfem.hh>
#include <dune/pdelab/finiteelementmap/qkfem.hh>
#include <dune/pdelab/constraints/common/constraints.hh>
#include <dune/pdelab/constraints/common/constraintsparameters.hh>
#include <dune/pdelab/constraints/conforming.hh>
#include <dune/pdelab/gridfunctionspace/vectorgridfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/subspace.hh>
#include <dune/pdelab/gridfunctionspace/vtk.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include <dune/pdelab/gridfunctionspace/interpolate.hh>
#include <dune/pdelab/gridoperator/onestep.hh>
#include <dune/pdelab/gridoperator/gridoperator.hh>

#include <dune/pdelab/common/function.hh>
#include <dune/pdelab/common/vtkexport.hh>

#include <dune/pdelab/backend/istl.hh>
//#include <dune/pdelab/backend/istl/bcrsmatrixbackend.hh>
//#include <dune/pdelab/backend/istlmatrixbackend.hh>
//#include <dune/pdelab/backend/istlsolverbackend.hh>
//#include <dune/pdelab/localoperator/cg_stokes.hh>
#include <dune/pdelab/localoperator/taylorhoodnavierstokes.hh>
#include <dune/pdelab/newton/newton.hh>

template<typename GV, typename IF, typename PRM>
void driver(const GV& gv, std::string filename, PRM & parameters, IF & bdry_solution)
{
    static const unsigned int dim = GV::dimension;

    typedef double RF;

    Dune::Timer timer;
    std::cout << "=== Initialize:" << timer.elapsed() << std::endl;
    timer.reset();

    ///////////////////////////////////////////////////////
    // Construct grid function spaces
    // make finite element map
    typedef typename GV::Grid::ctype DF;
    const int k = 2;
    const int q = 2 * k;  // preciznost integracijske formule lokalnog operatora

    // Taylor-Hoodovi elementi -- P2 za brzinu P1 za tlak 
    typedef Dune::PDELab::PkLocalFiniteElementMap<GV, DF, RF, k>      V_FEM;  // komponenta brzine
    typedef Dune::PDELab::PkLocalFiniteElementMap<GV, DF, RF, k - 1 > P_FEM;  // tlak
    V_FEM vFem(gv);
    P_FEM pFem(gv);

    typedef Dune::PDELab::ConformingDirichletConstraints CDC;
    typedef Dune::PDELab::istl::VectorBackend<Dune::PDELab::istl::Blocking::none, 1> VB;
    // Ova klasa direktno konstruira vektorske elemente u R^dim. 
    // Prostor mrežnih funkcija za brzinu (vektorski):  V_h 
    typedef Dune::PDELab::VectorGridFunctionSpace<GV, V_FEM, dim, VB, VB, CDC> PGFS_V_GFS;
    PGFS_V_GFS powerVGfs(gv, vFem);
    powerVGfs.name("velocity");
    // Prostor mrežnih funkcija za tlak (skalarni): W_h
    typedef Dune::PDELab::GridFunctionSpace<GV, P_FEM, CDC, VB> P_GFS;
    P_GFS pGfs(gv, pFem);
    pGfs.name("pressure");
    // Prostor V_h x W_h 
    typedef Dune::PDELab::CompositeGridFunctionSpace
            <VB, Dune::PDELab::LexicographicOrderingTag, PGFS_V_GFS, P_GFS> GFS;
    GFS gfs(powerVGfs, pGfs);
    ///////////////////////////////////////////////////////

    // Primjena Dirichletovih ograničenja
    typedef typename GFS::template ConstraintsContainer<RF>::Type C;
    C cg;
    cg.clear();

    // create Taylor-Hood constraints from boundary-type
    typedef Dune::PDELab::StokesVelocityDirichletConstraints<PRM>                                  ScalarVelocityConstraints;
    typedef Dune::PDELab::PowerConstraintsParameters<ScalarVelocityConstraints, dim>               VelocityConstraints;
    typedef Dune::PDELab::StokesPressureDirichletConstraints<PRM>                                  PressureConstraints;
    typedef Dune::PDELab::CompositeConstraintsParameters<VelocityConstraints, PressureConstraints> Constraints;

    ScalarVelocityConstraints scalarvelocity_constraints(parameters);
    VelocityConstraints       velocity_constraints(scalarvelocity_constraints);
    PressureConstraints       pressure_constraints(parameters);
    Constraints               constraints(velocity_constraints, pressure_constraints);

    Dune::PDELab::constraints(constraints, gfs, cg);

    // Make grid function operator
    typedef Dune::PDELab::TaylorHoodNavierStokes<PRM> LOP;
    LOP lop(parameters, q);

    typedef Dune::PDELab::istl::BCRSMatrixBackend<> MBE;
    MBE mbe(5); // Maximal number of nonzeroes per row can be cross-checked by patternStatistics().

    //Grid Operator
    typedef Dune::PDELab::GridOperator<GFS, GFS, LOP, MBE, RF, RF, RF, C, C> GO;
    GO go(gfs, cg, gfs, cg, lop, mbe);

    // Vektor koeficijenata i interpolacija rubnog uvjeta
    typedef typename GO::Traits::Domain V;
    V x0(gfs);
    x0 = 0.0;
    Dune::PDELab::interpolate(bdry_solution, gfs, x0);
    std::cout << "=== Finished interpolation:" << timer.elapsed() << std::endl;
    timer.reset();

    // Set non constrained dofs to zero
    Dune::PDELab::set_shifted_dofs(cg, 0.0, x0);

    // Linear solver
    typedef Dune::PDELab::ISTLBackend_SEQ_SuperLU LinearSolver;
    LinearSolver ls(false);

    // Solve (possibly) nonlinear problem 
    timer.reset();
    typedef Dune::PDELab::Newton<GO,LinearSolver,V> PDESOLVER;
    PDESOLVER newton(go, x0, ls);
    newton.setReassembleThreshold(0.0);
    newton.setVerbosityLevel(2);
    newton.setMaxIterations(25);
    newton.setLineSearchMaxIterations(30);
    newton.apply();

    Dune::SubsamplingVTKWriter<GV> vtkwriter(gv, 2);
    Dune::PDELab::addSolutionToVTKWriter(vtkwriter, gfs, x0);
    vtkwriter.write(filename, Dune::VTK::ascii);

    // Check residual
    V r(gfs);
    r = 0.;
    go.residual(x0, r);
    std::cout << "Final Residual: " << r.two_norm() << std::endl;

  
}

#endif	/* DRIVER_HH */

