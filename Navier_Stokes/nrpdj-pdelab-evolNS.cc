#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <dune/grid/io/file/gmshreader.hh>
#include <dune/alugrid/grid.hh>
#include <dune/grid/uggrid.hh>

#include <dune/common/parametertreeparser.hh>

#include "bc_extension.hh"
#include "driver.hh"

int main(int argc, char **argv) {
/*
#if !HAVE_SUPERLU
  std::cerr << "Error: These examples work only if SuperLU is available."
            << std::endl;
  exit(1);
#endif*/
  Dune::MPIHelper::instance(argc, argv);
  // Initialize Navier Stokes parameter class from file
  Dune::ParameterTree params;

  std::string input_filename("src_dir/nrpdj-pdelab-evolNS.input");

  std::cout << "Reading input file \"" << input_filename << "\"" << std::endl;
  try {
    Dune::ParameterTreeParser::readINITree(input_filename, params);
  } catch (...) {
    std::cerr << "The configuration file \"navierstokes.ini\" "
                 "could not be read. Exiting..."
              << std::endl;
    exit(1);
  }

  // 1. Konstrukcija mreže. Koristimo  UG Grid i L-shape domenu u 2D
  //Nap: Ja sam prebacio na ALUGrid, jer kod mene na radi UGGrid
  const int dim = 2;
  typedef Dune::ALUGrid<dim,dim,Dune::simplex,Dune::conforming> GridType;
  //ZA UGGrid
  // using GridType = Dune::UGGrid<dim>;
  typedef GridType::LeafGridView GV;
  typedef double RF;

  // Pročitaj sve parametre iz ulazne datoteke
  int refinement = params.get<int>("domain.level");
  std::string outf = params.get<std::string>("OutputFile"); // izlazni VTK file
  RF mu = params.get<double>("physics.mu");
  RF rho = params.get<double>("physics.rho");


  // Mreža se čita iz .msh datoteke
  GridType *pgrid = Dune::GmshReader<GridType>::read(
      "src_dir/grids/mesh.msh", true, false);
  pgrid->globalRefine(refinement);

  const GV &gv = pgrid->leafGridView();

  // 2. Konstrukcija klase paramatara
  // Dune::PDELab::NavierStokesDefaultParameters.
  //
  typedef ZeroFunction<GV, RF, 1> BdryPressure;
  typedef Velocity<GV, RF, 2> BdryVelocity;

  typedef Dune::PDELab::CompositeGridFunction<BdryVelocity, BdryPressure>
      BdrySolution;

  BdryVelocity bdry_velocity(gv);
  BdryPressure bdry_pressure(gv);
  BdrySolution bdry_solution(bdry_velocity, bdry_pressure);

  BCTypeParam bdry_type;

  typedef ZeroFunction<GV, RF, dim> SourceFunction;
  typedef SourceFunction
      NeumannFlux; // artificijelna definicija. Nemamo Neumannovog rubnog uvjeta
  NeumannFlux neumann_flux(gv); // Ne koristi se ali mora biti prisutan
  SourceFunction source_function(gv);

  const bool navier =
      true; // Treba li asemblirati nelinearni konvektivni član ?
  const bool tensor =
      true; // Treba li raditi sa simetriziranim gradijentom ili ne
            // (utječe na interpretaciju rubnih uvjeta)
  typedef Dune::PDELab::NavierStokesDefaultParameters<
      GV, RF, SourceFunction, BCTypeParam, BdrySolution, NeumannFlux, navier,
      tensor>
      LOPParameters;

  LOPParameters parameters(mu, rho, source_function, bdry_type, bdry_solution,
                           neumann_flux);

  // solve problem
  driver<GV, BdrySolution, LOPParameters>(gv, outf, parameters, bdry_solution);

  return 0;
}
