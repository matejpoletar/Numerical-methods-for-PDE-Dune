#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/grid/yaspgrid.hh>
#include <iostream>
#include <vector>

#include "driver.hh"

int main(int argc, char **argv) {
  Dune::MPIHelper::instance(argc, argv);
  if (argc != 2) {
      printf ("Usage: %s [0/1] (0 - verifikacija 1 - problem)\n", argv[0]);
      return 1;
  }
  int p = std::stoi(argv[1]);
  const int dim = 2;
  Dune::FieldVector<double, dim> L(1.0);
  Dune::array<int, dim> N{100, 100};
  std::bitset<dim> periodic(false);
  int overlap = 0;
  Dune::YaspGrid<dim> grid(L, N, periodic, overlap);
  // grid.globalRefine(1);
  using GV = Dune::YaspGrid<dim>::LeafGridView;
  auto const &gv = grid.leafGridView();
  driver(gv, p,  "Yasp2d_rt0q");

  return 0;
}
