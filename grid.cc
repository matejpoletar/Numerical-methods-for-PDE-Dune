//Kretanje po mrezi - Zadatak 4
//Određivanje najmanjeg i najvećeg kuta u mreži

#ifdef HAVE_CONFIG_H
# include "config.h"     // Datoteka koju kreira configure skripta i koja služi
#endif                   // adaptaciji koda na okolinu u kojoj se kompilira.

#include <iostream>     // Ulazno-izlazna biblioteka
#include <cstdlib>
#include <cmath>
#include "dune/common/parallel/mpihelper.hh"  // Inicijalizacija MPI sustava (za paralelno izvršavanje programa)

#include <dune/alugrid/grid.hh>           // koristimo ALUGrid
#include <dune/grid/io/file/gmshreader.hh> // GmshReader klasa
#include <dune/grid/common/gridinfo.hh>    // Informacije o mreži.
#include <dune/grid/io/file/vtk.hh>

#define PI 3.14159265

// glavni program
// Ime datoteke tipa .msh očekujemo kao prvi argument komandne linije.
int main(int argc, char** argv)
{
    Dune::MPIHelper::instance(argc, argv);

    if(argc < 3){
       std::cerr << "Usage: " << argv[0] << " ime_grid_datoteke.msh br_profinjenja" << std::endl;
       std::exit(1);
    }

    const int dim = 2 ;
    using GridType =  Dune::ALUGrid<dim,dim,Dune::simplex,Dune::conforming>;
    using GridView = GridType::LeafGridView;
    // statička metoda  Dune::GmshReader<GridType>::read vraća pokazivač na kreirani grid.
    GridType* pgrid = Dune::GmshReader<GridType>::read(argv[1]);
    int no_r = std::stoi(argv[2]);
    pgrid->globalRefine(no_r);     // profini mrežu

    auto gridView = pgrid->leafGridView();

    // Dimenzija mreže
    std::cout << "Dimenzija mreze = " << dim << std::endl;
    // Tip globalne koordinate
    using GlobalCoordinate = GridView::template Codim<0>::Geometry::GlobalCoordinate;

    double amin = 180, amax = 0;
    for(auto const & element : elements(gridView))
    {
        // Geometrijski tip elementa (tip Dune::GeometryType)
        auto gt = element.type();
        auto geom = element.geometry();
        // Broj vrhova elementa
        int n_v = geom.corners();

        // Uzmi koordinate svih vrhova elemenata
        std::vector<decltype(geom.corner(0))> coo_v(n_v);
        for(unsigned int i=0; i < n_v; ++i) coo_v[i] = geom.corner(i);

	//koordinate vrhova trokuta
        auto A = coo_v[0];
        auto B = coo_v[1];
        auto C = coo_v[2];

	//duljine stranica u trokutu
        double a = sqrt ((B[0]-C[0])*(B[0]-C[0])+(B[1]-C[1])*(B[1]-C[1]));
        double b = sqrt ((A[0]-C[0])*(A[0]-C[0])+(A[1]-C[1])*(A[1]-C[1]));
        double c = sqrt ((A[0]-B[0])*(A[0]-B[0])+(A[1]-B[1])*(A[1]-B[1]));

	//kutevi u trokutu 
        double alfa = acos ((b*b+c*c-a*a)/(2*b*c)) * 180.0/PI;
        double beta = acos ((a*a+c*c-b*b)/(2*a*c)) * 180.0/PI;
        double gamma = acos ((a*a+b*b-c*c)/(2*a*b)) * 180.0/PI;
   	
	if (alfa < amin ) amin = alfa;
	if (beta < amin ) amin = beta;
	if (gamma < amin ) amin = gamma; 
	
	if (alfa > amax ) amax = alfa;
	if (beta > amax ) amax = beta;
	if (gamma > amax ) amax = gamma;

    }

    std::cout << "Najmanji kut u mrezi je " << amin << std::endl;
    std::cout << "Najveci kut u mrezi je " << amax << std::endl;		

    std::cout << std::endl;

    Dune::VTKWriter<GridView> writer(pgrid->leafGridView());
    writer.write("alu_gmsh");
    return 0;
}
