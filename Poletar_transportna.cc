/*
* Rješavanje transportne jednadžbe u_t + q u_x = f na kružnoj domeni
Domena se učitava iz datoteke circle.msh
*/

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include <iostream>
#include <cstdlib>
#include <string>
#include <cmath>

#include <dune/common/parallel/mpihelper.hh>

#include <dune/alugrid/grid.hh>           // koristimo ALUGrid
#include <dune/grid/io/file/gmshreader.hh> // GmshReader klasa
#include <dune/grid/common/gridinfo.hh>    // Informacije o mreži.
#include <dune/grid/io/file/vtk.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/io/file/vtk/vtksequencewriter.hh>

#include <dune/istl/bvector.hh>
#include <dune/istl/bcrsmatrix.hh>


constexpr int dim = 2;  // dimenzija mreže (mora biti konstantna)
using Vector = Dune::FieldVector<double,1>;
using BlockVector = Dune::BlockVector<Vector>;

double omega = 0.2;

template <typename v, typename w>
double distance (v &x, w &y){
   double d = 0.0;
   for (int i = 0; i < dim; ++i)
        d += (x[i]-y[i])* (x[i]-y[i]);
   return std::sqrt(d);
}

template <typename GV>
void initialize (GV const &gv, BlockVector &c0){
    auto &set = gv.indexSet();
    std::vector <double> S(2);
    S[0]=0.5; S[1]=0.0;
    int broj_el = set.size(0);
    c0.resize(broj_el);
    int index_elementa;
    for(auto const & element : elements(gv)){
        index_elementa = set.index(element);
        //std::cout << "index el" << index_elementa << std::endl;
        auto centar = element.geometry().center();
        //std::cout << "udaljenost: " << distance (centar, S) << std::endl;

        if (distance (centar, S) < 0.25)
            c0[index_elementa]=1;
        else
            c0[index_elementa]=0;
    }
}

template <typename GV>
void nextStep (GV const &gv, double dt, double* dt_new, BlockVector & c_stari, BlockVector & c){
    Dune::FieldVector <double,2> q;
    auto &set = gv.indexSet();
    c.resize(set.size(0));
    double tmp;
    *dt_new = 1;

    for (auto const & element : elements(gv)){
       // petlja po svim stranicama elementa

        double vol_E = element.geometry().volume();
        int index_elementa = set.index(element);
        c[index_elementa] = c_stari[index_elementa];

        tmp = 0;
        for (auto const & side : intersections(gv, element))
        {
            if (side.boundary()== 0){ // Jesmo li na granici domene?

                int index_susjednog = set.index(side.outside());
                //vanjska normala na elementu
                auto N = side.centerUnitOuterNormal();
                N *= (side.geometry().volume()/vol_E);
                auto r = side.geometry().center();
                q[0] = -r[1] * omega;
                q[1] = r[0] * omega;
                double Nq = N*q;

                if (Nq >= 0){
                    tmp +=Nq;
                    c[index_elementa] -=dt*Nq*c_stari[index_elementa];
                }
                else
                    c[index_elementa] -=dt*Nq*c_stari[index_susjednog];
            }
        }
        if (*dt_new > 1/tmp) *dt_new = 1/tmp;
    }
    std::cout << "dt_new = " << *dt_new << std::endl;
}


int main(int argc, char** argv)
{
    Dune::MPIHelper::instance(argc, argv);

    if(argc != 3){
       std::cerr << "Usage: " << argv[0] << " ime_grid_datoteke.msh broj_koraka" << std::endl;
       std::exit(1);
    }

    int N = std::stoi(argv[2]); //broj koraka
    using GridType =  Dune::ALUGrid<2,2,Dune::simplex,Dune::conforming>;
    using GridView = GridType::LeafGridView;
    GridType* pgrid = Dune::GmshReader<GridType>::read(argv[1]);

    GridView gridView = pgrid->leafGridView();
    auto writer = std::make_shared<Dune::VTKWriter<GridView> > (gridView);
    Dune::VTKSequenceWriter<GridView> vtkSequnceWriter(writer,"transportna");

    BlockVector c0, c1;
    double time = 0.0;
    initialize (gridView, c0);
    writer->addCellData(c0,"koncentracija");
    vtkSequnceWriter.write(time); //pocetni trenutak
    double delta_t = 0.1;
    double delta_t_new;

    for (int i = 0 ; i < N ; i++)
    {
        nextStep (gridView, delta_t, &delta_t_new, c0, c1);
        c0 = c1;
        time += delta_t;
        vtkSequnceWriter.write(time);
        delta_t = delta_t_new;
    }
    return 0;
}
