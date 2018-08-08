#ifndef _BCTYPE_HH_
#define _BCTYPE_HH_

#include<dune/common/fvector.hh>

#include<dune/pdelab/common/function.hh>
#include<dune/pdelab/constraints/common/constraintsparameters.hh>

#define PI  3.14159265359

/*  Klasa koja određuje vrijednost Dirichletovog rubnog uvjeta i 
    njegovo proširenje na čitavu domenu. 
    Template parametri:
       GV = GridView
       RF = Range Field Type (tip kojim su predstavljeni elementi slike funkcije)
       
       Ovo je rubni uvjet ua skalarnu komponentu.
    */
template<typename GV, typename RF>
class BCExtension
  : public Dune::PDELab::GridFunctionBase<Dune::PDELab::
           GridFunctionTraits<GV,RF,1,Dune::FieldVector<RF,1> >, BCExtension<GV,RF> > {
  const GV& gv;
public:
  typedef Dune::PDELab::GridFunctionTraits<GV,RF,1,Dune::FieldVector<RF,1> > Traits;

  // Konstruktor samo uzima referencu na  GridView objekt. 
  BCExtension (const GV& gv_) : gv(gv_) {}

  // Izračunaj Dirichletovu vrijednost na elementu. Ako točka nije na 
  // Dirichletovoj granici, onda funkcija daje proširenje Dirichletovog rubnog
  // uvjeta na čitavu domenu. To je proširenje u osnovi proizvoljno. 
  // e     = element 
  // xlocal = lokalne koordinate točke u kojoj se računa Dirichletova vrijednost
  // y      = izračunata Dirichletova vrijednost
  inline void evaluate (const typename Traits::ElementType& e,
                        const typename Traits::DomainType& xlocal,
                        typename Traits::RangeType& y) const
  {
//    const int dim = Traits::GridViewType::Grid::dimension;
//    typedef typename Traits::GridViewType::Grid::ctype ctype;

    // Pretvori lokalne koordinate u globalne
    //Dune::FieldVector<ctype,dim> x = e.geometry().global(xlocal);
    y = 0.0;
    return;
  }

  // Vrati referencu na GridView
  inline const GV& getGridView () {return gv;}
};


// Ovo je klasa za skalarnu komponentu
// što znači da pretpostavljamo isti tip uvjeta za sve komponente.
template <typename GV>
class BCTypeParam : public Dune::PDELab::DirichletConstraintsParameters
{
  const GV& gv;
public:
//  typedef GV GridViewType;

  BCTypeParam( const GV& gv_ ) : gv(gv_) { }

public:
  //  intersection = stranica elementa (u 3D) ili brid elementa (u 2D)
  //  coord        = lokalne koordinate točke na "intersectionu" koja se ispituje
  //  povratna vrijednost: true ako je točka na Dirichletovoj granici
  //                       false ako nije. 
  template<typename I>
  bool isDirichlet(const I & intersection,
                   const Dune::FieldVector<typename I::ctype, I::dimension-1> & coord
                   ) const
  {
//  Globalne koordinate točke (uočite da su dimenzije lokalne i globalne točke različite )
    Dune::FieldVector<typename I::ctype, I::dimension> xg = intersection.geometry().global( coord );
//
    if( xg[0] < 1E-6) return true;

    return false;  
  }

  //! get a reference to the grid view
  inline const GV& getGridView () 
  {
    return gv;
  }

};

#endif
