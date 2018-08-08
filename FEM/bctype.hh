#ifndef BCTYPE_HH
#define BCTYPE_HH

#include <cmath>
#include <dune/pdelab/common/function.hh>
#include <dune/pdelab/constraints/common/constraints.hh>
#include <dune/pdelab/constraints/common/constraintsparameters.hh>

// Na tlak nema nikakvog rubnog uvjeta ugrađenog u prostor.
// Dirichletov rubni uvjet na tlak se zadovoljava varijacijski
// kroz funkciju g (vidi u operator.hh)
class BCTypePressure
  : public Dune::PDELab::DirichletConstraintsParameters
{
public:

  template<typename I>
  bool isDirichlet(
                   const I & intersection,
                   const Dune::FieldVector<typename I::ctype, I::dimension-1> & coord
                   ) const
  {

        return false;
  }
};


class BCTypeVelocity : public Dune::PDELab::DirichletConstraintsParameters
{
public:

  template<typename I>
  bool isDirichlet(const I & ig
                   , const Dune::FieldVector<typename I::ctype, I::dimension-1> & coord
                   ) const
  {
      auto xg = ig.geometry().global( coord );

      if (xg[0] == 0 || xg[0] == 1)
        return true;
      else
        return false;

  }

  template<typename I>
  bool isNeumann(const I & ig,
                 const Dune::FieldVector<typename I::ctype, I::dimension-1> & coord
                 ) const
  {
    return !isDirichlet( ig, coord );
  }

};

// Dirichlet rubni uvjet za tlak.
template <typename R, int dim>
R g (Dune::FieldVector<R,dim> const & x_global)
{

    if (x_global[0] == 1)
        return 1;
    else
        return 0;

}

template <typename GV>
class VelocityExtension
  : public Dune::PDELab::AnalyticGridFunctionBase<Dune::PDELab::AnalyticGridFunctionTraits<GV,double,GV::dimension>,
                                                  VelocityExtension<GV> >
{
public:
  typedef Dune::PDELab::AnalyticGridFunctionTraits<GV, double, GV::dimension> Traits;
  typedef Dune::PDELab::AnalyticGridFunctionBase<Traits,VelocityExtension<GV> > BaseT;


  //! constructor
  VelocityExtension(const typename Traits::GridViewType& gv_)
    : BaseT(gv_), gv(gv_)
  {}

  inline void evaluateGlobal (const typename Traits::DomainType& x,
                              typename Traits::RangeType& y) const
  {

        y = g(x);
  }

  inline void evaluate (const typename Traits::ElementType& e,
                        const typename Traits::DomainType& x,
                        typename Traits::RangeType& y) const
  {

      auto xglobal = e.geometry().global(x);
      y = g(xglobal);
  }

  inline const GV& getGridView () const { return gv; }

private:
  const GV & gv;
};





template<typename GV>
class PressureExtension
: public Dune::PDELab::GridFunctionBase<Dune::PDELab::GridFunctionTraits<GV, double, 1, Dune::FieldVector<double,1> >,
                                        PressureExtension<GV> >
{
public:
typedef Dune::PDELab::GridFunctionTraits<GV,double,1,Dune::FieldVector<double,1> > Traits;

PressureExtension(const GV& gv_) : gv(gv_) {}

inline void evaluateGlobal (const typename Traits::DomainType& x,
                            typename Traits::RangeType& y) const
{
  y = g(x);
}

inline void evaluate (const typename Traits::ElementType& e,
                      const typename Traits::DomainType& x,
                      typename Traits::RangeType& y) const
{
    auto xglobal = e.geometry().global(x);
    y = g(xglobal);
}

inline const GV & getGridView () const { return gv; }

private:
const GV & gv;
};


#endif // BCTYPE_HH
