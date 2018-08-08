#ifndef _EXACT_HH_
#define _EXACT_HH_

#include<dune/common/fvector.hh>

#include<dune/pdelab/common/function.hh>
#include<dune/pdelab/constraints/common/constraintsparameters.hh>

//Egzaktno rješenje
template <typename GV>
class Exact : public Dune::PDELab::AnalyticGridFunctionBase
          <
               Dune::PDELab::AnalyticGridFunctionTraits<GV,double,3>, // 1 = skalarna funkcija
           Exact<GV>
             >
{
public:
   typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,double,3> Traits;
   typedef Dune::PDELab::AnalyticGridFunctionBase<Traits, Exact<GV> > BaseT;

   Exact(GV const & gv) : BaseT(gv) {}
   void evaluateGlobal(typename Traits::DomainType const & x, typename Traits::RangeType & y) const
   {
       //y = DiffusionParameter<GV>::exact(x);
       y[0] =x[0] * x[0];
       y[1] =x[0] * x[1];
       y[2] = x[0] * x[2];
   }
};

#endif
