#ifndef _EXACT_HH_
#define _EXACT_HH_

#include<dune/common/fvector.hh>

#include<dune/pdelab/common/function.hh>
#include<dune/pdelab/constraints/common/constraintsparameters.hh>
# define M_PI           3.14159265358979323846

//Klasa koja implementira artificijelno rješenje
// u (x,y) = x + x * (1-x) * cos (pi * y)

template <typename GV>
class ExactPressure : public Dune::PDELab::AnalyticGridFunctionBase
          <
               Dune::PDELab::AnalyticGridFunctionTraits<GV,double,1>, // 1 = skalarna funkcija
           ExactPressure<GV>
             >
{
public:
   typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,double,1> Traits;
   typedef Dune::PDELab::AnalyticGridFunctionBase<Traits, ExactPressure<GV> > BaseT;

   ExactPressure(GV const & gv) : BaseT(gv) {}
   void evaluateGlobal(typename Traits::DomainType const & x, typename Traits::RangeType & y) const
   {
    y = x[0] + x[0] * (1-x[0]) * std::cos(M_PI*x[1]);
   }
};

template <typename GV>
class ExactVelocity : public Dune::PDELab::AnalyticGridFunctionBase
          <
               Dune::PDELab::AnalyticGridFunctionTraits<GV,double,2>, // 1 = skalarna funkcija
           ExactVelocity<GV>
             >
{
public:
   typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,double,2> Traits;
   typedef Dune::PDELab::AnalyticGridFunctionBase<Traits, ExactVelocity<GV> > BaseT;

   ExactVelocity(GV const & gv) : BaseT(gv) {}
   void evaluateGlobal(typename Traits::DomainType const & x, typename Traits::RangeType & y) const
   {
       y[0] = (2 * x[0] -1) * std::cos(M_PI * x[1]) - 1;
       y[1] = M_PI * std::sin(M_PI * x[1]) * x[0] * (1-x[0]);

   }
};

#endif
