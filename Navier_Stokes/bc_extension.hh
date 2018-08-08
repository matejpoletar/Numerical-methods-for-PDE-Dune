#ifndef CG_STOKES_INITIAL_HH
#define CG_STOKES_INITIAL_HH

#include <dune/common/fvector.hh>

#include <dune/pdelab/localoperator/stokesparameter.hh>
#include <dune/pdelab/common/function.hh>

//===============================================================
// Define parameter functions f,g,j and \partial\Omega_D/N
//===============================================================


// Klasa koja određuje tip granice
class BCTypeParam
{
private:
    double time;
public:
    // Ova klasa daje indekse: DoNothing i VelocityDirichlet i StressNeumann 
  typedef Dune::PDELab::StokesBoundaryCondition BC;

  struct Traits
  {
    typedef BC::Type RangeType;
  };

  // Domena se po osi x prostire od -5 do 15. Na granici x=15 ne postavljamo nikakve uvjete
  // (= homogeni Neumannovi uvjeti). Na ostatku granice imamo Dirichletove uvjete na brzinu. 
  template<typename I>
  inline void evaluate (const I & intersection,   
                        const Dune::FieldVector<typename I::ctype, I::dimension-1> & coord,
                        BC::Type& y) const
  {
    Dune::FieldVector<typename I::ctype, I::dimension>
        xg = intersection.geometry().global( coord );
    if(xg[0] > 15 - 1e-6)
      y = BC::DoNothing;
    else
      y = BC::VelocityDirichlet;
  }
  template <typename T>
  void setTime(T t){
    time = t;
  }
};


// Inicijalni i Dirichletov  rubni uvjet za brzinu 
template<typename GV, typename RF, int dim>
class Velocity :
  public Dune::PDELab::AnalyticGridFunctionBase<
                             Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,dim>,
                             Velocity<GV,RF,dim> 
                                               >
{
private:
  RF time;

public:
  typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,dim> Traits;
  typedef Dune::PDELab::AnalyticGridFunctionBase<Traits,Velocity<GV,RF,dim> > BaseT;

  typedef typename Traits::DomainType DomainType;
  typedef typename Traits::RangeType RangeType;

  Velocity(const GV & gv) : BaseT(gv) {
    time = 0.0;
  }

  inline void evaluateGlobal(const DomainType & x, RangeType & y) const
  {
      y[1] = 0;  // vertikalna brzina je svugdje nula

      if(x[0] < -1)
        y[0] = 1;
      else
        y[0]=0;

  }

  template <typename T>
  void setTime(T t){
    time = t;
  }

};

// Vektorska funkcija jednaka nuli. Ovdje je koristimo  za tlak (dim_range=1)
// i za Neumannov rubni uvjet (dim_range=dim) koji nije prisutan pa ga stavljamo na nulu.
// Rubni uvjet za tlak jednako tako ne postoji i stoga vrijednost tlaka stavljamo na nulu.
template<typename GV, typename RF, std::size_t dim_range>
class ZeroFunction :
  public Dune::PDELab::AnalyticGridFunctionBase<
                            Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,dim_range>,
                            ZeroFunction<GV,RF,dim_range> 
                                               >,
  public Dune::PDELab::InstationaryFunctionDefaults
{
public:
  typedef Dune::PDELab::AnalyticGridFunctionTraits<GV,RF,dim_range> Traits;
  typedef Dune::PDELab::AnalyticGridFunctionBase<Traits, ZeroFunction> BaseT;

  typedef typename Traits::DomainType DomainType;
  typedef typename Traits::RangeType RangeType;

  ZeroFunction(const GV & gv) : BaseT(gv) {}

  inline void evaluateGlobal(const DomainType & x, RangeType & y) const
  {
    y=0;
  }
};



#endif
