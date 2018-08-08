#ifndef  __OPERATOR_HH_IS_INCLUDED__
#define  __OPERATOR_HH_IS_INCLUDED__

#include<dune/geometry/quadraturerules.hh>
#include<dune/geometry/referenceelements.hh>

#include<dune/pdelab/localoperator/defaultimp.hh>
#include<dune/pdelab/localoperator/flags.hh>
#include<dune/pdelab/localoperator/pattern.hh>
#include <dune/pdelab/finiteelement/localbasiscache.hh>
#include <dune/pdelab/localoperator/variablefactories.hh> // za makeJacobianContainer
#include "LE_exact.hh"
/** Lokalni operator za zadaću :
 *
 *
 * \tparam BCType tip rubnog uvjeta
 * \tparam FEM skalarni prostor konačnih elemenata
 */

template<typename BCType, typename FEM>
class ElasticityLocalOperator : // derivacijska lista -- jakobijan i pattern računa PDELab
  public Dune::PDELab::NumericalJacobianApplyVolume  <ElasticityLocalOperator<BCType,FEM>>,
  public Dune::PDELab::NumericalJacobianVolume       <ElasticityLocalOperator<BCType,FEM>>,
  public Dune::PDELab::NumericalJacobianApplyBoundary<ElasticityLocalOperator<BCType,FEM>>,
  public Dune::PDELab::NumericalJacobianBoundary     <ElasticityLocalOperator<BCType,FEM>>,
  public Dune::PDELab::FullVolumePattern,
  public Dune::PDELab::LocalOperatorDefaultFlags
{
public:
  // Zastavice koje signaliziraju da na svakom elementu treba zvati: 
  enum { doPatternVolume = true };  // metodu za računanje patterna (iz volumnih doprinosa)
  enum { doAlphaVolume = true };    // alpha_volume
  enum { doAlphaBoundary = true };  // alpha_boundary         

  ElasticityLocalOperator(const BCType& bctype_, // boundary cond.type
                          double mu_, double lambda_, double g_vert_, double rho_, bool exact_,
                          unsigned int intorder_=2) :
    bctype( bctype_ ), mu(mu_), lambda(lambda_), g_vert(g_vert_), rho(rho_), use_exact(exact_), intorder( intorder_ )
  {}

  // volume integral depending on test and ansatz functions
  // eg = element 
  // lfsu = lokalni prostor funkcija za rješenje
  // lfsv =  lokalni prostor funkcija za test funkciju
  // x    = vektor koeficijenata rješenja 
  // r    = lokalni rezidual
  template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
  void alpha_volume(const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
  {
    // dimensions
    const int dim = EG::Geometry::mydimension;
    const int dimw = EG::Geometry::coorddimension;

    // Koristimo činjenicu da je LFSU = LFSV
    // Tipovi skalarnih prostora 
//    using LFSU0 = typename LFSU::template Child<0>::Type;
//    using LFSU1 = typename LFSU::template Child<1>::Type;

     // uobičajene tipove uzimamo od skalarnog prostora
    //typedef typename LFSU0::Traits::FiniteElementType::Traits::LocalBasisType::Traits::DomainFieldType DF;
    //typedef typename LFSU0::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeFieldType  RF;
    //typedef typename LFSU0::Traits::FiniteElementType::Traits::LocalBasisType::Traits::JacobianType    Jacobian;
    //typedef typename LFSU0::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeType       Range;
    typedef Dune::FieldVector<double,dimw> Gradient;
    typedef typename LFSU::Traits::SizeType size_type;

    // select quadrature rule
    auto gt = eg.geometry().type();
    const auto& rule = Dune::QuadratureRules<double,dim>::rule(gt,intorder);

       // skalarni prostori
    auto const & lfsu0 = lfsu.template child<0>();
    auto const & lfsu1 = lfsu.template child<1>();
    auto const & lfsu2 = lfsu.template child<2>();

       // loop over quadrature points
       for (auto const & qp : rule)
         {
           // Imamo onoliko vektora baznih funkcija koliko ima komponenti,
           // ali svi prostori su isti pa su iste i bazne funkcije
//           std::vector<Range> phi0(lfsu0.size());
//           lfsu0.finiteElement().localBasis().evaluateFunction(qp.position(),phi0);
           auto& phi0 = cache.evaluateFunction(qp.position(),lfsu0.finiteElement().localBasis());

           // Izračunajmo sve komponenete u integracijskoj točki.
           double u_0=0.0, u_1=0.0,u_2 = 0.0;
           for (size_type i=0; i<lfsu0.size(); ++i) u_0 += x(lfsu0,i)*phi0[i];
           for (size_type i=0; i<lfsu1.size(); ++i) u_1 += x(lfsu1,i)*phi0[i];
          for (size_type i=0; i<lfsu2.size(); ++i) u_2 += x(lfsu2,i)*phi0[i];
           // Gradijenti baznih funkcija na ref elementu.
//           std::vector<Jacobian> js0(lfsu0.size());
//           lfsu0.finiteElement().localBasis().evaluateJacobian(qp.position(),js0);
           auto const & js0 = cache.evaluateJacobian(qp.position(),lfsu0.finiteElement().localBasis());
           // Gradijenti baznih funkcija na fizičkom elementu
           //const Dune::FieldMatrix<double,dimw,dim> & jac = ...
           const auto &jac = eg.geometry().jacobianInverseTransposed(qp.position());
           std::vector<Gradient> gradphi0(lfsu0.size());
           for (size_type i=0; i<lfsu0.size(); i++) jac.mv(js0[i][0],gradphi0[i]);

           // Gradijent komponenti rješenja (pomaka).
           Gradient gradu_0(0.0),  gradu_1(0.0),gradu_2(0.0);
           for (size_type i=0; i<lfsu0.size(); ++i) gradu_0.axpy(x(lfsu0,i), gradphi0[i]);
           for (size_type i=0; i<lfsu1.size(); ++i) gradu_1.axpy(x(lfsu1,i), gradphi0[i]);
           for (size_type i=0; i<lfsu2.size(); ++i) gradu_2.axpy(x(lfsu2,i), gradphi0[i]);

           // evaluate parameters;
           // Dune::FieldVector<RF,dim>
           auto x = eg.geometry().global(qp.position());
           // eg je ElementGeometry, zato moramo zvati entity metodu.
           Gradient f(0.0);
           f[0] = - 4 * lambda - 6 * mu;
           f[1] = 0;
           f[2] = 0;


           double divu = gradu_0[0]+gradu_1[1]+gradu_2[2];
           double D12u =  0.5*(gradu_1[0] + gradu_0[1] );
           double D13u =  0.5*(gradu_2[0] + gradu_0[2] );
           double D23u =  0.5*(gradu_2[1] + gradu_1[2] );

           // integrate grad u * grad phi_i + a*u*phi_i - f phi_i
           double factor = qp.weight()*eg.geometry().integrationElement(qp.position());

           for (size_type i=0; i<lfsu0.size(); ++i)
             r.accumulate(lfsu0, i, (  2*mu * ( gradu_0[0]*gradphi0[i][0] + D12u * gradphi0[i][1] +
                                        D13u * gradphi0[i][2] )+

                                       lambda *( divu * gradphi0[i][0] )
                                      - f[0]*phi0[i]
                                    ) * factor);
           for (size_type i=0; i<lfsu1.size(); ++i)
             r.accumulate(lfsu1, i, (  2*mu * ( gradu_1[1]*gradphi0[i][1] +  D12u * gradphi0[i][0] +
                                         D23u * gradphi0[i][2]) +

                                      lambda * ( divu * gradphi0[i][1] ) +
                                      - f[1]*phi0[i]
                                    ) * factor);
           for (size_type i=0; i<lfsu2.size(); ++i)
             r.accumulate(lfsu2, i, (  2*mu * ( gradu_2[2]*gradphi0[i][2] + D23u * gradphi0[i][1] +
                                     D13u * gradphi0[i][0]) +

                                       lambda * ( divu * gradphi0[i][2] )+
                                      - f[2]*phi0[i]
                                    ) * factor);
         }


  }

  // boundary integral
  // ig = intersection (= stranica elementa)
  // lfsu_s = lokalni prostor funkcija na stranici za rješenje
  // lfsu_v = lokalni prostor funkcija na stranici za test funkciju 
  // x_s    = vektor koeficijenata rješenja (na stranici)
  // r_s    = rezidual (na stranici)
  template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
  void alpha_boundary (const IG& ig, const LFSU& lfsu_s, const X& x_s,
                       const LFSV& lfsv_s, R& r_s) const
  {
    const int dim = IG::dimension;
    // Tipovi skalarnih prostora 
//    using LFSU0 = typename LFSU:: template Child<0>::Type;
//    using LFSU1 = typename LFSU:: template Child<dim-1>::Type;
    // skalarni prostori
    auto const & lfsu0 = lfsu_s.template child<0>();
    auto const & lfsu1 = lfsu_s.template child<1>();
    auto const & lfsu2 = lfsu_s.template child<2>();

    // some types
//    typedef typename LFSU0::Traits::FiniteElementType::Traits::LocalBasisType::Traits::DomainFieldType DF;
//    typedef typename LFSU0::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeFieldType RF;
//    typedef typename LFSU0::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeType Range;
    typedef typename LFSU::Traits::SizeType size_type;

    // dimensions

    // select quadrature rule for face
    auto gtface = ig.geometryInInside().type();
    auto const & rule = Dune::QuadratureRules<double,dim-1>::rule(gtface,intorder);

    // loop over quadrature points and integrate normal flux
    for (auto const & qp : rule)
      {
        // skip rest if we are on Dirichlet boundary
        if ( bctype.isDirichlet( ig, qp.position() ) )
          continue;
        // Global position
        auto globalpos = ig.geometry().global(qp.position());
        // Opterećenje imamo samo na gornjoj granici y=2.
       // if(globalpos[dim-1] > 2.0 - 1E-5)
         //     continue;
        // position of quadrature point in local coordinates of element
        auto local = ig.geometryInInside().global(qp.position());

        // Imamo onoliko vektora baznih funkcija koliko ima komponeanta, ali sve su "iste"
        auto& phi0 = cache.evaluateFunction(local,lfsu0.finiteElement().localBasis());
        // evaluate flux boundary condition

        // integrate j
        auto factor = qp.weight()*ig.geometry().integrationElement(qp.position());

        auto x = globalpos;
        Dune::FieldVector<double,dim> flux;
        double eps = 1E-6;
       if ( x[1] < eps) {          //y == 0
        flux[0] = - mu * x[1];
        flux[1] = (- 2 *mu - 4 * lambda) * x[0];
        flux[2] = 0.0;
       }
       else if ( std::abs(x[1]-2) < eps) {          //y == 1
        flux[0] =  mu * x[1];
        flux[1] = (  2 *mu +4 * lambda) * x[0];
        flux[2] = 0.0;
       }
       else if ( std::abs(x[0]-20) < eps) {          //x == 1
        flux[0] =  ( 4 *mu +4 * lambda) * x[0];
        flux[1] =  mu * x[1];
        flux[2] = mu * x[2];
       }
       else if ( std::abs(x[2]-2) < eps) {          //z == 1
        flux[0] =  mu * x[2];
        flux[1] =  0.0;
        flux[2] = (2 * mu + 4 * lambda) * x[0];
       }
       else if ( x[2] < eps) {          //z == 0
        flux[0] =  -mu * x[2];
        flux[1] =  0.0;
        flux[2] = -(2 * mu + 4 * lambda) * x[0];
       }



        for (size_type i=0; i<lfsu0.size(); ++i)
          r_s.accumulate(lfsu0,i, -flux[0] * phi0[i] * factor);

        for (size_type i=0; i<lfsu1.size(); ++i)
          r_s.accumulate(lfsu1,i, -flux[1] * phi0[i] * factor);

        for (size_type i=0; i<lfsu1.size(); ++i)
          r_s.accumulate(lfsu2,i, -flux[2] * phi0[i] * factor);

      }
  }

private:
  const BCType & bctype;
  double mu;
  double lambda;
  double g_vert;
  double rho;
  bool   use_exact;
  unsigned int intorder;
  typedef typename FEM::Traits::FiniteElementType::Traits::LocalBasisType LocalBasis;
  Dune::PDELab::LocalBasisCache<LocalBasis> cache;
};



#endif  
