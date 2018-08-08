#ifndef DUNE_PDELAB_LOCALOPERATOR_DIFFUSIONMIXED_HH
#define DUNE_PDELAB_LOCALOPERATOR_DIFFUSIONMIXED_HH
 
#include <cstddef>
#include <vector>
 
#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
 
#include <dune/geometry/type.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/geometry/referenceelements.hh>

#include <dune/pdelab/common/quadraturerules.hh>
#include <dune/pdelab/common/referenceelements.hh>
 
#include <dune/pdelab/localoperator/defaultimp.hh>
#include <dune/pdelab/localoperator/pattern.hh>
#include <dune/pdelab/localoperator/flags.hh>

#include "bctype.hh"
 
# define M_PI           3.14159265358979323846


// Matrica A(x) kao funkcija globalne koordinate.

template <typename R, int dim>
Dune::FieldMatrix<R,dim,dim> A1 (Dune::FieldVector<R,dim> const & x_global)
{   double x = x_global[0];
    double y = x_global[1];
    double eps = 1e-6;
  Dune::FieldMatrix<R,dim,dim> mat;

      if ((std::abs(x - 0.5) < 0.2 + eps) && (std::abs(y - 0.5) < 0.2 + eps)) {
          for (std::size_t i=0; i<dim; i++)
            for (std::size_t j=0; j<dim; j++)
                mat[i][j] = (i==j) ? 100 : 0;
        }
      else{
          for (std::size_t i=0; i<dim; i++)
            for (std::size_t j=0; j<dim; j++)
                mat[i][j] = (i==j) ? 1 : 0;
    }
  return mat;
}


//za testiranje artificijelnog rjesenja
template <typename R, int dim>
Dune::FieldMatrix<R,dim,dim> A0 (Dune::FieldVector<R,dim> const & x_global)
{
  Dune::FieldMatrix<R,dim,dim> mat;
  for (std::size_t i=0; i<dim; i++)
    for (std::size_t j=0; j<dim; j++)
      mat[i][j] = (i==j) ? 1 : 0;

  return mat;
}


// slobodni koeficijent
template <typename R, int dim>
R a0(Dune::FieldVector<R,dim> const & x_global){
    return 0.0;
}

//desna strana
template <typename R, int dim>
R f1(Dune::FieldVector<R,dim> const & x_global){
    return 0;
}


// desna strana - za testiranje aritificijelnog rješenja
template <typename R, int dim>
R f0(Dune::FieldVector<R,dim> const & x_global){
    double x = x_global[0];
    double y = x_global[1];

    return (2 - x*(x-1)*(M_PI*M_PI)) * std::cos(M_PI * y);
}



     // Lokalni operator za skalarnu eliptičku rubnu zadaću:
     //     div p +a_0 u = f         in \Omega,
     //                p = -A grad u in \Omega,
     //                u = g         on \partial\Omega_D
     //              p.n = h         on \partial\Omega_N
     // sa H(div) elementima u mješovitoj formulaciji
     //
     // param.bctype : grid function type selecting boundary condition
     template<typename BCType>
     class DiffusionMixed : public Dune::PDELab::NumericalJacobianApplyVolume<DiffusionMixed<BCType> >,
                            public Dune::PDELab::NumericalJacobianVolume<DiffusionMixed<BCType> >,
                            public Dune::PDELab::FullVolumePattern,
                            public Dune::PDELab::LocalOperatorDefaultFlags
     {
 
//       using BCType = typename Dune::PDELab::ConvectionDiffusionBoundaryConditions::Type;
 
     public:
       // pattern assembly flags
       enum { doPatternVolume = true };
 
       // residual assembly flags
       enum { doAlphaVolume = true };
       enum { doLambdaVolume = true };
       enum { doLambdaBoundary = true };
 
       DiffusionMixed (const BCType& bctype_, int qorder_v_=2, int qorder_p_=1, int problem_ = 0)
         : bctype(bctype_), qorder_v(qorder_v_), qorder_p(qorder_p_), problem(problem_)
       {
       }
 
       template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
       void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
       {
         // Define types
         using VelocitySpace = typename LFSU::template Child<0>::Type;
         using PressureSpace = typename LFSU::template Child<1>::Type;

         using DF = typename VelocitySpace::Traits::FiniteElementType::Traits::LocalBasisType::Traits::DomainFieldType;
         using RF = typename VelocitySpace::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeFieldType;
         using VelocityJacobianType = typename VelocitySpace::Traits::FiniteElementType::Traits::LocalBasisType::Traits::JacobianType;
         using VelocityRangeType    = typename VelocitySpace::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeType;
         using PressureRangeType    = typename PressureSpace::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeType;
 
         // Prostori komponenti
         const auto& velocityspace = lfsu.template child<0>(); //child(lfsu,_0);
         const auto& pressurespace = lfsu.template child<1>();

         const int dim = EG::Geometry::mydimension;
//         const auto& cell = eg.entity();
         auto geo = eg.geometry();
 
         // Pretpostavljamo da je geometrisko preslikavanje afino pa B_K^{-tau} možemo izračunati bio gdje
         Dune::FieldVector<DF,dim> pos;
         pos=0.0;
         auto jac = geo.jacobianInverseTransposed(pos);
         jac.invert();   // B_K^\tau
         auto det = geo.integrationElement(pos);
 
         // Računamo matricu A
         auto tensor = A1(geo.center());
         if (problem == 0){
            tensor = A0(geo.center());
         }

            tensor.invert(); // B = A^{-1}


 
         // Vektori baze i transformirane baze (množene s B_K)
         std::vector<VelocityRangeType> v_basis(velocityspace.size());
         std::vector<VelocityRangeType> v_transformed_basis(velocityspace.size());

         // Vektorski dio rješenja
         VelocityRangeType p;  // rješenje p
         VelocityRangeType Bp; // Bp (B = A^{-1})
         std::vector<VelocityJacobianType> v_jac_basis(velocityspace.size());
         std::vector<PressureRangeType>    p_basis(pressurespace.size());
         std::vector<RF>                   divergence(velocityspace.size(),0.0);
 
 
         //  član Bp*q
         for (const auto& ip : quadratureRule(geo,qorder_v))
           {
             // Bazne funkcije na referentnom element: phi_i
             velocityspace.finiteElement().localBasis().evaluateFunction(ip.position(),v_basis);
 
             // Transformirane bazne funkcije B_K phi_i
             for (std::size_t i=0; i<velocityspace.size(); i++)
               {
                 v_transformed_basis[i] = 0.0;
                 jac.umtv(v_basis[i],v_transformed_basis[i]); // vtransformedbasis[i] += jac^t * v_basis[i] (= B_K * v_basis[i])
               }
 
             // vektorski dio rješenja, transformirani, odnosno množen s B_K
             p=0.0;
             for (std::size_t i=0; i<velocityspace.size(); i++)
               p.axpy(x(velocityspace,i), v_transformed_basis[i]);
 
             // B * p
             tensor.mv(p,Bp);
 
             // integrate  Bp * phi_i
             auto factor = ip.weight() / det;
             for (std::size_t i=0; i<velocityspace.size(); i++)
               r.accumulate(velocityspace,i,(Bp*v_transformed_basis[i])*factor);
           }
 
         //  Član u div q  i član (div p + a0*u)v
         for (const auto& ip : quadratureRule(geo, qorder_p))
           {
             // evaluate shape functions at ip (this is a Galerkin method)
             velocityspace.finiteElement().localBasis().evaluateJacobian(ip.position(), v_jac_basis);
             pressurespace.finiteElement().localBasis().evaluateFunction(ip.position(), p_basis);
 
             // Skalarni dio rješenja
             PressureRangeType u;
             u=0.0;
             for (std::size_t i=0; i<pressurespace.size(); i++)
               u.axpy(x(pressurespace,i), p_basis[i]);

             RF factor = ip.weight();
    /*
             // evaluate Helmholtz term (reaction term)
             auto a0value = a0(geo.global(ip.position()));
 
             // član a0 * u * v u drugoj jednadžbi

             for (std::size_t i=0; i<pressurespace.size(); i++)
               r.accumulate(pressurespace,i, -a0value * u * p_basis[i] * factor * det);
               //r.accumulate(pressurespace,i, -a0value * u * p_basis[i] * factor);
*/
             // divergencija baznih funkcija
             for (std::size_t i=0; i<velocityspace.size(); i++){
               divergence[i] = 0;
               for (int j=0; j<dim; j++)
                 divergence[i] += v_jac_basis[i][j][j];
             }
 
             // član u * div(q) u prvoj jednadžbi
             for (std::size_t i=0; i<velocityspace.size(); i++)
               r.accumulate(velocityspace,i, -u * divergence[i] * factor);
 
             // div(p)
             RF div_p = 0.0;
             for (std::size_t i=0; i<velocityspace.size(); i++)
               div_p += x(velocityspace,i)*divergence[i];
 
             // Član div p * v u drugoj jednadžbi
             for (std::size_t i=0; i<pressurespace.size(); i++)
               r.accumulate(pressurespace,i, - div_p * p_basis[i] * factor);
           }
       }
 
       // volumni integral desne strane u drugoj jednadžbi.
       template<typename EG, typename LFSV, typename R>
       void lambda_volume (const EG& eg, const LFSV& lfsv, R& r) const
       {
         // Define types
         using PressureSpace = typename LFSV::template Child<1>::Type;
         using PressureRangeType = typename PressureSpace::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeType;

         const auto& pressurespace = lfsv.template child<1>(); //child(lfsv,_1);
 //        const auto& cell = eg.entity();
         auto geo = eg.geometry();

         std::vector<PressureRangeType> p_basis(pressurespace.size());
 
         for (const auto& ip : quadratureRule(geo,qorder_p))
           {
             // bazne funkcije
             pressurespace.finiteElement().localBasis().evaluateFunction(ip.position(), p_basis);
 
             // desna strana
             auto ff = f1(geo.global(ip.position()));
             if (problem == 0)
                ff = f0(geo.global(ip.position()));


 
             // integrate f
             auto factor = ip.weight() * geo.integrationElement(ip.position());
             for (std::size_t i=0; i<pressurespace.size(); i++)
               r.accumulate(pressurespace,i, ff * p_basis[i] * factor);
           }
       }
 
       // integral po rubu koji dolazi od Dirichletovog uvjeta
       template<typename IG, typename LFSV, typename R>
       void lambda_boundary (const IG& ig, const LFSV& lfsv, R& r) const
       {
         using VelocitySpace = typename LFSV::template Child<0>::Type;
         using DF = typename VelocitySpace::Traits::FiniteElementType::Traits::LocalBasisType::Traits::DomainFieldType;
         using VelocityRangeType = typename VelocitySpace::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeType;
 
         const auto& velocityspace = lfsv.template child<0>(); //child(lfsv,_0);
         const int dim = IG::dimension;
         // Referenca na unutarnji element
         const auto& cell_inside = ig.inside();
 
         // geometrija stranice
         auto geo = ig.geometry();
         // geometrija elementa
         auto geo_inside = cell_inside.geometry();
 
         // Geometrija stranice u lokalnim koordinatama elementa
         auto geo_in_inside = ig.geometryInInside();
 
         // Pretpostavka: g_K je afino pa B_K izračunavamo bilo gdje
         Dune::FieldVector<DF,dim> pos;
         pos = 0.0;
         auto jac = geo_inside.jacobianInverseTransposed(pos);  // B_K^{-tau}
         jac.invert();                                          // B_K^tau
         auto det = geo_inside.integrationElement(pos);
 
         std::vector<VelocityRangeType> v_basis(velocityspace.size());
         std::vector<VelocityRangeType> v_transformed_basis(velocityspace.size());
 
         // loop over quadrature points and integrate normal flux
         for (const auto& ip : quadratureRule(geo,qorder_v))
           {
             // evaluate boundary condition type
             //auto bctype = param.bctype(ig.intersection(),ip.position());
 
             // Preskoči Neumannovu granicu
             if (bctype.isNeumann(ig.intersection(),ip.position()))
               continue;
 
             // pozicija kvadraturne točke u lokalnim koordinatama elementa
             auto local = geo_in_inside.global(ip.position());
 
             // Vektorske bazne funkcije
             velocityspace.finiteElement().localBasis().evaluateFunction(local,v_basis);
 
             // transformacija baznih funkcija
             for (std::size_t i=0; i<velocityspace.size(); i++)
               {
                 v_transformed_basis[i] = 0.0;
                 jac.umtv(v_basis[i],v_transformed_basis[i]);
               }
 
             // Vrijednost Dirichletovog rubnog uvjeta u integracijskoj točki stranice.
             auto gvalue = g(geo_inside.global(local));
 
             // integrate g v*normal
             auto factor = ip.weight()*geo.integrationElement(ip.position())/det;
             for (std::size_t i=0; i<velocityspace.size(); i++)
               r.accumulate(velocityspace,i,gvalue * (v_transformed_basis[i]*ig.unitOuterNormal(ip.position()))*factor);
           }
       }
 
     private:
       const BCType& bctype;
       int qorder_v;
       int qorder_p;
       int problem;
     };
 

 
 #endif // DUNE_PDELAB_LOCALOPERATOR_DIFFUSIONMIXED_HH
