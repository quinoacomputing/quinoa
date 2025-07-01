#include <iostream>
#include <cmath>

//! Create struct of original EOS class that only contains
//! member variable and a method to compute min_eff_pressure
//! which is the only computation needed for Volume integral
//! Varaint style polymorphism is not friendly, and recreating
//! union-style polymorphism that contains all the variables, and
//! member function is not required for Volume integral

namespace EOS_Kokkos {

    struct GodunovRomenski {
        tk::real m_gamma, m_mu, m_rho0, m_alpha, m_K0;

        //! need default constructor for View initilization
        KOKKOS_INLINE_FUNCTION
        GodunovRomenski() = default;

        //! basic constructor
        KOKKOS_INLINE_FUNCTION
        GodunovRomenski(tk::real gamma, tk::real mu, tk::real rho0,
            tk::real alpha, tk::real K0) : m_gamma(gamma), m_mu(mu), m_rh0(rh0), 
                m_alpha(alpha), m_K0(K0);
        
        //! Directly copied frm GodunovRomenski.cpp, with some subsitution to
        //! avoid defining more functions in the struct
        //! Conglomerate pressure_cold_compr and DpccDrho
        KOKKOS_INLINE_FUNCTION
        tk::real min_eff_pressure(tk::real min, tk::real arho, tk::real alpha ) const {
        // minimum pressure is constrained by zero soundspeed.
        auto rho = arho/alpha;

        auto rrho0a = std::pow(rho/m_rho0, m_alpha);

        return min - rho/(m_gamma+1.0) * (m_K0/(m_rho0*m_alpha) *
        ((2.0*m_alpha+1.0)*(rrho0a*rrho0a) - (m_alpha+1.0)*rrho0a))
        + (m_K0/m_alpha * (rrho0a*rho/m_rho0) * (rrho0a-1.0));
        }
    }

    struct WL {
        
        tk::real m_w, m_cv, m_rho0, m_de, m_rhor, m_tr, m_pr, m_a, m_b, m_r1, m_r2;
        //! need default constructor for View initilization
        KOKKOS_INLINE_FUNCTION
        JWL() = default;

        //! basic constructor
        KOKKOS_INLINE_FUNCTION
        JWL( tk::real w, tk::real cv, tk::real rho0, tk::real de, tk::real rhor,
        tk::real tr, tk::real pr, tk::real A, tk::real B, tk::real R1,
        tk::real R2 ) : m_w(w), m_cv(cv), m_rho0(rho0), m_de(de), m_rhor(rhor),
            m_tr(tr), m_pr(pr), m_a(A), m_b(B), m_r1(R1), m_r2(R2) {};

        KOKKOS_INLINE_FUNCTION
        tk::real min_eff_pressure(tk::real min, tk::real arho, tk::real alpha ) const {

            alpha = std::max(1e-14,alpha);
            auto co1 = m_rho0*alpha*alpha/(arho*arho);
            auto co2 = alpha*(1.0+m_w)/arho;

            // minimum pressure is constrained by zero soundspeed.
            auto apr = -(m_a*(m_r1*co1 - co2) * exp(-m_r1*alpha*m_rho0/arho)
                        + m_b*(m_r2*co1 - co2) * exp(-m_r2*alpha*m_rho0/arho))
                * arho/(1.0+m_w);
            return apr/alpha + min;
        }
    }

    struct SmallShearSolid {
        tk::real m_gamma, m_pstiff, m_cv, m_mu;

        //! need default constructor for View initilization
        KOKKOS_INLINE_FUNCTION
        SmallShearSolid() = default;

        KOKKOS_INLINE_FUNCTION
        SmallShearSolid(tk::real gamma, tk::real pstiff, tk::real cv,tk::real mu ) :
            m_gamma(gamma), m_pstiff(pstiff), m_cv(cv), m_mu(mu) {};
        
        KOKKOS_INLINE_FUNCTION
        tk::real min_eff_pressure(tk::real min, tk::real, tk::real ) const {
        // minimum pressure is constrained by zero soundspeed.
            return (min - m_pstiff);
        }
    }

    struct StiffenedGas {
    
        tk::real m_gamma, m_pstiff, m_cv;

        //! need default constructor for View initilization
        KOKKOS_INLINE_FUNCTION
        StiffenedGas() = default;

        //! basic constructor
        KOKKOS_INLINE_FUNCTION
        StiffenedGas(tk::real gamma, tk::real pstiff, tk::real cv ) :
            m_gamma(gamma), m_pstiff(pstiff), m_cv(cv) {};
        
        KOKKOS_INLINE_FUNCTION
        tk::real min_eff_pressure(tk::real min, tk::real, tk::real ) const {
            // minimum pressure is constrained by zero soundspeed.
            return (min - m_pstiff);
        }
    }
    struct ThermallyPerfectGas {

        tk::real m_R;
        Kokkos::Array<Kokkos::Array<tk::real, 8>, 3> m_cp_coeff = {};
        Kokkos::Array<tk::real, 4> m_t_range = {};
        tk::real m_dH_ref;

        KOKKOS_INLINE_FUNCTION
        ThermallyPerfectGas() = default();

        KOKKOS_INLINE_FUNCTION
        ThermallyPerfectGas(tk::real R, std::vector< std::vector< tk::real >> cp_coeff,
            std::vector< tk::real > t_range, tk::real dH_ref)  {
            m_R = R;
            m_dH_ref = dH_ref;
            for (int i = 0; i < 3;i++) {
                for(int j = 0; j < 8;j++) {
                    m_cp_coeff[i][j] = cp_coeff[i][j]
                }
                m_t_range[i] = t_range[i];
            }
            m_t_range[3] = t_range[3];
            
        }
        KOKKOS_INLINE_FUNCTION
        tk::real min_eff_pressure(tk::real min, tk::real, tk::real ) const
            { return min; }
    }

    struct WilkinsAluminum {
        tk::real m_gamma, m_cv, m_mu, m_rho0;
        
        KOKKOS_INLINE_FUNCTION WilkinsAluminum(tk::real gamma, tk::real cv, tk::real mu ) :
        m_gamma(gamma), m_cv(cv), m_mu(mu) {};

        KOKKOS_INLINE_FUNCTION
        tk::real min_eff_pressure(tk::real min, tk::real, tk::real ) const {
        // minimum pressure is constrained by zero soundspeed.
            return min;
        }
    }
    class EOS_Kokkos {
        private:
            //! Using union style instead of variant
            enum class EOSType = {StiffenedGas
                        , JWL
                        , SmallShearSolid
                        , WilkinsAluminum
                        , GodunovRomenski
                        , ThermallyPerfectGas};

            EOSType type;
            union EOSUnion {
                StiffiedGas stiffenedGas;
                JWl jwl;
                SmallShearSolid smallShearSolid;
                WilkinsAluminum wilkinsAluminum;
                GodunovRomenski godunovRomenski;
                ThermallyPerfectGas thermallyPerfectGas;

                KOKKOS_INLINE_FUNCTION
                EOSUnion() {};

                KOKKOS_INLINE_FUNCTION
                ~EOSUnion() {};

                KOKKOS_INLINE_FUNCTION
                EOSUnion(const StiffenedGas& s) : stiffenedGas(s) {};

                KOKKOS_INLINE_FUNCTION
                EOSUnion(const JWL& s) : jwl(s) {};

                KOKKOS_INLINE_FUNCTION
                EOSUnion(const SmallShearSolid& s) : smallShearSolid(s) {};

                KOKKOS_INLINE_FUNCTION
                EOSUnion(const WilkinsAluminum& s) : wilkinsAluminum(s) {};

                KOKKOS_INLINE_FUNCTION
                EOSUnion(const GodunovRomenskis& s) : godunovRomenski(s) {};

                KOKKOS_INLINE_FUNCTION
                EOSUnion(const ThermallyPerfectGas& s) : thermallyPerfectGas(s) {};

            } m_material;

        public:

            KOKKOS_INLINE_FUNCTION
            EOS_Kokkos(): type(EOSType::StiffenedGas), m_material(StiffenedGas()) {};

            KOKKOS_INLINE_FUNCTION
            EOS_Kokkos(const StiffenedGas& s) : type(EOSType::StiffenedGas), m_material(s) {};

            KOKKOS_INLINE_FUNCTION
            EOS_Kokkos(const JWL& s) : type(EOSType::JWL), m_material(s) {};
            
            KOKKOS_INLINE_FUNCTION
            EOS_Kokkos(const SmallShearSolid& s) : type(EOSType::SmallShearSolid), m_material(s) {};
            
            KOKKOS_INLINE_FUNCTION
            EOS_Kokkos(const WilkinsAluminum& s) : type(EOSType::WilkinsAluminum), m_material(s) {};
            
            KOKKOS_INLINE_FUNCTION
            EOS_Kokkos(const GodunovRomenski& s) : type(EOSType::GodunovRomenski), m_material(s) {};
            
            KOKKOS_INLINE_FUNCTION
            EOS_Kokkos(const ThermallyPerfectGas& s) : type(EOSType::ThermallyPerfectGas), m_material(s) {};
    }
}