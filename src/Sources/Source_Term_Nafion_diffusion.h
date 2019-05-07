/****************************************************************************
* Copyright (c) 2018, CEA
* All rights reserved.
*
* Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
* 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
* 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
* 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
* IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
* OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*
*****************************************************************************/
/////////////////////////////////////////////////////////////////////////////
//
// File      : Source_Term_Nafion_diffusion.h
// Directory : $PEMFC_ROOT/src/Sources
//
/////////////////////////////////////////////////////////////////////////////

#ifndef Source_Term_Nafion_diffusion_included
#define Source_Term_Nafion_diffusion_included

#include <Source_base.h>
#include <Probleme_base.h>
#include <Matrice_Morse.h>
#include <DoubleTab.h>
#include <Param.h>
#include <Ref_Domaine.h>
#include <Champ_Fonc.h>
#include <Discretisation_base.h>
#include <Ref_Zone_VF.h>
#include <Ref_Champ_Inc.h>
#include <Champ_Don.h>

/////////////////////////////////////////////////////////////////////////////
//
// .DESCRIPTION : class Source_Term_Nafion_diffusion
//  Modele de source de dissolution des especes (H2, N2, O2, H20 ou vap) de type:
//    Sa_i = D_i_naf * gamma_CL / e_naf * (Ceq_i - C_i) / ((1-por_naf)eps_naf)
//  avec
//    Ceq_i = c_i * R * T * H_i					avec i = H2 O2 N2
//    Ceq_i = C_SO3 * ld_eq(c_i*R*T/P_sat(T)) 	avec i = H2O
//    ld_eq (a) = 0.043 + 17.81*a -39.85a^2+36a^3
//    P_sat(T) = exp(23.1961-3816.44/(T-46.13))
// <Description of class Source_Term_Nafion_diffusion>
//
/////////////////////////////////////////////////////////////////////////////

class Source_Term_Nafion_diffusion : public Source_base
{

  Declare_instanciable( Source_Term_Nafion_diffusion ) ;

public :
  DoubleTab& ajouter(DoubleTab& ) const;
  DoubleTab& calculer(DoubleTab& ) const;
  void mettre_a_jour(double temps);
  void contribuer_a_avec(const DoubleTab&, Matrice_Morse&) const;
  void associer_zones(const Zone_dis& ,const Zone_Cl_dis& );
  void associer_pb(const Probleme_base& );
  void completer();
  void set_param(Param& param);
  void discretiser(const Discretisation_base& dis);
protected :
  Nom nom_espece_;						// nom d'espece dans la liste { H2 O2 N2 vap H2O}
  Nom nom_domaine_;
  Nom nom_ssz_CLa_;
  Nom nom_ssz_CLc_;
  Nom nom_pb_C_;
  Nom nom_champ_C_;
  Nom nom_champ_D_;
  Nom nom_pb_T_;
  Nom nom_champ_T_;
  Nom nom_pb_ci_cathode_;
  Nom nom_champ_ci_cathode_;
  Nom nom_pb_ci_anode_;
  Nom nom_champ_ci_anode_;

  REF(Domaine) dom_;						// domaine
  REF(Sous_Zone) CL_a_;						// sous_zone anode
  REF(Sous_Zone) CL_c_;						// sous_zone cathode
  REF(Champ_base)  ch_T_;   				// Champ Temperature de conduction de la chaleur   -> couple
  REF(Champ_base)  ch_D_i_naf_;				// champ conductivite Da -> get_champ
  REF(Champ_base)  ch_C_;					// Champ_Inc dissolved concentration -> inconnu()
  REF(Champ_base)  ch_ci_cathode_;			// Champ concentration de  diffusion des multi-especes -> couple
  REF(Champ_base)  ch_ci_anode_;			// Champ concentration de  diffusion des multi-especes -> couple

  Champ_Don eps_naf_;					// ionomer proportion
  Champ_Don por_naf_; 					// porosity
  Champ_Don gamma_CL_; 				    // specific surface m2/m3
  double T_0_;
  double C_SO3_;

  DoubleTab diffu_;
  DoubleTab C_;
  DoubleTab ci_;
  DoubleTab T_;

  REF(Zone_VF) la_zone_;

  Champ_Fonc ch_S_;			// champ de terme source P0

  double eval_f(double diffu, double Ci, double ci, double T, double por, double eps, double gamma) const;
  double eval_derivee_f(double diffu, double por, double eps, double gamma) const;

  inline double f_nd(double C) const;
  inline double f_Psat(double T) const;
  inline double f_lambda(double a) const;
  inline double f_Henry_H2(double T) const;
  inline double f_Henry_O2(double T) const;
  inline double f_Henry_N2(double T) const;
};
const double R = 8.314;
const double F = 96500;
const double a_lim = 1e-3;				// Activit minimale en dessous de laquelle son effet est linearise

inline double Source_Term_Nafion_diffusion::f_Psat(double T) const
{
  return exp(23.1961-3816.44/(T-46.13));
}

inline double Source_Term_Nafion_diffusion::f_lambda(double a) const
{
  return 0.043+17.81*a-39.85*a*a+36*a*a*a;
}

inline double Source_Term_Nafion_diffusion::f_Henry_H2(double T) const
{
  return (1./1.01325e5)*exp(9.05e3/(R*T));
}

inline double Source_Term_Nafion_diffusion::f_Henry_O2(double T)  const
{
  return (1/1.01325e5)*exp(5.88e3/(R*T));
}

inline double Source_Term_Nafion_diffusion::f_nd(double C) const
{
  double ld = C / C_SO3_;
  return 1. + 0.0028*ld + 0.0026*ld*ld;
}

inline double Source_Term_Nafion_diffusion::f_Henry_N2(double T) const
{
  return 6.4e-6*exp(1300*(1/T-1/298.15));
}
#endif /* Source_Term_Nafion_diffusion_included */
