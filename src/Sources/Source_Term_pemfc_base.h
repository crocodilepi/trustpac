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
// File      : Source_Term_pemfc_base.h
// Directory : $PEMFC_ROOT/src/Sources
//
/////////////////////////////////////////////////////////////////////////////

#ifndef Source_Term_pemfc_base_included
#define Source_Term_pemfc_base_included

#include <Source_base.h>
#include <Ref_Champ_base.h>
#include <Probleme_base.h>
#include <Matrice_Morse.h>
#include <DoubleTab.h>
#include <Zone_dis.h>
#include <Zone_Cl_dis.h>
#include <Param.h>
#include <Ref_Domaine.h>
#include <Ref_Operateur_base.h>

/////////////////////////////////////////////////////////////////////////////
//
// .DESCRIPTION : class Source_Term_pemfc_base
//  Modele de source de dissolution des especes (H2, N2, O2, H20 ou vap) de type:
//  Sa_i = D_i_naf * gamma_CL / e_naf * (Ceq_i - C_i) + S_i
//  avec
//  Ceq_i = c_i * R * T * H_i		avec i = H3 O2 N2
//  Ceq_i = C_SO3 * ld_eq(c_i*R*T/P_sat(T)) 	avec i = H2O
//  ld_eq (a) = 0.043 + 17.81*a -39.85a^2+36a^3
//  P_sat(T) = exp(23.1961-3816.44/(T-46.13))
//  S_H2 = -ir / (2F)
//  S_O2 = (ir+ip) / (4F)
//  S_N2 = 0
//  S_H20 = -(ir+ip) / (2F) dans CL_c et null dans les zones restants
//  En particulier avec H20, un terme source supplementaire :
//  S_H20_supp = -nd/F*div(kappa.grad(phi)) = -nd/F*Operateur_diff(0)
// <Description of class Source_Term_pemfc_base>
//
/////////////////////////////////////////////////////////////////////////////

class Source_Term_pemfc_base : public Source_base
{

  Declare_base( Source_Term_pemfc_base ) ;

public :
  virtual DoubleTab& ajouter(DoubleTab& ) const =0;
  virtual DoubleTab& calculer(DoubleTab& ) const;
  virtual void mettre_a_jour(double temps);
  virtual void contribuer_a_avec(const DoubleTab&, Matrice_Morse&) const =0;
  //virtual void associer_zones(const Zone_dis& ,const Zone_Cl_dis& )=0;
  //virtual void associer_pb(const Probleme_base& )=0;
  virtual void completer();
  //void set_param(Param& param);
protected :

  Nom nom_espece_;						// nom d'espece dans la liste { H2 O2 N2 vap H2O}
  Nom nom_domaine_;
  Nom nom_ssz_CLa_;
  Nom nom_ssz_CLc_;
  Nom nom_champ_D_;
  Nom nom_pb_T_;
  Nom nom_champ_T_;
  Nom nom_pb_ci_cathode_;
  Nom nom_champ_ci_cathode_;
  Nom nom_pb_ci_anode_;
  Nom nom_champ_ci_anode_;
  Nom nom_pb_phi_;
  Nom nom_champ_ir_;
  Nom nom_champ_ip_;
  Nom nom_op_diff_;

  REF(Domaine) dom_;					// domaine
  REF(Sous_Zone) CL_a_;					// sous_zone anode
  REF(Sous_Zone) CL_c_;					// sous_zone cathode
  REF(Champ_base)  ch_T_;   				// Champ Temperature de conduction de la chaleur   -> couple
  REF(Champ_base)  ch_D_i_naf_;				// champ conductivite Da -> get_champ
  REF(Champ_Inc)   ch_C_;					// Champ_Inc dissolved concentration -> inconnu()
  REF(Champ_base)  ch_ci_cathode_;			// Champ concentration de  diffusion des multi-especes -> couple
  REF(Champ_base)  ch_ci_anode_;			// Champ concentration de  diffusion des multi-especes -> couple
  REF(Champ_base)  ch_ir_;	// couple
  REF(Champ_base)  ch_ip_;	// couple
  REF(Champ_base)  ch_op_;	// couple

  double eps_naf_;					// ionomer proportion
  double por_naf_; 					// porosity
  double gamma_CL_; 				// specific surface m2/m3
  double T_0_;
  double C_SO3_;

  DoubleTab diffu_;
  DoubleTab C_;
  DoubleTab ci_;
  DoubleTab T_;
  DoubleTab ir_;
  DoubleTab ip_;
  DoubleTab op_;

  DoubleVect volumes_;

  virtual void remplir_volumes()=0;
  double eval_f(double diffu, double Ci, double ci, double T) const;
  double eval_derivee_f(double diffu) const;

  double f_nd(double C) const;
  double f_Psat(double T) const;
  double f_lambda(double a) const;
  double f_Henry_H2(double T) const;
  double f_Henry_O2(double T) const;
  double f_Henry_N2(double T) const;
};

const double R = 8.314;
const double F = 96500;

inline double Source_Term_pemfc_base::f_Psat(double T) const
{
  return exp(23.1961-3816.44/(T-46.13));
}

inline double Source_Term_pemfc_base::f_lambda(double a) const
{
  return 0.043+17.81*a-39.85*a*a+36*a*a*a;
}

inline double Source_Term_pemfc_base::f_Henry_H2(double T) const
{
  return (1./1.01325e5)*exp(9.05e3/(R*T));
}

inline double Source_Term_pemfc_base::f_Henry_O2(double T)  const
{
  return (1/1.01325e5)*exp(5.88e3/(R*T));
}

inline double Source_Term_pemfc_base::f_nd(double C) const
{
  double ld = C / C_SO3_;
  return 1. + 0.0028*ld + 0.0026*ld*ld;
}

inline double Source_Term_pemfc_base::f_Henry_N2(double T) const
{
  return 6.4e-6*exp(1300*(1/T-1/298.15));
}

#endif /* Source_Term_pemfc_base_included */
