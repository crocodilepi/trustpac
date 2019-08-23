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
// File      : Young_Todd_Source_multi.h
// Directory : $PEMFC_ROOT/src/Sources
//
/////////////////////////////////////////////////////////////////////////////

#ifndef Young_Todd_Source_multi_included
#define Young_Todd_Source_multi_included

#include <Source_base.h>
#include <Probleme_base.h>
#include <Matrice_Morse.h>
#include <DoubleTab.h>
#include <Param.h>
#include <Ref_Domaine.h>
#include <Champ_Fonc.h>
#include <Discretisation_base.h>
#include <Ref_Zone_VF.h>
#include <Champ_Don.h>

/////////////////////////////////////////////////////////////////////////////
//
// .DESCRIPTION : class Young_Todd_Source_multi
// Impliciter le terme source S = { Sa_X2, Sa_vap, Sa_N2 } dans la sous zone (CL)
// avec
// Sa_X2 = D_X2.gamma_CL/(ep_naf) (C_X2 - c_X2.R.T.H_X2)
// Sa_vap = D_vap.gamma_CL/(ep_naf) (C_vap - C_SO3.f_ld_eq(c_vap.R.T/f_Psat(T))
// Sa_N2 = D_N2.gamma_CL/(ep_naf) (C_N2 - c_N2.R.T.H_N2)
// ep_naf = (1-por_naf)eps_naf/gamma_CL
// La formulation pour le terme source implicite est la suivante:
// S_(t+1) = S_t + dS/dt(c_(t+1) - c_t) donc
// terme1 = -dS/dt.c_(t+1) dans contribuer_a_avec()
// terme2 = S_t - dS/dt.c_t dans ajouter()
// <Description of class Young_Todd_Source_multi>
//
/////////////////////////////////////////////////////////////////////////////

class Young_Todd_Source_multi : public Source_base
{

  Declare_instanciable( Young_Todd_Source_multi ) ;

public :
  DoubleTab& ajouter(DoubleTab& ) const;
  DoubleTab& calculer(DoubleTab& ) const;
  void mettre_a_jour(double temps);
  void contribuer_a_avec(const DoubleTab&, Matrice_Morse&) const;
  void associer_zones(const Zone_dis& ,const Zone_Cl_dis& );
  void associer_pb(const Probleme_base& );
  void completer();
  void set_param(Param& param);
protected :
  REF(Zone_VF) la_zone_;
  Nom nom_pb_ci_;
  Nom nom_pb_C_X2_;			// pb de dissous X2 = O2/H2 dans Nafion
  Nom nom_pb_C_vap_;			// pb de dissous vap dans Nafion
  Nom nom_pb_C_N2_;			    // pb de dissous N2 dans Nafion
  Nom nom_pb_T_;
  Nom nom_champ_ci_;		    //
  Nom nom_champ_C_X2_;		    // champ source couple
  Nom nom_champ_D_X2_;			// diffusivity X2 dans Nafion
  Nom nom_champ_C_vap_;		    // champ source couple
  Nom nom_champ_D_vap_;			// diffusivity vap dans Nafion
  Nom nom_champ_C_N2_;		    // champ source couple
  Nom nom_champ_D_N2_;			// diffusivity N2 dans Nafion
  Nom nom_champ_T_;
  Nom nom_dom_;					// domaine
  Nom nom_ssz_CLa_;				// sous_zone
  Nom nom_ssz_CLc_;				// sous_zone

  Champ_Don gamma_CL_;			// P0
  Champ_Don por_naf_;			// P0
  Champ_Don eps_naf_;			// P0

  REF(Domaine) dom_;			// domaine
  REF(Sous_Zone) CL_a_;			// sous_zone CL anode -> X2 = H2
  REF(Sous_Zone) CL_c_;			// sous_zone CL cathode -> X2 = O2
  REF(Champ_base) ch_ci_;		// concentration multi especes { X2 vap N2 } P1NC 3 composants
  REF(Champ_base) ch_C_X2_;		// P1NC
  REF(Champ_base) ch_D_X2_;		// P0
  REF(Champ_base) ch_C_vap_;    // P1NC
  REF(Champ_base) ch_D_vap_;	// P0
  REF(Champ_base) ch_C_N2_;     // P1NC
  REF(Champ_base) ch_D_N2_;		// P0
  REF(Champ_base) ch_T_;		// P1NC

  DoubleTab ci_;		// P0 3 composants { c_X2, c_vap, c_N2}
  DoubleTab D_X2;		// P0
  DoubleTab C_X2;		// P0
  DoubleTab D_vap;		// P0
  DoubleTab C_vap;		// P0
  DoubleTab D_N2;		// P0
  DoubleTab C_N2;		// P0
  DoubleTab T_;			// P0

  DoubleTab S_;			        // P0, 3 composants S_X2, S_vap, S_N2
  DoubleTab dSdc_;				// P0, 3 composants dSdc_X2, dSdc_vap, dSdc_N2

  inline double f_Psat(double T) const;
  inline double f_lambda(double a) const;
  inline double f_derivee_lambda(double a) const;
  inline double f_Henry_H2(double T) const;
  inline double f_Henry_O2(double T) const;
  inline double f_Henry_N2(double T) const;

  double f_S_H2(double diffu, double ci, double C, double T, double por, double eps, double gamma) const;
  double f_dSdc_H2(double diffu, double ci, double C, double T, double por, double eps, double gamma) const;
  double f_S_O2(double diffu, double ci, double C, double T, double por, double eps, double gamma) const;
  double f_dSdc_O2(double diffu, double ci, double C, double T, double por, double eps, double gamma) const;
  double f_S_N2(double diffu, double ci, double C, double T, double por, double eps, double gamma) const;
  double f_dSdc_N2(double diffu, double ci, double C, double T, double por, double eps, double gamma) const;
  double f_S_vap(double diffu, double ci, double C, double T, double por, double eps, double gamma) const;
  double f_dSdc_vap(double diffu, double ci, double C, double T, double por, double eps, double gamma) const;
};

const double R = 8.314;
const double F = 96500;
const double C_SO3 = 2036.;

inline double Young_Todd_Source_multi::f_Psat(double T) const
{
  double dT = T - 46.13;
  double ds = 23.1961 - 3816.44/dT;
  double Psat = exp(ds);
  return Psat;
  //return exp(23.1961 - 3816.44/(T-46.13));
}

inline double Young_Todd_Source_multi::f_lambda(double a) const
{
  return 0.043+17.81*a-39.85*a*a+36*a*a*a;
}

inline double Young_Todd_Source_multi::f_derivee_lambda(double a) const
{
  return 17.81-39.85*2.*a+36*3.*a*a;
}

inline double Young_Todd_Source_multi::f_Henry_H2(double T) const
{
  return (1./1.01325e5)*exp(9.05e3/(R*T));
}

inline double Young_Todd_Source_multi::f_Henry_O2(double T)  const
{
  return (1/1.01325e5)*exp(5.88e3/(R*T));
}

inline double Young_Todd_Source_multi::f_Henry_N2(double T) const
{
  return 6.4e-6*exp(1300*(1/T-1/298.15));
}

#endif /* Young_Todd_Source_multi_included */
