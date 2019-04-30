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
// File      : Loi_Fermeture_Diffusion_Nafion.h
// Directory : $PEMFC_ROOT/src/Loi_Fermeture
//
/////////////////////////////////////////////////////////////////////////////

#ifndef Loi_Fermeture_Diffusion_Nafion_included
#define Loi_Fermeture_Diffusion_Nafion_included

#include <Loi_Fermeture_base.h>
#include <Ref_Champ_base.h>
#include <Ref_Champ_Inc.h>
#include <Champ_Fonc.h>
#include <Ref_Equation_base.h>
#include <DoubleTab.h>
#include <Param.h>

/////////////////////////////////////////////////////////////////////////////
//
// .DESCRIPTION : class Loi_Fermeture_Diffusion_Nafion
// Cette classe est couple avec un probleme de diffusion de l'espece dans Nafion
// (1-por)eps_naf.dC/dt = -div(D_i_eff.grad(C)) + Source_term_pemfc_.
// Loi_Fermeture_Diffusion_Nafion permet de creer et mettre a jour
// + champ de coefficient de diffusion 'effective' des especes dans Nafion
//   D_i_eff = (1-por)epsilon_ionomer/tor^2*D_i_naf(T)
//   avec i = H2, O2, N2
//   ou
//   D_i_eff = -D_w avec H20
// + champ de diffusion des especes dans Nafion i = H2, O2, N2
//   D_i_naf(T)
//   ou bien le champ de diffusion de H20 dans Nafion
//   Dw(T, C_H20)
// + champ de flux de diffusion des especes dans Nafion i = H2, O2, N2
//   N_i_naf = -D_i_eff * grad(C_i)
//   ou bien champ de flux de diffusion de H20 dans Nafion
//   N_H20_naf = nd/F*I_ion - Dw*grad(C_H20)
//   avec I_ion est le champ de courant ionique
// <Description of class Loi_Fermeture_Diffusion_Nafion>
//
/////////////////////////////////////////////////////////////////////////////

class Loi_Fermeture_Diffusion_Nafion : public Loi_Fermeture_base
{

  Declare_instanciable( Loi_Fermeture_Diffusion_Nafion ) ;

public :
  // associer avec le probleme de diffusion
  void associer_pb_base(const Probleme_base&);
  // discretiser le champ diffusivite
  void discretiser(const Discretisation_base& );
  // recuperer tous Reference
  void completer();
  // mettre a jour une fois avant du calcul iteratif
  void preparer_calcul();
  // mettre a jour le tableau de valeurs
  void mettre_a_jour(double temps);
  //
  void set_param(Param& param);
protected :
  Champ_Fonc diffu_;		// champ de diffusion 'effective' de l'espece
  Champ_Fonc D_i_naf_;		// champ de diffusion de l'espece
  Champ_Fonc N_i_naf_;		// champ de flux de diffusion de l'espece

  REF(Equation_base) ref_equation_;
  const Equation_base& equation() const
  {
    return ref_equation_.valeur();
  } ;

  // Motcles nom_especes_compris_ = {"H2", "O2", "H2O", "N2", "vap"};

  Nom nom_espece_, nom_pb_T_, nom_champ_T_, nom_pb_phi_, nom_champ_I_;	// pour readOn
  double T_0_;				// dans le cas T constant
  double CSO3_;				// au cas ou "H2O" ou "vap"
  double por_naf_;			// porosite de Nafion
  double eps_naf_;			// ionomer proportionnel de Nafion
  double tor_naf_;			// tortuosite de Nafion

  REF(Champ_base) ch_T_;		// champ reference pour temperature (optionnel)
  REF(Champ_base) ch_I_;		// champ reference pour le courant ionique I = -kappa.grad(phi) (optionnel)
  REF(Champ_Inc) ch_C_;			// champ reference du champ inconnu C

  DoubleTab T_, C_, I_;			// tableau des valeurs du champ T, C, I (P0)
  double eval_D_i_naf(double T, double C);
  double eval_diffu_(double T, double C);
};

#endif /* Loi_Fermeture_Diffusion_Nafion_included */
