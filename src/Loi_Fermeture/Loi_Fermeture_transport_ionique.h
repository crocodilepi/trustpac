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
// File      : Loi_Fermeture_transport_ionique.h
// Directory : $PEMFC_ROOT/src/Loi_Fermeture
//
/////////////////////////////////////////////////////////////////////////////

#ifndef Loi_Fermeture_transport_ionique_included
#define Loi_Fermeture_transport_ionique_included

#include <Loi_Fermeture_base.h>
#include <Param.h>
#include <Discretisation_base.h>
#include <Champ_Fonc.h>
#include <Equation_base.h>
#include <Ref_Champ_base.h>
#include <Ref_Champ_Inc.h>
#include <Ref_Operateur_base.h>
#include <DoubleTab.h>
#include <Ref_Sous_Zone.h>


/////////////////////////////////////////////////////////////////////////////
//
// .DESCRIPTION : class Loi_Fermeture_transport_ionique
// cette clase permet de creer / mettre a jour
// - champ scalaire de coefficient de conductivite inonique kappa
//   qui est dependant a la temperature T et l'absoption de l'eau C_H20:
//   kappa = f_kappa(T,ld) avec ld = C_H20 / C_SO3
// - champ de courant ionique Ii = - kappa*grad(phi)
//   avec phi est le champ inconnu de l'equation de transport ionique
//   div(-kappa.grad(phi)) = source_term_pemfc_transport_ionique
// <Description of class Loi_Fermeture_transport_ionique>
//
/////////////////////////////////////////////////////////////////////////////

class Loi_Fermeture_transport_ionique : public Loi_Fermeture_base
{

  Declare_instanciable( Loi_Fermeture_transport_ionique ) ;

public :
  // associer avec le probleme de diffusion
  //void associer_pb_base(const Probleme_base&);
  // discretiser le champ diffusivite
  void discretiser(const Discretisation_base& );
  // recuperer tous Reference
  void completer();
  // mettre a jour une fois avant du calcul iteratif
  //void preparer_calcul();
  // mettre a jour le tableau de valeurs
  void mettre_a_jour(double temps);
  //
  void set_param(Param& param);
protected :
  Champ_Fonc kappa_;	// champ de conductivite (scalaire) P0
  Champ_Fonc I_i_;		// champ de courant (vectoriel) P0

  REF(Equation_base) ref_equation_;
  REF(Sous_Zone) MB_;						// sous_zone membrane
  REF(Sous_Zone) CL_a_;						// sous_zone anode
  REF(Sous_Zone) CL_c_;						// sous_zone cathode
  const Equation_base& equation() const
  {
    return ref_equation_.valeur();
  } ;

  // pour readOn
  Nom nom_ssz_MB_;
  Nom nom_ssz_CLa_;
  Nom nom_ssz_CLc_;
  Nom nom_pb_T_;
  Nom nom_champ_T_;
  Nom nom_pb_C_;
  Nom nom_champ_C_;
  double T_0_;				// dans le cas T constant
  double C_SO3_;			//
  double por_naf_;			// porosite de Nafion
  double eps_naf_;			// ionomer proportionnel de Nafion
  double tor_naf_;			// tortuosite de Nafion

  REF(Champ_Inc) ch_phi_;   // champ reference pour potentiel ionique
  REF(Champ_base) ch_T_;	// champ reference pour temperature (optionnel)
  REF(Champ_base) ch_C_;	// champ reference pour le courant ionique I = -kappa.grad(phi) (optionnel)

  DoubleTab T_, C_, phi_;			// tableau des valeurs du champ T, C, phi (P0)
  double f_kappa(double T, double C) const;
};

inline double Loi_Fermeture_transport_ionique::f_kappa(double T, double C) const
{
  double ld = C / C_SO3_;
  double sigma = exp(1268*(1./303.-1./T))*(-0.326+0.5139*ld);
  return sigma*(1-por_naf_)*eps_naf_/(tor_naf_*tor_naf_);
}

#endif /* Loi_Fermeture_transport_ionique_included */
