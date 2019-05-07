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
// File      : Source_Terme_Nafion_electro_osmotic.h
// Directory : $PEMFC_ROOT/src/Sources
//
/////////////////////////////////////////////////////////////////////////////

#ifndef Source_Terme_Nafion_electro_osmotic_included
#define Source_Terme_Nafion_electro_osmotic_included

#include <Source_base.h>
#include <Param.h>
#include <Ref_Zone_VEF.h>
#include <Ref_Zone_Cl_VEF.h>
#include <Ref_Domaine.h>
#include <Ref_Sous_Zone.h>
#include <DoubleTab.h>
#include <Ref_Champ_base.h>
#include <Ref_Champ_Inc.h>
#include <Champ_base.h>
#include <Champ_Inc.h>
#include <Champ_Don.h>

/////////////////////////////////////////////////////////////////////////////
//
// .DESCRIPTION : class Source_Terme_Nafion_electro_osmotic
//  Terme source de type electro-osmotic pour le dissous de H2O dans Nafion:
//  Source = -div(nd/F*I_i) avec I_i = -kappa.grad(phi)
//  -> champ couple: C_H2O, I_i
//  Donc la formulation integrale
//  <Source>w = nd/F*(I_i.face_normale)
// <Description of class Source_Terme_Nafion_electro_osmotic>
//
/////////////////////////////////////////////////////////////////////////////

class Source_Terme_Nafion_electro_osmotic : public Source_base
{

  Declare_instanciable( Source_Terme_Nafion_electro_osmotic ) ;

public :
  DoubleTab& ajouter(DoubleTab& ) const;
  DoubleTab& calculer(DoubleTab& ) const;
  void mettre_a_jour(double temps);
  void associer_zones(const Zone_dis& ,const Zone_Cl_dis& );
  void associer_pb(const Probleme_base& );
  void completer();
  void set_param(Param& param);
protected :
  void remplir_volumes();

  REF(Zone_VEF) la_zone_VEF;
  REF(Zone_Cl_VEF) la_zcl_VEF;
  DoubleVect volumes_;
  Nom nom_pb_phi_;				// transport ionique pour recuperer le champ electro chimique
  Nom nom_champ_I_;			    // courant ionique
  REF(Champ_base) ch_I_;
  REF(Champ_Inc) ch_C_;
  DoubleTab I_, C_;
  Champ_Don por_naf_;			// porosite de Nafion
  Champ_Don eps_naf_;			// ionomer proportionnel de Nafion
  inline double f_nd(double C) const;
};

const double C_SO3 = 2036.; 	// [mol/m^3], Concentration en sites sulfones dans le Nafion
const double F = 96500;			// Faraday constant C/mol
inline double Source_Terme_Nafion_electro_osmotic::f_nd(double C) const
{
  double ld = C / C_SO3;
  return 1. + 0.0028*ld + 0.0026*ld*ld;
}

#endif /* Source_Terme_Nafion_electro_osmotic_included */
