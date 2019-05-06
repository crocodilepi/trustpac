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
// File      : Source_Term_Diffusion_Multiespeces.h
// Directory : $PEMFC_ROOT/src/Sources
//
/////////////////////////////////////////////////////////////////////////////

#ifndef Source_Term_Diffusion_Multiespeces_included
#define Source_Term_Diffusion_Multiespeces_included

#include <Source_base.h>
#include <Ref_Zone_VF.h>
#include <Ref_Domaine.h>
#include <Ref_Sous_Zone.h>
#include <Param.h>
#include <DoubleTab.h>

/////////////////////////////////////////////////////////////////////////////
//
// .DESCRIPTION : class Source_Term_Diffusion_Multiespeces
// Cette classe represente un term source couple:
// Source = - Si
//   avec Si champ de taux de gas dissous dans Nafion:
//	Si = Da_i*gamma_CL/eps_naf()Ceq_i - C_i)
// <Description of class Source_Term_Diffusion_Multiespeces>
//
/////////////////////////////////////////////////////////////////////////////

class Source_Term_Diffusion_Multiespeces : public Source_base
{

  Declare_instanciable( Source_Term_Diffusion_Multiespeces ) ;

public :
  DoubleTab& ajouter(DoubleTab& ) const;
  DoubleTab& calculer(DoubleTab& ) const;
  void mettre_a_jour(double temps);
  void associer_zones(const Zone_dis& ,const Zone_Cl_dis& );
  void associer_pb(const Probleme_base& );
  void completer();
  void set_param(Param& param);
protected :
  REF(Zone_VF) la_zone_;

  Nom nom_pb_X2_;			// pb de dissous X2 = O2/H2 dans Nafion
  Nom nom_pb_H2O_;			// pb de dissous X2 = O2/H2 dans Nafion
  Nom nom_pb_N2_;			// pb de dissous X2 = O2/H2 dans Nafion
  Nom nom_champ_S_X2_;		    // champ source couple
  Nom nom_champ_S_H2O_;		    // champ source couple
  Nom nom_champ_S_N2_;		    // champ source couple

  REF(Domaine) dom_;			// domaine
  REF(Champ_base) ch_S_X2_;		// champ source couple P0
  REF(Champ_base) ch_S_H2O_;	// champ source couple P0
  REF(Champ_base) ch_S_N2_;		// champ source couple P0

  DoubleTab S_;			// 3 composants, stockage de valeurs interpoles depuis ch_S_X2, ch_S_vap, ch_S_N2_

};

#endif /* Source_Term_Diffusion_Multiespeces_included */
