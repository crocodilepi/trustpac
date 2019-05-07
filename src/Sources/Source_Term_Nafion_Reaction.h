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
// File      : Source_Term_Nafion_Reaction.h
// Directory : $PEMFC_ROOT/src/Sources
//
/////////////////////////////////////////////////////////////////////////////

#ifndef Source_Term_Nafion_Reaction_included
#define Source_Term_Nafion_Reaction_included

#include <Source_base.h>
#include <Ref_Zone_VF.h>
#include <Ref_Domaine.h>
#include <Ref_Sous_Zone.h>
#include <Param.h>
#include <DoubleTab.h>
#include <Champ_Don.h>

/////////////////////////////////////////////////////////////////////////////
//
// .DESCRIPTION : class Source_Term_Nafion_Reaction
//  Source terme Si qui vient de la reaction electro-chimique dans la couche active CLa, CLc
//    S_H2 = -ir / (2F)
//    S_O2 = (ir+ip) / (4F)
//    S_N2 = 0
//    S_H20 = -(ir+ip) / (2F)
// <Description of class Source_Term_Nafion_Reaction>
//
/////////////////////////////////////////////////////////////////////////////

class Source_Term_Nafion_Reaction : public Source_base
{

  Declare_instanciable( Source_Term_Nafion_Reaction ) ;

public :
  DoubleTab& ajouter(DoubleTab& ) const;
  DoubleTab& calculer(DoubleTab& ) const;
  void mettre_a_jour(double temps);
  void associer_zones(const Zone_dis& ,const Zone_Cl_dis& );
  void associer_pb(const Probleme_base& );
  void completer();
  void set_param(Param& param);
protected :
  double eval_f(double jr, double jp) const;

  REF(Zone_VF) la_zone_;

  Nom nom_espece_;
  Nom nom_ssz_CLc_;			// le sous zone situant le terme source
  Nom nom_ssz_CLa_;			// le sous zone situant le terme source
  Nom nom_pb_phi_;			// transport ionique
  Nom nom_champ_ir_;
  Nom nom_champ_ip_;

  REF(Domaine) dom_;
  REF(Sous_Zone) CL_a_;
  REF(Sous_Zone) CL_c_;
  REF(Champ_base) ch_ir_;
  REF(Champ_base) ch_ip_;

  DoubleTab ir_, ip_;	// P0

  Champ_Don por_naf_;			// porosite de Nafion
  Champ_Don eps_naf_;			// ionomer proportionnel de Nafion

};

#endif /* Source_Term_Nafion_Reaction_included */
