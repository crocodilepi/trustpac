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
// File      : Source_Term_pemfc_transport_ionique.h
// Directory : $PEMFC_ROOT/src/Sources
//
/////////////////////////////////////////////////////////////////////////////

#ifndef Source_Term_pemfc_transport_ionique_included
#define Source_Term_pemfc_transport_ionique_included

#include <Source_base.h>
#include <Ref_Zone_VF.h>
#include <Ref_Zone_VEF.h>
#include <Ref_Zone_Cl_VEF.h>
#include <Ref_Domaine.h>
#include <Ref_Sous_Zone.h>
#include <Param.h>
#include <DoubleTab.h>
#include <Champ_Fonc.h>
#include <Champ_Don.h>
#include <Discretisation_base.h>

/////////////////////////////////////////////////////////////////////////////
//
// .DESCRIPTION : class Source_Term_pemfc_transport_ionique
// Cette classe a pour but de calculer le terme source electro-chimique de la pile:
// Source = ie pour le transport ionique
// Source = -ie pour le transport electrique
// avec:
// + ie = ir+ip 	pour CA
// + ie = 0 		pour MB
// Les champs suivants sont associes avec:
// + E_rev 	champ de potentiel standard
// + eta 	psi - phi - E_rev
// + i0 	exchange current density (formula dans CA anode et cathode est different)
// + ir		terme source de reaction (le courant de reaction total a la surface du catalyseur, <0 a la cathode)
// + ip		terme source de permeation (le courant de permeation)
// + Q_reac Chaleur liee a la reaction
// + Q_perm Chaleur liee au courant de permeation
// <Description of class Source_Term_pemfc_transport_ionique>
//
/////////////////////////////////////////////////////////////////////////////

class Source_Term_pemfc_transport_ionique : public Source_base
{

  Declare_instanciable( Source_Term_pemfc_transport_ionique ) ;

public :
  DoubleTab& ajouter(DoubleTab& ) const;
  DoubleTab& calculer(DoubleTab& ) const;
  void contribuer_a_avec(const DoubleTab&, Matrice_Morse&) const;
  void mettre_a_jour(double temps);
  void associer_zones(const Zone_dis& ,const Zone_Cl_dis& );
  void associer_pb(const Probleme_base& );
  void completer();
  void set_param(Param& param);
protected :

  REF(Zone_VF) la_zone_;

  Nom nom_ssz_CLc_;			// le sous zone situant le terme source
  Nom nom_ssz_CLa_;			// le sous zone situant le terme source
  Nom nom_pb_phi_;			// transport ionique
  Nom nom_champ_ir_;
  Nom nom_champ_ip_;
  Nom nom_champ_cdl_;
  double sign_;
  Champ_Don ch_cdl_;
  Nom nom_champ_dir_dphi_;

  REF(Domaine) dom_;		// domaine de transport ionique (CL_c + MB + CL_a)
  REF(Sous_Zone) CL_a_;		// sous_zone anode
  REF(Sous_Zone) CL_c_;		// sous_zone cathode
  REF(Champ_base) ch_ir_;
  REF(Champ_base) ch_ip_;
  REF(Champ_base) ch_dir_dphi_;


  DoubleTab ir_, ip_, dir_dphi_;
};

#endif /* Source_Term_pemfc_transport_ionique_included */
