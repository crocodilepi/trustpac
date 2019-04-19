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
#include <Ref_Zone_VEF.h>
#include <Ref_Zone_Cl_VEF.h>

/////////////////////////////////////////////////////////////////////////////
//
// .DESCRIPTION : class Source_Term_pemfc_transport_ionique
//
// <Description of class Source_Term_pemfc_transport_ionique>
//
/////////////////////////////////////////////////////////////////////////////

class Source_Term_pemfc_transport_ionique : public Source_base
{

  Declare_instanciable( Source_Term_pemfc_transport_ionique ) ;

public :
  	DoubleTab& ajouter(DoubleTab& ) const;
    DoubleTab& calculer(DoubleTab& ) const;
    void mettre_a_jour(double temps);
    void associer_zones(const Zone_dis& ,const Zone_Cl_dis& );
    void associer_pb(const Probleme_base& );
    void completer();
    void set_param(Param& param);
    void discretiser(const Discretisation_base& dis);
protected :

    void eval_erev(double T, double a_H2, double a_O2, double a_vap, double a_ionH = 1);
    void eval_eta(double psi, double phi, double erev);
    void eval_jo(double a_H2, double a_O2, double a_vap, double a_ionH = 1);
    void eval_jr(double jo, double eta, double T);

    void remplir_volumes();
	double eval_f(double jr, double jp) const;
	REF(Zone_VEF) la_zone_VEF;
	REF(Zone_Cl_VEF) la_zcl_VEF;
	DoubleVect volumes_;

	Nom nom_domaine_;			// pas necessaire car associe a l'equation
	Nom nom_ssz_;
	Nom nom_pb_psi_;
	Nom nom_champ_psi_;
	Nom nom_pb_T_;
	Nom nom_champ_T_;
	Nom nom_pb_c_;
	Nom nom_champ_a_;		// champ activite

	REF(Domaine) dom_;
	REF(Sous_Zone) ssz_;
	REF(Champ_base) ch_T_;
	REF(Champ_base) ch_a_;
	REF(Champ_base) ch_phi_;
	REF(Champ_base) ch_psi_;

	DoubleTab T_;			// P0 pour VDF, P1NC pour VEF
	DoubleTab a_;			// P0 pour VDF, P1NC pour VEF
	DoubleTab psi_;			// P0 pour VDF, P1NC pour VEF
	DoubleTab& phi_;			// P0 pour VDF, P1NC pour VEF		// ionic potential

	Champ_Fonc ch_Erev_;		// P0 pour VDF, P1NC pour VEF reference potential
	Champ_Fonc ch_eta_;		// P0 pour VDF, P1NC pour VEF eta
	Champ_Fonc ch_jo_;			// P0 pour VDF, P1NC pour VEF exchange current density
	Champ_Fonc ch_jr_;		// P0 pour VDF, P1NC pour VEF champ electrochemique courant jr = jo.gamma_CL.[exp()-exp()]
	Champ_Fonc ch_jp_;  	// P0 pour VDF, P1NC pour VEF champ electrochemique courant jp = N_H2.2.F/e_CL
};

#endif /* Source_Term_pemfc_transport_ionique_included */
