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
// Cette classe represente le terme source dans le transport ionique dans
// la couche active et la membrane de pemfc S = jr + jp
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

    double eval_erev(double T, double a_H2, double a_O2, double a_vap, double a_ionH=1);
    double eval_eta(double psi, double phi, double erev);
    double eval_jo(double T, double a_H2, double a_O2, double a_vap, double a_ionH=1);
    double eval_jr(double jo, double eta, double T);

    void remplir_volumes();

	REF(Zone_VEF) la_zone_VEF;
	REF(Zone_Cl_VEF) la_zcl_VEF;
	DoubleVect volumes_;

	Nom nom_domaine_;			// pas necessaire car le terme source est associe a l'equation
	Nom nom_ssz_CLc_;			// le sous zone situant le terme source
	Nom nom_ssz_CLa_;			// le sous zone situant le terme source
	Nom nom_pb_psi_;			// transport electrique
	Nom nom_champ_psi_;			// champ potentiel electrique
	Nom nom_pb_T_;				// conduction de la chaleur
	Nom nom_champ_T_;			// champ temperature
	Nom nom_pb_dissolve_H2_;	// dissolve H2
	Nom nom_champ_a_H2_;		// champ activite H2 a_H2 = C_H2 / (Henry_H2 * Pref)
	Nom nom_pb_dissolve_O2_;	// dissolve 02
	Nom nom_champ_a_O2_;		// champ activite 02 a_02 = C_02 / (Henry_O2 * Pref)
	Nom nom_pb_dissolve_H20_;	// dissolve H2O
	Nom nom_champ_a_H20_;		// champ activite H2O a_H2O = lambda / lambda_eq = C_H2O*R*T / P_sat

	REF(Domaine) dom_;
	REF(Sous_Zone) CL_a_;
	REF(Sous_Zone) CL_c_;
	REF(Champ_base) ch_T_;
	REF(Champ_base) ch_a_H2_;
	REF(Champ_base) ch_a_O2_;
	REF(Champ_base) ch_a_vap_;
	REF(Champ_base) ch_phi_;
	REF(Champ_base) ch_psi_;

	DoubleTab T_;			// P0 pour VDF, P1NC pour VEF
	DoubleTab a_H2_;		// P0 pour VDF, P1NC pour VEF
	DoubleTab a_O2_;		// P0 pour VDF, P1NC pour VEF
	DoubleTab a_vap_;		// P0 pour VDF, P1NC pour VEF
	DoubleTab psi_;			// P0 pour VDF, P1NC pour VEF
	DoubleTab& phi_;		// P0 pour VDF, P1NC pour VEF		// ionic potential

	REF(Champ_Don) alpha;

	Champ_Fonc ch_Erev_;	// P0 pour VDF, P1NC pour VEF reference potential
	Champ_Fonc ch_eta_;		// P0 pour VDF, P1NC pour VEF eta
	Champ_Fonc ch_jo_;		// P0 pour VDF, P1NC pour VEF exchange current density
	Champ_Fonc ch_jr_;		// P0 pour VDF, P1NC pour VEF champ electrochemique courant jr = jo.gamma_CL.[exp()-exp()]
	//Champ_Fonc ch_jp_;  	// P0 pour VDF, P1NC pour VEF champ electrochemique courant jp = N_H2.2.F/e_CL
};

#endif /* Source_Term_pemfc_transport_ionique_included */
