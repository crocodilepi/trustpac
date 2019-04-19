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
// File      : Source_Term_pemfc_transport_ionique.cpp
// Directory : $PEMFC_ROOT/src/Sources
//
/////////////////////////////////////////////////////////////////////////////

#include <Source_Term_pemfc_transport_ionique.h>
#include <Zone_VEF.h>
#include <Zone_Cl_VEF.h>
#include <Probleme_base.h>

Implemente_instanciable( Source_Term_pemfc_transport_ionique, "Source_Term_pemfc_transport_ionique_VEF_P1NC", Source_base ) ;

Sortie& Source_Term_pemfc_transport_ionique::printOn( Sortie& os ) const
{
  Source_base::printOn( os );
  return os;
}

Entree& Source_Term_pemfc_transport_ionique::readOn( Entree& is )
{
	Source_base::readOn( is );
  	Cerr << " Source_Term_pemfc_transport_ionique::readOn " << finl  ;
	Param param(que_suis_je());
	set_param(param);
	param.lire_avec_accolades(is);

	dom_ = equation().probleme().domaine();
	ssz_ = dom_.valeur().ss_zone(nom_ssz_);
	Probleme_base& pb_T = ref_cast(Probleme_base,interprete().objet(nom_pb_T_));
	ch_T_ = pb_T.get_champ(nom_champ_T_);
	assert(ch_T_.que_suis_je().find("P1NC") !=-1);

	Probleme_base& pb_c = ref_cast(Probleme_base,interprete().objet(nom_pb_c_));
	ch_a_ = pb_c.get_champ(nom_champ_a_);
	assert(ch_a_.que_suis_je().find("P1NC") !=-1);

	Probleme_base& pb_psi = ref_cast(Probleme_base,interprete().objet(nom_pb_psi_));
	ch_psi_ = pb_psi.get_champ(nom_champ_psi_);
	assert(ch_psi_.que_suis_je().find("P1NC") !=-1);

	ch_phi_ = equation().inconnue();
	assert(ch_phi_.que_suis_je().find("P1NC") !=-1);
	phi_.ref(ch_phi_.valeur().valeurs());

	int dim = 1;		// all fields are scalar
	T_.resize(0, dim);
	a_.resize(0, dim);
	psi_.resize(0, dim);
	la_zone_VEF.valeur().creer_tableau_faces(T_);
	la_zone_VEF.valeur().creer_tableau_faces(a_);
	la_zone_VEF.valeur().creer_tableau_faces(psi_);

  return is;
}

void Source_Term_pemfc_transport_ionique::set_param(Param& param) {
	//param.ajouter("nom_domaine", &nom_domaine_, Param::REQUIRED);		// quang: pas necessaire car associer(equation)
	param.ajouter("nom_sous_zone", &nom_ssz_, Param::REQUIRED);
	param.ajouter("nom_pb_psi", &nom_pb_psi_, Param::REQUIRED);
	param.ajouter("nom_champ_psi", &nom_champ_psi_, Param::REQUIRED);
	param.ajouter("nom_pb_T", &nom_pb_T_, Param::REQUIRED);
	param.ajouter("nom_champ_T", &nom_champ_T_, Param::REQUIRED);
	param.ajouter("nom_pb_c", &nom_pb_c_, Param::REQUIRED);
	param.ajouter("nom_champ_c", &nom_champ_a_, Param::REQUIRED);
}

void Source_Term_pemfc_transport_ionique::associer_pb(const Probleme_base& pb) {
	// recuperer le milieu_base -> nothing to do
	Cerr << " Source_Term_pemfc_transport_ionique::associer_pb " << finl ;
	  assert(pb.que_suis_je() == "Pb_Conduction");
	  int ok = 0;
	  const Equation_base& eqn = pb.equation(0);
	  if  (eqn.que_suis_je() == "Conduction")
	    {
	      associer_zones(eqn.zone_dis(),eqn.zone_Cl_dis());
	      ok = 1;
	    }
	  if (!ok)
	    {
	      Cerr << "Erreur TRUST dans Source_Term_pemfc_transport_ionique::associer_pb()" << finl;
	      Cerr << "On ne trouve pas d'equation de conduction dans le probleme" << finl;
	      exit();
	    }

}

void Source_Term_pemfc_transport_ionique::completer() {
	// get all references to the coupling fields
	Source_base::completer();
}

void Source_Term_pemfc_transport_ionique::associer_zones(
		const Zone_dis& zone_dis, const Zone_Cl_dis& zcl_dis) {
	la_zone_VEF = ref_cast(Zone_VEF,zone_dis.valeur());
	  la_zcl_VEF = ref_cast(Zone_Cl_VEF,zcl_dis.valeur());
	  remplir_volumes();
}

void Source_Term_pemfc_transport_ionique::discretiser(const Discretisation_base& dis)
{
	  dis.discretiser_champ("temperature",equation().zone_dis().valeur(),"Erev","unit", 1 , 0. , ch_Erev_);
	  champs_compris_.ajoute_champ(ch_Erev_);
	  dis.discretiser_champ("temperature",equation().zone_dis().valeur(),"eta","unit", 1 , 0. , ch_eta_);
	  champs_compris_.ajoute_champ(ch_eta_);
	  dis.discretiser_champ("temperature",equation().zone_dis().valeur(),"jo","unit", 1 , 0. , ch_jo_);
	  champs_compris_.ajoute_champ(ch_jo_);
	  dis.discretiser_champ("temperature",equation().zone_dis().valeur(),"jr","unit", 1 , 0. , ch_jr_);
	  champs_compris_.ajoute_champ(ch_jr_);
	  dis.discretiser_champ("temperature",equation().zone_dis().valeur(),"jp","unit", 1 , 0. , ch_jp_);
	  champs_compris_.ajoute_champ(ch_jp_);

}

void Source_Term_pemfc_transport_ionique::mettre_a_jour(double temps) {
	// mettre a jour 4 tableaux de valeurs du champ couple
	ch_T_.valeur().mettre_a_jour(temps);
	ch_a_.valeur().mettre_a_jour(temps);
	ch_psi_.valeur().mettre_a_jour(temps);
	ch_phi_.valeur().mettre_a_jour(temps);

	const DoubleTab& xv=la_zone_VEF.valeur().xv(); // centre de gravite des faces pour P1NC

	// interpolation vers P1NC
	ch_T_.valeur().valeur_aux( xv, T_ );
	ch_a_.valeur().valeur_aux( xv, a_ );
	ch_psi_.valeur().valeur_aux( xv, psi_ );
	phi_.ref(ch_phi_.valeur().valeurs());

	// mettre a jour des champs compris
	DoubleTab& val_eta = ch_eta_.valeurs();
	DoubleTab& val_erev_ = ch_Erev_.valeurs();
	DoubleTab& val_jo = ch_jo_.valeurs();
	DoubleTab& val_jr = ch_jr_.valeurs();
	DoubleTab& val_jo = ch_jp_.valeurs();

	IntTab faces_ssz;	// faces belong to the sous_zone -> flag = 1, if not, flag = 0
	la_zone_VEF.valeur().creer_tableau_faces(faces_ssz);
	faces_ssz = 0;		// init with no flag (all faces are unchecked)

	for (int poly = 0; poly < ssz_.valeur().nb_elem_tot(); poly++)
	{
	  int elem = ssz_.valeur()(poly);
	  for (int f = 0; f < la_zone_VEF.valeur().zone().nb_faces_elem(0); ++f)
		{
		  int face = la_zone_VEF.valeur().elem_faces(elem, f);
		  if(!faces_ssz(face))
			{
			  val_eta(face) = psi_(face) - phi_(face) - val_erev_(face);
			  // necessaire (source*porosite_surf(num_face));
			  faces_ssz(face) = 1;		// marquer comme deja traite
			}
		}
	}
}

void Source_Term_pemfc_transport_ionique::eval_erev(double T, double a_H2, double a_O2, double a_vap, double a_ionH = 1)
{

}

void Source_Term_pemfc_transport_ionique::remplir_volumes() {
	volumes_.ref(la_zone_VEF.valeur().volumes_entrelaces());
}
