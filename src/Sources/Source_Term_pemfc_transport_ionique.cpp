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
	CL_a_ = dom_.valeur().ss_zone(nom_ssz_CLc_);
	CL_c_ = dom_.valeur().ss_zone(nom_ssz_CLc_);
	Probleme_base& pb_T = ref_cast(Probleme_base,interprete().objet(nom_pb_T_));
	ch_T_ = pb_T.get_champ(nom_champ_T_);
	assert(ch_T_.que_suis_je().find("P1NC") !=-1);

	Probleme_base& pb_H2 = ref_cast(Probleme_base,interprete().objet(nom_pb_dissolve_H2_));
	ch_a_H2_ = pb_H2.get_champ(nom_champ_a_H2_);
	assert(ch_a_H2_.que_suis_je().find("P1NC") !=-1);

	Probleme_base& pb_O2 = ref_cast(Probleme_base,interprete().objet(nom_pb_dissolve_O2_));
	ch_a_O2_ = pb_O2.get_champ(nom_champ_a_O2_);
	assert(ch_a_O2_.que_suis_je().find("P1NC") !=-1);

	Probleme_base& pb_H2O = ref_cast(Probleme_base,interprete().objet(nom_pb_dissolve_H20_));
	ch_a_vap_ = pb_H2O.get_champ(nom_champ_a_H20_);
	assert(ch_a_vap_.que_suis_je().find("P1NC") !=-1);

	Probleme_base& pb_psi = ref_cast(Probleme_base,interprete().objet(nom_pb_psi_));
	ch_psi_ = pb_psi.get_champ(nom_champ_psi_);
	assert(ch_psi_.que_suis_je().find("P1NC") !=-1);

	ch_phi_ = equation().inconnue();
	assert(ch_phi_.que_suis_je().find("P1NC") !=-1);
	phi_.ref(ch_phi_.valeur().valeurs());

	int dim = 1;		// all fields are scalar
	T_.resize(0, dim);
	a_H2_.resize(0, dim);
	a_O2_.resize(0, dim);
	a_vap_.resize(0, dim);
	psi_.resize(0, dim);
	la_zone_VEF.valeur().creer_tableau_faces(T_);
	la_zone_VEF.valeur().creer_tableau_faces(a_H2_);
	la_zone_VEF.valeur().creer_tableau_faces(a_O2_);
	la_zone_VEF.valeur().creer_tableau_faces(a_vap_);
	la_zone_VEF.valeur().creer_tableau_faces(psi_);

  return is;
}

void Source_Term_pemfc_transport_ionique::set_param(Param& param) {
	//param.ajouter("nom_domaine", &nom_domaine_, Param::REQUIRED);		// quang: pas necessaire car associer(equation)
	param.ajouter("nom_ssz_CLc", &nom_ssz_CLc_, Param::REQUIRED);
	param.ajouter("nom_ssz_CLa", &nom_ssz_CLa_, Param::REQUIRED);
	param.ajouter("nom_pb_psi", &nom_pb_psi_, Param::REQUIRED);
	param.ajouter("nom_champ_psi", &nom_champ_psi_, Param::REQUIRED);
	param.ajouter("nom_pb_T", &nom_pb_T_, Param::REQUIRED);
	param.ajouter("nom_champ_T", &nom_champ_T_, Param::REQUIRED);
	param.ajouter("nom_pb_dissolve_H2", &nom_pb_dissolve_H2_, Param::REQUIRED);
	param.ajouter("nom_champ_a_H2", &nom_champ_a_H2_, Param::REQUIRED);
	param.ajouter("nom_pb_dissolve_O2", &nom_pb_dissolve_O2_, Param::REQUIRED);
	param.ajouter("nom_champ_a_O2", &nom_champ_a_O2_, Param::REQUIRED);
	param.ajouter("nom_pb_dissolve_H2O", &nom_pb_dissolve_H20_, Param::REQUIRED);
	param.ajouter("nom_champ_a_H2O", &nom_champ_a_H20_, Param::REQUIRED);
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
	  //dis.discretiser_champ("temperature",equation().zone_dis().valeur(),"jp","unit", 1 , 0. , ch_jp_);
	  //champs_compris_.ajoute_champ(ch_jp_);

}

void Source_Term_pemfc_transport_ionique::mettre_a_jour(double temps) {
	// mettre a jour 4 tableaux de valeurs du champ couple
	ch_T_.valeur().mettre_a_jour(temps);
	ch_a_H2_.valeur().mettre_a_jour(temps);
	ch_a_O2_.valeur().mettre_a_jour(temps);
	ch_a_vap_.valeur().mettre_a_jour(temps);
	ch_psi_.valeur().mettre_a_jour(temps);
	ch_phi_.valeur().mettre_a_jour(temps);

	const DoubleTab& xv=la_zone_VEF.valeur().xv(); // centre de gravite des faces pour P1NC

	// interpolation vers P1NC
	ch_T_.valeur().valeur_aux( xv, T_ );
	ch_a_H2_.valeur().valeur_aux( xv, a_H2_ );
	ch_a_O2_.valeur().valeur_aux( xv, a_O2_ );
	ch_a_vap_.valeur().valeur_aux( xv, a_vap_ );
	ch_psi_.valeur().valeur_aux( xv, psi_ );
	phi_.ref(ch_phi_.valeur().valeurs());

	// mettre a jour des champs compris
	DoubleTab& val_eta = ch_eta_.valeurs();
	DoubleTab& val_erev_ = ch_Erev_.valeurs();
	DoubleTab& val_jo = ch_jo_.valeurs();
	DoubleTab& val_jr = ch_jr_.valeurs();
	//DoubleTab& val_jp = ch_jp_.valeurs();

	IntTab faces_ssz;	// faces belong to the sous_zone -> flag = 1, if not, flag = 0
	la_zone_VEF.valeur().creer_tableau_faces(faces_ssz);
	faces_ssz = 0;		// init with no flag (all faces are unchecked)

	for (int poly = 0; poly < CL_a_.valeur().nb_elem_tot(); poly++)
	{
	  int elem = CL_a_.valeur()(poly);
	  for (int f = 0; f < la_zone_VEF.valeur().zone().nb_faces_elem(0); ++f)
		{
		  int face = la_zone_VEF.valeur().elem_faces(elem, f);
		  if(!faces_ssz(face))
			{
			  double erev = eval_erev(T_(face),a_O2_(face), a_vap_(face), 1);
			  double eta  = eval_eta(psi_(face), phi_(face), erev);
			  double jo   = eval_jo(T_(face), a_H2_(face), a_O2_(face), a_vap_(face), 1);
			  val_erev_(face) 	= erev;
			  val_eta(face) 	= eta;
			  val_jo(face)		= jo;
			  val_jr(face) 		= eval_jr(jo, eta, T_(face));
			  // necessaire (source*porosite_surf(num_face));
			  faces_ssz(face) = 1;		// marquer comme deja traite
			}
		}
	}
}

double Source_Term_pemfc_transport_ionique::eval_erev(double T, double a_H2, double a_O2, double a_vap, double a_ionH = 1)
{
	double nF = 2. * 96500.;
	double RTsurnF = 8.314 * T / nF;
	double RTsurnFlog = RTsurnF * a_H2 * sqrt(a_O2) / a_vap; // pow(a_H2, 1) * pow(a_O2, 0.5) * pow(a_vap, -1) * pow(a_ionH,-2)
	double delta_H = 25e3; 			// TO-DO: dHox0_a   = 25e3; dHox0_c   = 167.9e3;
	double delta_S = -172; 			// TO-DO: dSox0_a   = -172; dSox0_c   = -205.6;
	double delta_G = delta_H - T * delta_S;
	return -delta_G/nF + RTsurnFlog;
}

double Source_Term_pemfc_transport_ionique::eval_eta(double psi, double phi, double erev)
{
	return psi - phi - erev;
}



/************************ ATTENTION: WRONG FORMULE: Anode # Cathode **************************/
double Source_Term_pemfc_transport_ionique::eval_jo(double T, double a_H2, double a_O2, double a_vap, double a_ionH = 1)
{
	double kB = 1.38064852e-23;	 	// Boltzman constant [J/K]
	double s0 = 6.41e-20;			// average Pt surface per reaction site [m2]
	double NA = 6.022140857e23;		// Avogadro constant [/mol]
	double h  = 6.62607015e-34;		// Planck constant [Js]
	double k0 = kB*T/(s0*NA*h);
	double RT = 8.314 * T;
	double delta_Hox = 24.36e3;		// [J/mol]	-> TO-DO: dHox0_a   = 24.36e3; dHox0_c   = 81.52e3;
	double delta_Sox = -172.3;		// [J/mol/K]-> TO-DO: dSox0_a   = -172.3; dSox0_c   = -285;
	double delta_Gox = delta_Hox - T*delta_Sox;
	double delta_H = 25e3; 			// TO-DO: dHox0_a   = 25e3; dHox0_c   = 167.9e3;
	double delta_S = -172; 			// TO-DO: dSox0_a   = -172; dSox0_c   = -205.6;
	double delta_G = delta_H - T * delta_S;
	double delta_Gred = delta_Gox + delta_G;
	double k0_ox = k0*exp(delta_Gox/RT);
	double k0_red=k0*exp(delta_Gred/RT);
	double nF = 2.*96500.;
	double alpha = 0.5;				// TO-DO: alpha_a = 0.5; alpha_c = 0.216
	double i0 = nF*pow(k0_ox, (1-alpha))*pow(k0_red,alpha);

	// VERIFIER
	// Anode i0*pow(a_H2, (1-alpha_a)*1)*pow(a_ionH, alpha_a*2)
	// Cathode i0*pow(a_O2, alpha_c*0.5)*pow(a_ionH, alpha_c*2)*pow(a_vap, (1-alpha_c)*1)
	return i0*sqrt(a_H2);
	// return i0*pow(a_O2,alpha*0.5)*pow(a_vap,1-alpha)
}

double Source_Term_pemfc_transport_ionique::eval_jr(double jo, double eta, double T)
{
	double nFsurRT = 2 * 96500 / (8.314 * T);
	double gamma_CL = 1.67e7;
	double alpha = 0.5;				// TO-DO: alpha_a = 0.5; alpha_c = 0.216
	double res = 0;
	double x1 = alpha * nFsurRT * eta;
	res += exp(x1);
	double x2 = -(1 - alpha) * nFsurRT * eta;
	res -= exp(x2);
	res *= jo * gamma_CL;
	return res;
}

void Source_Term_pemfc_transport_ionique::remplir_volumes() {
	volumes_.ref(la_zone_VEF.valeur().volumes_entrelaces());
}
