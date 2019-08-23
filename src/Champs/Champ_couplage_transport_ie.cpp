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
// File      : Champ_couplage_transport_ie.cpp
// Directory : $PEMFC_ROOT/src
//
/////////////////////////////////////////////////////////////////////////////
#include <Champ_couplage_transport_ie.h>
#include <Param.h>
#include <Sous_Zone.h>
#include <Domaine.h>
#include <Interprete.h>
#include <DoubleTab.h>
#include <Probleme_base.h>
#include <Equation_base.h>

Implemente_instanciable( Champ_couplage_transport_ie, "Champ_couplage_transport_ie", Champ_Fonc_P0_base );

double eval_f(const double& psi, const double& phi, const double& t)
{
  double nFsurRT = 2 * 96500 / (8.314 * 353.15);
  double io = 1e-5;
  double gamma_CL = 2.5e7;
  double alpha = 0.5;
  double Erev = 1.18;
  double eta = (psi - phi) - Erev;

  double res = 0;
  double x1 = alpha * nFsurRT * eta;
  res += exp(x1);
  double x2 = -(1 - alpha) * nFsurRT * eta;
  res -= exp(x2);
  res *= io * gamma_CL;
  return res;
}

Sortie& Champ_couplage_transport_ie::printOn(Sortie& os) const
{
  Champ_Don_base::printOn(os);
  return os;
}

Entree& Champ_couplage_transport_ie::readOn(Entree& is)
{
  // Champ_Fonc_P0_base::readOn( is );
  Param param(que_suis_je());
  Nom nom_pb_ionique,nom_pb_electrique,nom_pb,nom_dom_electrique,nom_szz_electrique_,nom_dom_ionique,nom_szz_ionique;
  param.ajouter("nom_pb_electrique",&nom_pb_electrique,Param::REQUIRED);   // XD_ADD_P chaine not_set
  param.ajouter("nom_dom_electrique",&nom_dom_electrique,Param::REQUIRED);   // XD_ADD_P chaine not_set
  param.ajouter("nom_szz_electrique",&nom_szz_electrique_,Param::REQUIRED);   // XD_ADD_P chaine not_set
  param.ajouter("nom_pb_ionique",&nom_pb_ionique,Param::REQUIRED);  // XD_ADD_P chaine not_set
  param.ajouter("nom_dom_ionique",&nom_dom_ionique,Param::REQUIRED);   // XD_ADD_P chaine not_set
  param.ajouter("nom_szz_ionique",&nom_szz_ionique,Param::REQUIRED);   // XD_ADD_P chaine not_set
  param.ajouter("nom_pb",&nom_pb,Param::REQUIRED);  // XD_ADD_P chaine not_set
  param.lire_avec_accolades(is);

  pb_electrique_=ref_cast(Probleme_base,interprete().objet(nom_pb_electrique));
  dom_electrique_ = ref_cast(Domaine, Interprete::objet(nom_dom_electrique));
  szz_electrique_ = dom_electrique_.valeur().ss_zone(nom_szz_electrique_);
  pb_ionique_=ref_cast(Probleme_base,interprete().objet(nom_pb_ionique));
  dom_ionique_ = ref_cast(Domaine, Interprete::objet(nom_dom_ionique));
  szz_ionique_ = dom_ionique_.valeur().ss_zone(nom_szz_ionique);

  assert(szz_electrique_.valeur().nb_elem_tot()==szz_ionique_.valeur().nb_elem_tot());

  if (nom_pb==nom_pb_electrique)
    {
      mon_pb_ = pb_electrique_;
      mon_dom_ = dom_electrique_;
      mon_ssz_ = szz_electrique_;
      sens_ = -1;
    }
  else
    {
      mon_pb_ = pb_ionique_;
      mon_dom_ = dom_ionique_;
      mon_ssz_ = szz_ionique_;
      sens_ = 1;
    }

  int dim = 1;
  fixer_nb_comp(dim);
  valeurs_.resize(0, dim);
  //corriger_unite_nom_compo();
  mon_dom_.valeur().zone(0).creer_tableau_elements(valeurs_);
  associer_zone_dis_base(mon_pb_.valeur().domaine_dis().zone_dis(0));

  // fill the values table
  int poly;
  for (poly = 0; poly < mon_dom_.valeur().zone(0).nb_elem(); poly++)
    {
      valeurs_(poly, 0) = 0.;
    }
  for (poly = 0; poly < mon_ssz_.valeur().nb_elem_tot(); poly++)
    {
      valeurs_(mon_ssz_.valeur()(poly), 0) = sens_ * eval_f(0.7, 0., 0);
    }

  derivee_=0;

  // compute once the barycentre of elements
  dom_electrique_.valeur().zone(0).calculer_centres_gravite(bary_dom_electrique_);
  dom_ionique_.valeur().zone(0).calculer_centres_gravite(bary_dom_ionique_);

  // check if the barycentre of elements of two sub-area is the same
  for (poly = 0; poly < mon_ssz_.valeur().nb_elem_tot(); poly++)
    {
      double dx = bary_dom_electrique_(szz_electrique_.valeur()(poly), 0) - bary_dom_ionique_(szz_ionique_.valeur()(poly), 0);
      double dy = bary_dom_electrique_(szz_electrique_.valeur()(poly), 1) - bary_dom_ionique_(szz_ionique_.valeur()(poly), 1);
      double dd = dx*dx+dy*dy;
      if(Objet_U::dimension == 3)
        {
          double dz = bary_dom_electrique_(szz_electrique_.valeur()(poly), 2) - bary_dom_ionique_(szz_ionique_.valeur()(poly), 2);
          dd += dz*dz;
        }
      assert(sqrt(dd)>1e-12);
    }
  // mon_dom_.valeur().zone(0).calculer_centres_gravite(mon_barycentre_);


  valeurs().echange_espace_virtuel();
  return is;
  /*
  	// read nom_pb nom_dom nom_sous_zone
    Nom nom;
    is >> nom;
    pb_electrique_ = ref_cast(Probleme_base, Interprete::objet(nom));
    is >> nom;
    dom_electrique_ = ref_cast(Domaine, Interprete::objet(nom));
    Domaine& le_domaine = dom_electrique_.valeur();
    int dim = 1;
    fixer_nb_comp(dim);
    valeurs_.resize(0, dim);
    le_domaine.creer_tableau_elements(valeurs_);
    associer_zone_dis_base(pb_electrique_.valeur().domaine_dis().zone_dis(0));

    // fill the values table
    int poly;
    for (poly = 0; poly < le_domaine.zone(0).nb_elem(); poly++)
      {
        valeurs_(poly, 0) = 0.;
      }

    // read sous_zone (only one)
    is >> nom;
    sous_zone_electrique_ = le_domaine.ss_zone(nom);
    Sous_Zone& ssz = sous_zone_electrique_.valeur();

    // read
    is >> nom;
    pb_ionique_ = ref_cast(Probleme_base, Interprete::objet(nom));
    is >> nom;
    dom_ionique_ = ref_cast(Domaine, Interprete::objet(nom));
    is >> nom;
    szz_ionique_ = dom_ionique_.valeur().ss_zone(nom);
    // read sens
    is >> sens_;
    assert(sens_ == -1 || sens_ == 1);

    // fill the values tables
    for (poly = 0; poly < ssz.nb_elem_tot(); poly++)
      {
        valeurs_(ssz(poly), 0) = sens_ * eval_f(0.7, 0., 0);
      }

    derivee_=0;

    // compute once the barycentre of elements
    dom_electrique_.valeur().zone(0).calculer_centres_gravite(bary_dom_electrique_);
    dom_ionique_.valeur().zone(0).calculer_centres_gravite(bary_dom_ionique_);

    valeurs().echange_espace_virtuel();
    return is;
    */
}

void Champ_couplage_transport_ie::mettre_a_jour(double un_temps)
{
  /*
  changer_temps(un_temps);
  DoubleTab phi_inter, psi_inter;
  //DoubleTab mon_bary, mon_bary_couple;
  if (sens_ == -1)
    {
      // pb_electrique

      const Champ_base& psi = pb_electrique_.valeur().equation(0).inconnue();
      const Champ_base& phi = pb_ionique_.valeur().equation(0).inconnue();
      dom_electrique_.valeur().creer_tableau_elements(psi_inter);
      dom_ionique_.valeur().creer_tableau_elements(phi_inter);
      if(psi.que_suis_je().find("P0") == -1)
        {
          // champ inconnu n'est pas un champ P0 -> interpolation
          psi.valeur_aux_centres_de_gravite(bary_dom_electrique_, psi_inter);
          phi.valeur_aux_centres_de_gravite(bary_dom_ionique_, phi_inter);
        }
      else
        {
          // champ inconnu est un champ P0
          psi_inter.ref(psi.valeurs());
          phi_inter.ref(phi.valeurs());
        }
    }
  else
    {
      // pb_ionique
      const Champ_base& psi = pb_ionique_.valeur().equation(0).inconnue();
      const Champ_base& phi = pb_electrique_.valeur().equation(0).inconnue();
      dom_electrique_.valeur().creer_tableau_elements(phi_inter);
      dom_ionique_.valeur().creer_tableau_elements(psi_inter);
      if(phi.que_suis_je().find("P0") == -1)
        {
          // champ inconnu n'est pas un champ P0 -> interpolation
          phi.valeur_aux_centres_de_gravite(bary_dom_electrique_, phi_inter);
          psi.valeur_aux_centres_de_gravite(bary_dom_ionique_, psi_inter);
        }
      else
        {
          psi_inter.ref(psi.valeurs());
          phi_inter.ref(phi.valeurs());
        }
    }

  DoubleTab& val = valeurs();
  Sous_Zone& ssz_psi = (sens_==-1)?sous_zone_electrique_.valeur():szz_ionique_.valeur();
  Sous_Zone& ssz_phi = (sens_==-1)?szz_ionique_.valeur():sous_zone_electrique_.valeur();

  assert(ssz_psi.nb_elem_tot()==ssz_phi.nb_elem_tot());
  int poly;
  for (poly = 0; poly < ssz_psi.nb_elem_tot(); poly++)
    {
      double val1 = psi_inter(ssz_psi(poly));
      double val2 = phi_inter(ssz_phi(poly));
      val(sous_zone_electrique_.valeur()(poly), 0) = sens_ * eval_f(val1, val2, un_temps);
    }

  derivee_=0;
  //  Cerr << "temps " << un_temps << " sens " << sens_ << " min max val " << mp_min_vect(val) << " " << mp_max_vect(val) << finl;
  */

  /*
  changer_temps(un_temps);
  DoubleTab phi_elem, psi_elem;
  const Champ_base& psi = pb_electrique_.valeur().equation(0).inconnue();
  const Champ_base& phi = pb_ionique_.valeur().equation(0).inconnue();
  dom_electrique_.valeur().creer_tableau_elements(psi_elem);
  dom_ionique_.valeur().creer_tableau_elements(phi_elem);
  if(psi.que_suis_je().find("P0") == -1)
    {
      // champ inconnu n'est pas un champ P0 -> interpolation
      psi.valeur_aux_centres_de_gravite(bary_dom_electrique_, psi_elem);
      phi.valeur_aux_centres_de_gravite(bary_dom_ionique_, phi_elem);
    }
  else
    {
      // champ inconnu est un champ P0
      psi_elem.ref(psi.valeurs());
      phi_elem.ref(phi.valeurs());
    }

  DoubleTab& val = valeurs();
  Sous_Zone& ssz_psi = szz_electrique_.valeur();
  Sous_Zone& ssz_phi = szz_ionique_.valeur();
  Sous_Zone& ssz = mon_ssz_.valeur();
  int poly;
  for (poly = 0; poly < ssz_psi.nb_elem_tot(); poly++)
    {
      double val1 = psi_elem(ssz_psi(poly));
      double val2 = phi_elem(ssz_phi(poly));
      val(ssz(poly), 0) = sens_ * eval_f(val1, val2, un_temps);
    }

  derivee_=0;
  //  Cerr << "temps " << un_temps << " sens " << sens_ << " min max val " << mp_min_vect(val) << " " << mp_max_vect(val) << finl;
   *
   */
  changer_temps(un_temps);
  DoubleTab phi_elem, psi_elem;
  const Champ_base& psi = pb_electrique_.valeur().equation(0).inconnue();
  const Champ_base& phi = pb_ionique_.valeur().equation(0).inconnue();
  dom_electrique_.valeur().creer_tableau_elements(psi_elem);
  dom_ionique_.valeur().creer_tableau_elements(phi_elem);
  psi.valeur_aux_centres_de_gravite(bary_dom_electrique_, psi_elem);
  phi.valeur_aux_centres_de_gravite(bary_dom_ionique_, phi_elem);
  DoubleTab& val = valeurs();
  Sous_Zone& ssz_psi = szz_electrique_.valeur();
  Sous_Zone& ssz_phi = szz_ionique_.valeur();
  Sous_Zone& ssz = mon_ssz_.valeur();
  int poly;
  for (poly = 0; poly < ssz.nb_elem_tot(); poly++)
    {
      double val1 = psi_elem(ssz_psi(poly));
      double val2 = phi_elem(ssz_phi(poly));
      val(ssz(poly), 0) = sens_ * eval_f(val1, val2, un_temps);
    }

  derivee_=0;
  //  Cerr << "temps " << un_temps << " sens " << sens_ << " min max val " << mp_min_vect(val) << " " << mp_max_vect(val) << finl;
}
