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
#include <Probleme_base.h>
#include <Equation_base.h>
#include <Interprete.h>
#include <Domaine.h>
#include <Zone_VF.h>

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
  CL_a_ = dom_.valeur().ss_zone(nom_ssz_CLa_);
  CL_c_ = dom_.valeur().ss_zone(nom_ssz_CLc_);
  dom_.valeur().creer_tableau_elements(ir_);
  dom_.valeur().creer_tableau_elements(ip_);
  return is;
}

void Source_Term_pemfc_transport_ionique::set_param(Param& param)
{
  //param.ajouter("nom_domaine", &nom_domaine_, Param::REQUIRED);		// DVQ: pas necessaire car associe a l'equation -> get domain
  param.ajouter("nom_ssz_CLc", &nom_ssz_CLc_, Param::REQUIRED);
  param.ajouter("nom_ssz_CLa", &nom_ssz_CLa_, Param::REQUIRED);
  param.ajouter("nom_pb_phi", &nom_pb_phi_, Param::REQUIRED);
  param.ajouter("nom_champ_ir", &nom_champ_ir_, Param::REQUIRED);
  param.ajouter("nom_champ_ip", &nom_champ_ip_, Param::REQUIRED);
  param.ajouter("sign", &sign_, Param::REQUIRED);
  param.ajouter("Cdl", &ch_cdl_, Param::REQUIRED);
}

void Source_Term_pemfc_transport_ionique::associer_pb(const Probleme_base& pb)
{
  // recuperer le milieu_base -> nothing to do
  Cerr << " Source_Term_pemfc_transport_ionique::associer_pb " << finl ;
  assert(pb.que_suis_je() == "Pb_Conduction");
}

void Source_Term_pemfc_transport_ionique::completer()
{
  // get all references to the coupling fields
  Source_base::completer();
  Probleme_base& pb_phi = ref_cast(Probleme_base,interprete().objet(nom_pb_phi_));
  ch_ir_ = pb_phi.get_champ(nom_champ_ir_);
  assert(ch_ir_.valeur().que_suis_je().find("P0") !=-1);
  ch_ip_ = pb_phi.get_champ(nom_champ_ip_);
  assert(ch_ip_.valeur().que_suis_je().find("P0") !=-1);
  //ch_cdl_ = pb_phi.get_champ(nom_champ_cdl_);
  //assert(ch_cdl_.valeur().que_suis_je().find("P0") !=-1);
}

void Source_Term_pemfc_transport_ionique::associer_zones(
  const Zone_dis& zone_dis, const Zone_Cl_dis& zcl_dis)
{
  la_zone_ = ref_cast(Zone_VF,zone_dis.valeur());
}

void Source_Term_pemfc_transport_ionique::mettre_a_jour(double temps)
{
  const DoubleTab& xp=la_zone_.valeur().xp(); // Recuperation des centre de gravite des elements pour P0

  // mettre a jour 4 tableaux de valeurs du champ couple
  ch_ir_.valeur().mettre_a_jour(temps);
  ch_ip_.valeur().mettre_a_jour(temps);
  //ch_cdl_.valeur().mettre_a_jour(temps);
  // interpolation vers P0
  ch_ir_.valeur().valeur_aux( xp, ir_ );			// ir
  ch_ip_.valeur().valeur_aux( xp, ip_ );		    // ip
  //ch_cdl_.valeur().valeur_aux( xp, cdl_ );		// ip
}


// source = ie = ir+ip si ionique, source = -ie = -(ir+ip) si electrique
DoubleTab& Source_Term_pemfc_transport_ionique::ajouter(DoubleTab& resu) const
{
  assert(resu.dimension(0)==la_zone_.valeur().nb_faces());

  DoubleTab Cdl = ch_cdl_.valeurs();
  assert(Cdl.size() == la_zone_.valeur().nb_elem());

  DoubleVect vol = la_zone_.valeur().volumes();
  for (int poly = 0; poly < CL_a_.valeur().nb_elem_tot(); ++poly)
    {
      int elem = CL_a_.valeur()(poly);
      int nb_face_elem = la_zone_.valeur().zone().nb_faces_elem(0);
      for (int f = 0; f < nb_face_elem; ++f)
        {
          int face = la_zone_.valeur().elem_faces(elem, f);
          resu(face) += sign_ * (ir_(elem) + ip_(elem))/Cdl(elem,0) * vol(elem) / nb_face_elem;
        }
    }
  for (int poly = 0; poly < CL_c_.valeur().nb_elem_tot(); ++poly)
    {
      int elem = CL_c_.valeur()(poly);
      int nb_face_elem = la_zone_.valeur().zone().nb_faces_elem(0);
      for (int f = 0; f < nb_face_elem; ++f)
        {
          int face = la_zone_.valeur().elem_faces(elem, f);
          resu(face) += sign_ * (ir_(elem) + ip_(elem))/Cdl(elem,0) * vol(elem) / nb_face_elem;
        }
    }
  return resu;
}

DoubleTab& Source_Term_pemfc_transport_ionique::calculer(DoubleTab& resu) const
{
  resu = 0;
  ajouter(resu);
  return resu;
}
