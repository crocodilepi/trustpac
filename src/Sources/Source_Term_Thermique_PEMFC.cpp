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
// File      : Source_Term_Thermique_PEMFC.cpp
// Directory : $PEMFC_ROOT/src/Sources
//
/////////////////////////////////////////////////////////////////////////////

#include <Source_Term_Thermique_PEMFC.h>
#include <Probleme_base.h>
#include <Equation_base.h>
#include <Matrice_Morse.h>
#include <Interprete.h>
#include <Domaine.h>
#include <Zone_VF.h>
#include <Champ_Uniforme.h>
#include <Milieu_base.h>

Implemente_instanciable( Source_Term_Thermique_PEMFC, "Source_Term_Thermique_PEMFC_VEF_P1NC", Source_base ) ;

Sortie& Source_Term_Thermique_PEMFC::printOn( Sortie& os ) const
{
  Source_base::printOn( os );
  return os;
}

Entree& Source_Term_Thermique_PEMFC::readOn( Entree& is )
{
  Source_base::readOn( is );
  Cerr << " Source_Term_Thermique_PEMFC::readOn " << finl  ;
  Param param(que_suis_je());
  set_param(param);
  param.lire_avec_accolades(is);

  dom_ = equation().probleme().domaine();
  dom_.valeur().creer_tableau_elements(Q_);
  dom_.valeur().creer_tableau_elements(dQdT_);
  return is;
  return is;
}

DoubleTab& Source_Term_Thermique_PEMFC::ajouter(DoubleTab& resu) const
{
  assert(resu.dimension(0)==la_zone_.valeur().nb_faces());

  DoubleVect vol = la_zone_.valeur().volumes();
  const DoubleTab& Cp = Cp_.valeur().valeurs();			// P0
  const DoubleTab& rho = rho_.valeur().valeurs();		// P0

  double invrhoCp = 1. / (Cp(0,0)*rho(0,0));

  int nb_elem = la_zone_.valeur().nb_elem_tot();
  for (int elem = 0; elem < nb_elem; ++elem)
    {
      double Q = Q_(elem)*invrhoCp;
      int nb_face_elem = la_zone_.valeur().zone().nb_faces_elem(0);

      for (int f = 0; f < nb_face_elem; ++f)
        {
          int face = la_zone_.valeur().elem_faces(elem, f);
          resu(face) += Q * vol(elem)/ nb_face_elem;
        }
    }

  return resu;
}

DoubleTab& Source_Term_Thermique_PEMFC::calculer(DoubleTab& resu) const
{
  resu = 0;
  ajouter(resu);
  return resu;
}


void Source_Term_Thermique_PEMFC::mettre_a_jour(double temps)
{
  Cerr << "Source_Term_Thermique_PEMFC::mettre_a_jour " << equation().probleme().le_nom() << finl;
  ch_Q_.valeur().mettre_a_jour(temps);
  ch_dQdT_.valeur().mettre_a_jour(temps);

  const DoubleTab& xp=la_zone_.valeur().xp(); // Recuperation des centre de gravite des elements pour P0
  ch_Q_.valeur().valeur_aux( xp, Q_ );
  ch_dQdT_.valeur().valeur_aux( xp, dQdT_ );
  Q_.echange_espace_virtuel();
  dQdT_.echange_espace_virtuel();
  Cerr << "Qreac min max " << mp_min_vect(Q_) << " " << mp_max_vect(Q_) << finl;
  Cerr << "DQreacDT min max " << mp_min_vect(dQdT_) << " " << mp_max_vect(dQdT_) << finl;
}

void Source_Term_Thermique_PEMFC::associer_zones(const Zone_dis& zone_dis, const Zone_Cl_dis& zcl_dis)
{
  la_zone_ = ref_cast(Zone_VF,zone_dis.valeur());
}

void Source_Term_Thermique_PEMFC::associer_pb(const Probleme_base& pb)
{
  Cerr << " Source_Term_Thermique_PEMFC::associer_pb " << finl ;
  assert(pb.que_suis_je() == "Pb_Conduction");
  const Champ_Don& Cp = pb.milieu().capacite_calorifique();
  const Champ_Don& rho = pb.milieu().masse_volumique();
  rho_ = rho;
  Cp_ = Cp;
  assert(sub_type(Champ_Uniforme,Cp_.valeur().valeur()));
  assert(sub_type(Champ_Uniforme,rho_.valeur().valeur()));
}

void Source_Term_Thermique_PEMFC::completer()
{
  Cerr << "Source_Term_Thermique_PEMFC::completer " << equation().probleme().le_nom() << finl;
  // get the reference to the coupling fields
  Source_base::completer();
  Probleme_base& pb = ref_cast(Probleme_base,interprete().objet(nom_pb_));
  ch_Q_ = pb.get_champ(nom_champ_);
  ch_dQdT_ = pb.get_champ(nom_derivee_);
  assert(ch_Q_.valeur().que_suis_je().find("P0") !=-1);
  assert(ch_dQdT_.valeur().que_suis_je().find("P0") !=-1);
}

void Source_Term_Thermique_PEMFC::contribuer_a_avec(const DoubleTab& inco, Matrice_Morse& mat) const
{
  DoubleVect vol = la_zone_.valeur().volumes();
  const DoubleTab& Cp = Cp_.valeur().valeurs();			// P0
  const DoubleTab& rho = rho_.valeur().valeurs();		// P0

  double invrhoCp = 1. / (Cp(0,0)*rho(0,0));

  int nb_elem = la_zone_.valeur().nb_elem_tot();
  for (int elem = 0; elem < nb_elem; ++elem)
    {
      double dQdT = dQdT_(elem)*invrhoCp;
      int nb_face_elem = la_zone_.valeur().zone().nb_faces_elem(0);

      for (int f = 0; f < nb_face_elem; ++f)
        {
          int face = la_zone_.valeur().elem_faces(elem, f);
          mat.coef(face,face) += - dQdT * vol(elem)/ nb_face_elem;			// -dQdT
        }
    }
}

void Source_Term_Thermique_PEMFC::set_param(Param& param)
{
  param.ajouter("nom_pb", &nom_pb_, Param::REQUIRED);
  param.ajouter("nom_champ", &nom_champ_, Param::REQUIRED);
  param.ajouter("nom_derivee", &nom_derivee_, Param::REQUIRED);
}
