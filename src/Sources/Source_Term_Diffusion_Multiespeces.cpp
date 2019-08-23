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
// File      : Source_Term_Diffusion_Multiespeces.cpp
// Directory : $PEMFC_ROOT/src/Sources
//
/////////////////////////////////////////////////////////////////////////////

#include <Source_Term_Diffusion_Multiespeces.h>
#include <Probleme_base.h>
#include <Equation_base.h>
#include <Interprete.h>
#include <Domaine.h>
#include <Zone_VF.h>

Implemente_instanciable( Source_Term_Diffusion_Multiespeces, "Source_Term_Diffusion_Multiespeces_VEF_P1NC", Source_base ) ;

Sortie& Source_Term_Diffusion_Multiespeces::printOn( Sortie& os ) const
{
  Source_base::printOn( os );
  return os;
}

Entree& Source_Term_Diffusion_Multiespeces::readOn( Entree& is )
{
  Source_base::readOn( is );
  Cerr << " Source_Term_Diffusion_Multiespeces::readOn " << finl  ;
  Param param(que_suis_je());
  set_param(param);
  param.lire_avec_accolades(is);
  return is;
}

DoubleTab& Source_Term_Diffusion_Multiespeces::ajouter(DoubleTab& resu) const
{
  assert(resu.dimension(0)==la_zone_.valeur().nb_faces());
  assert(resu.dimension(1)==S_.dimension(1));			// 3 composants

  DoubleVect vol = la_zone_.valeur().volumes();
  int nb_face_elem = la_zone_.valeur().zone().nb_faces_elem(0);
  for (int elem = 0; elem < la_zone_.valeur().nb_elem(); ++elem)
    {
      for (int f = 0; f < nb_face_elem; ++f)
        {
          int face = la_zone_.valeur().elem_faces(elem, f);

          for (int comp = 0; comp < resu.dimension(1); ++comp)
            {
              resu(face, comp) += S_(elem, comp) * vol(elem) / nb_face_elem;
            }
        }
    }

  return resu;
}

DoubleTab& Source_Term_Diffusion_Multiespeces::calculer(DoubleTab& resu) const
{
  resu = 0;
  ajouter(resu);
  return resu;
}

void Source_Term_Diffusion_Multiespeces::mettre_a_jour(double temps)
{
  Cerr << "Source_Term_Diffusion_Multiespeces::mettre_a_jour " << equation().probleme().le_nom() << finl;
  // mettre a jour les champs couples
  ch_S_X2_.valeur().mettre_a_jour(temps);
  ch_S_H2O_.valeur().mettre_a_jour(temps);
  ch_S_N2_.valeur().mettre_a_jour(temps);

  const DoubleTab& xp=la_zone_.valeur().xp(); // centre de gravite des faces pour P0

  // interpolation P0
  DoubleTab S_X2, S_vap, S_N2;
  la_zone_.valeur().zone().creer_tableau_elements(S_X2);
  la_zone_.valeur().zone().creer_tableau_elements(S_vap);
  la_zone_.valeur().zone().creer_tableau_elements(S_N2);
  ch_S_X2_.valeur().valeur_aux( xp, S_X2 );
  ch_S_H2O_.valeur().valeur_aux( xp, S_vap );
  ch_S_N2_.valeur().valeur_aux( xp, S_N2 );

  int nb_elem = la_zone_.valeur().nb_elem();
  for (int elem = 0; elem < nb_elem; ++elem)
    {
      S_(elem, 0) = -S_X2(elem);
      S_(elem, 1) = -S_vap(elem);
      S_(elem, 2) = -S_N2(elem);
    }

  Cerr << "champ de source de diffusion Sa_X2 min max " << mp_min_vect(S_X2) << " " << mp_max_vect(S_X2) << finl;
  Cerr << "champ de source de diffusion Sa_vap min max " << mp_min_vect(S_vap) << " " << mp_max_vect(S_vap) << finl;
  Cerr << "champ de source de diffusion Sa_N2 min max " << mp_min_vect(S_N2) << " " << mp_max_vect(S_N2) << finl;
}

void Source_Term_Diffusion_Multiespeces::associer_zones(const Zone_dis& zone_dis, const Zone_Cl_dis& zcl_dis)
{
  la_zone_ = ref_cast(Zone_VF,zone_dis.valeur());
}

void Source_Term_Diffusion_Multiespeces::associer_pb(const Probleme_base& pb)
{
  // recuperer le milieu_base -> nothing to do
  Cerr << " Source_Term_Diffusion_Multiespeces::associer_pb " << finl ;
  assert(pb.que_suis_je() == "Pb_Hydraulique_Concentration");
}

void Source_Term_Diffusion_Multiespeces::completer()
{
  Cerr << " Source_Term_Diffusion_Multiespeces::completer " << finl ;

  // get all references to the coupling fields
  Source_base::completer();

  dom_ = equation().probleme().domaine();
  int dim = equation().inconnue().valeur().nb_comp();			// =3, multi-especes = X2, vap, N2
  S_.resize(0, dim);
  la_zone_.valeur().zone().creer_tableau_elements(S_);

  Probleme_base& pb_X2 = ref_cast(Probleme_base,interprete().objet(nom_pb_X2_));
  ch_S_X2_ = pb_X2.get_champ(nom_champ_S_X2_);
  assert(ch_S_X2_.valeur().que_suis_je().find("P0") !=-1);

  Probleme_base& pb_H2O = ref_cast(Probleme_base,interprete().objet(nom_pb_H2O_));
  ch_S_H2O_ = pb_H2O.get_champ(nom_champ_S_H2O_);
  assert(ch_S_H2O_.valeur().que_suis_je().find("P0") !=-1);

  Probleme_base& pb_N2 = ref_cast(Probleme_base,interprete().objet(nom_pb_N2_));
  ch_S_N2_ = pb_N2.get_champ(nom_champ_S_N2_);
  assert(ch_S_N2_.valeur().que_suis_je().find("P0") !=-1);
}

void Source_Term_Diffusion_Multiespeces::set_param(Param& param)
{
  param.ajouter("nom_pb_X2", &nom_pb_X2_, Param::REQUIRED);
  param.ajouter("nom_champ_S_X2", &nom_champ_S_X2_, Param::REQUIRED);
  param.ajouter("nom_pb_H2O", &nom_pb_H2O_, Param::REQUIRED);
  param.ajouter("nom_champ_S_H2O", &nom_champ_S_H2O_, Param::REQUIRED);
  param.ajouter("nom_pb_N2", &nom_pb_N2_, Param::REQUIRED);
  param.ajouter("nom_champ_S_N2", &nom_champ_S_N2_, Param::REQUIRED);
}
