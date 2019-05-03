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
#include <Zone_VEF.h>
#include <Zone_Cl_VEF.h>
#include <Probleme_base.h>
#include <Equation_base.h>
#include <Interprete.h>
#include <Domaine.h>

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
  assert(resu.dimension(0)==volumes_.size());
  assert(resu.dimension(0)==S_.dimension(0));
  assert(resu.dimension(1)==S_.dimension(1));

  int nb_faces = la_zone_VEF.valeur().nb_faces();
  for (int face = 0; face < nb_faces; ++face)
    {
      // champ ir doit etre mettre_a_jour() avant ajouter()
      for (int comp = 0; comp < resu.dimension(1); ++comp)
        {
          resu(face, comp) += - S_(face, comp) * volumes_(face);					// -> NEEDS VERIFYING
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
  Cerr << " Source_Term_Diffusion_Multiespeces::mettre_a_jour " << finl ;
  // mettre a jour les champs couples
  ch_S_X2_.valeur().mettre_a_jour(temps);
  ch_S_H2O_.valeur().mettre_a_jour(temps);
  ch_S_N2_.valeur().mettre_a_jour(temps);

  const DoubleTab& xv=la_zone_VEF.valeur().xv(); // centre de gravite des faces pour P1NC

  // interpolation vers P1NC
  DoubleTab S_X2, S_vap, S_N2;
  la_zone_VEF.valeur().creer_tableau_faces(S_X2);
  la_zone_VEF.valeur().creer_tableau_faces(S_vap);
  la_zone_VEF.valeur().creer_tableau_faces(S_N2);
  ch_S_X2_.valeur().valeur_aux( xv, S_X2 );
  ch_S_H2O_.valeur().valeur_aux( xv, S_vap );
  ch_S_N2_.valeur().valeur_aux( xv, S_N2 );

  int nb_faces = la_zone_VEF.valeur().nb_faces();
  for (int face = 0; face < nb_faces; ++face)
    {
      S_(face, 0) = S_X2(face);
      S_(face, 1) = S_vap(face);
      S_(face, 2) = S_N2(face);
    }

}

void Source_Term_Diffusion_Multiespeces::associer_zones(const Zone_dis& zone_dis, const Zone_Cl_dis& zcl_dis)
{
  la_zone_VEF = ref_cast(Zone_VEF,zone_dis.valeur());
  la_zcl_VEF = ref_cast(Zone_Cl_VEF,zcl_dis.valeur());
  remplir_volumes();
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
  la_zone_VEF.valeur().creer_tableau_faces(S_);

  Probleme_base& pb_X2 = ref_cast(Probleme_base,interprete().objet(nom_pb_X2_));
  ch_S_X2_ = pb_X2.get_champ(nom_champ_S_X2_);
  assert(ch_S_X2_.valeur().que_suis_je().find("P1NC") !=-1);

  Probleme_base& pb_H2O = ref_cast(Probleme_base,interprete().objet(nom_pb_H2O_));
  ch_S_H2O_ = pb_H2O.get_champ(nom_champ_S_H2O_);
  assert(ch_S_H2O_.valeur().que_suis_je().find("P1NC") !=-1);

  Probleme_base& pb_N2 = ref_cast(Probleme_base,interprete().objet(nom_pb_N2_));
  ch_S_N2_ = pb_N2.get_champ(nom_champ_S_N2_);
  assert(ch_S_N2_.valeur().que_suis_je().find("P1NC") !=-1);
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

void Source_Term_Diffusion_Multiespeces::remplir_volumes()
{
  volumes_.ref(la_zone_VEF.valeur().volumes_entrelaces());
}
