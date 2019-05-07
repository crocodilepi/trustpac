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
// File      : Source_Terme_Nafion_electro_osmotic.cpp
// Directory : $PEMFC_ROOT/src/Sources
//
/////////////////////////////////////////////////////////////////////////////

#include <Source_Terme_Nafion_electro_osmotic.h>
#include <Probleme_base.h>
#include <Equation_base.h>
#include <Conduction.h>
#include <Zone_VEF.h>
#include <Zone_Cl_VEF.h>
#include <DoubleVect.h>
#include <Interprete.h>
#include <Domaine.h>

Implemente_instanciable( Source_Terme_Nafion_electro_osmotic, "Source_Terme_Nafion_electro_osmotic_VEF_P1NC", Source_base ) ;

Sortie& Source_Terme_Nafion_electro_osmotic::printOn( Sortie& os ) const
{
  Source_base::printOn( os );
  return os;
}

Entree& Source_Terme_Nafion_electro_osmotic::readOn( Entree& is )
{
  Source_base::readOn( is );
  Cerr << " Source_Terme_Nafion_electro_osmotic::readOn " << finl  ;
  Param param(que_suis_je());
  set_param(param);
  param.lire_avec_accolades(is);
  return is;
}

// integrale surfacique de (nd/F*I*n*dS)/((1-por)*eps)
DoubleTab& Source_Terme_Nafion_electro_osmotic::ajouter(DoubleTab& resu) const
{
  // ajouter un terme source de type: -nd/F*I_i
  const DoubleTab& face_norm = la_zone_VEF.valeur().face_normales();
  const DoubleTab& por = por_naf_.valeurs();
  const DoubleTab& eps = eps_naf_.valeurs();

  for (int elem = 0; elem < la_zone_VEF.valeur().nb_elem(); elem++)
    {
      double coef = (1-por(elem,0))*eps(elem,0);
      int nb_face_elem = la_zone_VEF.valeur().zone().nb_faces_elem(0);		// 3 pour triangle, 4 pour tetrahedre
      for (int f = 0; f < nb_face_elem; ++f)
        {
          int face = la_zone_VEF.valeur().elem_faces(elem, f);
          double face_surf = la_zone_VEF.valeur().face_surfaces(face);
          double res = 0;
          for (int j = 0; j < dimension; ++j)
            {
              res += I_(elem, j)*face_norm(face,j);							// produit scalaire I*norma
            }
          resu(face) +=  f_nd(C_(face))/F * res * face_surf / coef ;	     // A VERIFIER
        }
    }
  return resu;
}

DoubleTab& Source_Terme_Nafion_electro_osmotic::calculer(DoubleTab& resu) const
{
  resu = 0;
  ajouter(resu);
  return resu;
}

void Source_Terme_Nafion_electro_osmotic::mettre_a_jour(double temps)
{
  //ch_C_.valeur().mettre_a_jour(temps);
  C_.ref(ch_C_.valeur().valeurs());
  //ch_I_.valeur().mettre_a_jour(temps);
  const DoubleTab& xp = la_zone_VEF.valeur().xp(); // centre de gravite des elements pour P0
  ch_I_.valeur().valeur_aux( xp, I_ );

  Cerr << "Source_Terme_Nafion_electro_osmotic::mettre_a_jour" << finl;
  Cerr << "champ de flux d'eau ch_I min max " << mp_min_vect(I_) << " " << mp_max_vect(I_) << finl;
}

void Source_Terme_Nafion_electro_osmotic::associer_zones(const Zone_dis& zone_dis, const Zone_Cl_dis& zcl_dis)
{
  la_zone_VEF = ref_cast(Zone_VEF,zone_dis.valeur());
  la_zcl_VEF = ref_cast(Zone_Cl_VEF,zcl_dis.valeur());
  remplir_volumes();
}

void Source_Terme_Nafion_electro_osmotic::associer_pb(const Probleme_base& pb)
{
  Cerr << " Source_Term_Nafion_Reaction::associer_pb " << finl ;
  assert(pb.que_suis_je() == "Pb_Conduction");
}

void Source_Terme_Nafion_electro_osmotic::completer()
{
  Source_base::completer();
  // get the reference to the coupling fields
  ch_C_ = equation().inconnue();
  Probleme_base& pb_phi = ref_cast(Probleme_base,interprete().objet(nom_pb_phi_));
  ch_I_ = pb_phi.get_champ(nom_champ_I_);
  I_.resize(0, dimension);
  equation().probleme().domaine().creer_tableau_elements(I_);
}

void Source_Terme_Nafion_electro_osmotic::set_param(Param& param)
{
  param.ajouter("nom_pb_phi", &nom_pb_phi_, Param::REQUIRED);
  param.ajouter("nom_champ_i", &nom_champ_I_, Param::REQUIRED);
  param.ajouter("por_naf", &por_naf_, Param::REQUIRED);
  param.ajouter("eps_naf", &eps_naf_, Param::REQUIRED);
}

void Source_Terme_Nafion_electro_osmotic::remplir_volumes()
{
  volumes_.ref(la_zone_VEF.valeur().volumes_entrelaces());
}
