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
#include <Champ_P1NC.h>

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

// integrale volumique de -div(nd/F*I)/((1-por)*eps)
DoubleTab& Source_Terme_Nafion_electro_osmotic::ajouter(DoubleTab& resu) const
{


  // ajouter un terme source de type: -nd/F*I_i
  const DoubleTab& face_norm = la_zone_VEF.valeur().face_normales();
  //const DoubleTab& facette_norm = la_zone_VEF.valeur().facette_normales();
  const IntTab& face_voisins = la_zone_VEF.valeur().face_voisins();

  double surface, signe;

  const DoubleTab& por = por_naf_.valeurs();
  const DoubleTab& eps = eps_naf_.valeurs();

  //DoubleTab Ce;
  //equation().zone_dis().zone().creer_tableau_elements(Ce);
  //const DoubleTab& xp = la_zone_VEF.valeur().xp(); // centre de gravite des elements pour P0
  //ch_C_.valeur().valeur().valeur_aux(xp, Ce);
  //Ce.echange_espace_virtuel();

  const int nint=la_zone_VEF.valeur().premiere_face_int();
  const int nb_faces=la_zone_VEF.valeur().nb_faces();
  assert(nb_faces==resu.dimension(0));
  for(int f=nint; f<nb_faces; f++)
    {

      surface = la_zone_VEF.valeur().surface(f);
      for(int v=0; v<face_voisins.dimension(1); v++)
        {
          int e = face_voisins(f,v);
          assert(e!= -1);
          double coef = 1./(F*(1-por(e,0))*eps(e,0));
          signe = la_zone_VEF.valeur().oriente_normale(f,e);
          //double Ce1 = Ce(e);
          double C = C_(e);
          double nd = f_nd(C);
          double res = 0;
          for (int j = 0; j < dimension; ++j)
            {
              res += I_(e, j)*face_norm(f,j);							// produit scalaire I*norma
            }
          resu(f) +=  -signe*nd*res*surface*coef;
        }
    }

  /*
  DoubleTab nd_face, nd_elem, grad(0, 1, Objet_U::dimension);
  la_zone_VEF.valeur().creer_tableau_faces(nd_face);
  equation().zone_dis().zone().creer_tableau_elements(nd_elem);
  equation().zone_dis().zone().creer_tableau_elements(grad);
  int nb_faces = la_zone_VEF.valeur().nb_faces();
  for (int face = 0; face < nb_faces; ++face)
    {
      nd_face(face) = f_nd(C_(face));
    }
  Champ_P1NC::calcul_gradient(nd_face,grad,ref_cast(Zone_Cl_VEF,equation().zone_Cl_dis().valeur()));
  grad.echange_espace_virtuel();
  const DoubleTab& xp = la_zone_VEF.valeur().xp(); // centre de gravite des elements pour P0
  ch_C_.valeur().valeur().valeur_aux(xp, nd_elem);
  int nb_elem = la_zone_VEF.valeur().nb_elem();
  for (int elem = 0; elem < nb_elem; ++elem)
    {
      nd_elem(elem) = f_nd(nd_elem(elem));
    }

  // Q1 = - integral_volum{(grad(nd)*Ii)/(F*(1-por)*eps)}
  DoubleVect vol = la_zone_VEF.valeur().volumes();
  const IntTab& elem_faces = la_zone_VEF.valeur().elem_faces();
  int nb_face_elem = la_zone_VEF.valeur().zone().nb_faces_elem(0);
  for (int elem = 0; elem < nb_elem; elem++)
    {
      double coef = F*(1-por(elem,0))*eps(elem,0);
      double res = 0;
      for (int j = 0; j < dimension; ++j)
        {
          res += I_(elem, j)*grad(elem,0,j);							// produit scalaire I*grad(nd)
        }
      for (int f = 0; f < nb_face_elem; ++f)
        {
          int face = elem_faces(elem, f);
          resu(face) += - res / coef * vol(elem)/nb_face_elem;
        }
    }

  // Q2 = - integral_surf{(Ii*norma_face)*nd/(F*(1-por)*eps)}
  const IntTab& face_voisins = la_zone_VEF.valeur().face_voisins();

  for (int face = 0; face < nb_faces; ++face)
    {
      double face_surf = la_zone_VEF.valeur().face_surfaces(face);
      int n0 = face_voisins(face,0);
      int n1 = face_voisins(face,1);
      if (n0 != -1)
        {
          double coef = F*(1-por(n0,0))*eps(n0,0);
          double res = 0;
          for (int j = 0; j < dimension; ++j)
            {
              res += I_(n0, j)*face_norm(face,j);							// produit scalaire I*norma
            }
          resu(face) += -res * nd_elem(n0) / coef * face_surf;				// DVQ: corriger la signe " - "
        }
      if (n1 != -1)
        {
          double coef = F*(1-por(n1,0))*eps(n1,0);
          double res = 0;
          for (int j = 0; j < dimension; ++j)
            {
              res += I_(n1, j)*face_norm(face,j);							// produit scalaire I*norma
            }
          resu(face) += -res * nd_elem(n1) / coef * face_surf ;
        }
    }
    */
  /*
  for (int elem = 0; elem < la_zone_VEF.valeur().nb_elem(); elem++)
  {
    double coef = (1-por(elem,0))*eps(elem,0);
    double nd_sur_F = f_nd(C_(elem))/F;

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
        resu(face) +=  nd_sur_F * res * face_surf / coef ;	     // A VERIFIER
      }
  }
  */


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
  Cerr << "Source_Terme_Nafion_electro_osmotic::mettre_a_jour" << equation().probleme().le_nom() << finl;
  ch_C_.valeur().mettre_a_jour(temps);
  //C_.ref(ch_C_.valeur().valeurs());

  const DoubleTab& xp = la_zone_VEF.valeur().xp(); // centre de gravite des elements pour P0
  ch_C_.valeur().valeur().valeur_aux(xp, C_);
  C_.echange_espace_virtuel();

  if(ch_I_.non_nul())
    {
      ch_I_.valeur().mettre_a_jour(temps);
      ch_I_.valeur().valeur_aux( xp, I_ );
    }
  if(ionicCurrent_.non_nul())
    {
      ionicCurrent_.mettre_a_jour(temps);
      ionicCurrent_.valeur().valeur_aux( xp, I_ );		// note: Not really need ? Reference to values is enough if ionicCurrent is P0
    }
  I_.echange_espace_virtuel();

  Cerr << "champ de courant ionique ch_I min max " << mp_min_vect(I_) << " " << mp_max_vect(I_) << finl;
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
  equation().probleme().domaine().creer_tableau_elements(C_);

  if(nom_pb_phi_ != "??")
    {
      Probleme_base& pb_phi = ref_cast(Probleme_base,interprete().objet(nom_pb_phi_));
      ch_I_ = pb_phi.get_champ(nom_champ_I_);
    }

  I_.resize(0, dimension);
  equation().probleme().domaine().creer_tableau_elements(I_);
}

void Source_Terme_Nafion_electro_osmotic::set_param(Param& param)
{
  param.ajouter("nom_pb_phi", &nom_pb_phi_, Param::OPTIONAL);
  param.ajouter("nom_champ_i", &nom_champ_I_, Param::OPTIONAL);
  param.ajouter("ionicCurrent", &ionicCurrent_, Param::OPTIONAL);
  param.ajouter("por_naf", &por_naf_, Param::REQUIRED);
  param.ajouter("eps_naf", &eps_naf_, Param::REQUIRED);
}

void Source_Terme_Nafion_electro_osmotic::remplir_volumes()
{
  volumes_.ref(la_zone_VEF.valeur().volumes_entrelaces());
}
