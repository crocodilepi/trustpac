/****************************************************************************
* Copyright (c) 2015 - 2016, CEA
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
//////////////////////////////////////////////////////////////////////////////
//
// File:        Perte_Charge_Cas5_VEF_P1NC.cpp
// Directory:   $TRUST_ROOT/src/VEF/Sources
// Version:     /main/20
//
//////////////////////////////////////////////////////////////////////////////

#include <Perte_Charge_Cas5_VEF_P1NC.h>
#include <Zone_VEF.h>
#include <Fluide_Incompressible.h>
#include <Equation_base.h>
#include <Probleme_base.h>
#include <Domaine.h>
#include <Champ_Uniforme.h>
#include <Champ_Inc.h>
#include <Matrice_Morse.h>
#include <Param.h>

Implemente_instanciable(Perte_Charge_Cas5_VEF_P1NC,"Perte_Charge_Cas5_QC_VEF_P1NC",Perte_Charge_VEF_Face);

// XD Perte_Charge_Cas5 source_base Perte_Charge_Cas5 1 terme source pour le cas 5 -mu_sur_K*U

//// printOn
//

Sortie& Perte_Charge_Cas5_VEF_P1NC::printOn(Sortie& s ) const
{
  return s << que_suis_je() ;
}


//// readOn
//

Entree& Perte_Charge_Cas5_VEF_P1NC::readOn(Entree& s )
{
  Cerr << " Perte_Charge_Cas5_VEF_P1NC::readOn " << finl  ;
  Param param(que_suis_je());
  param.ajouter("mu_sur_K",&mu_sur_K_,Param::REQUIRED) ; // XD_ADD_P field_base mu_sur_K aux elems
  param.lire_avec_accolades(s);
  return s ;
}



DoubleTab& Perte_Charge_Cas5_VEF_P1NC::ajouter(DoubleTab& resu) const
{
  //Cerr << " Perte_Charge_Cas5_VEF_P1NC::ajouter " << finl;
  const Zone_VEF& zone_VEF = ref_cast(Zone_VEF,equation().zone_dis().valeur());
  const DoubleVect& volumes_entrelaces = zone_VEF.volumes_entrelaces();
  const DoubleVect& porosite_face = zone_VEF.porosite_face();
  //const DoubleTab& vit = la_vitesse->valeurs();
  const DoubleTab& vit= equation().probleme().get_champ("Um").valeurs();
  const IntTab& face_voisins = zone_VEF.face_voisins();

  const DoubleVect& visco = mu_sur_K_.valeurs();

  int coeff_unif=0;
  double mu_sur_k=1e99;
  if (sub_type(Champ_Uniforme,mu_sur_K_.valeur()))
    {
      const Champ_Uniforme& ch = ref_cast(Champ_Uniforme,mu_sur_K_.valeur());
      mu_sur_k = ch(0,0);
      coeff_unif = 1;
    }



  int premiere_face_int=zone_VEF.premiere_face_int();
  int nb_faces =zone_VEF.nb_faces();
  for (int numfa=0; numfa<nb_faces; numfa++)
    {

      if (!coeff_unif)
        {
          int n0 = face_voisins(numfa,0);
          int n1 = face_voisins(numfa,1);
          if (numfa < premiere_face_int)
            {
              if (n0 != -1)
                mu_sur_k = visco[n0];
              else
                mu_sur_k = visco[n1];
            }
          else
            mu_sur_k = 0.5*(visco[n0]+visco[n1]);
        }
      for (int direction=0; direction<dimension; direction++)
        resu(numfa,direction) -= mu_sur_k*vit(numfa,direction)*volumes_entrelaces(numfa)*porosite_face(numfa);

    }

  return resu;
}

void  Perte_Charge_Cas5_VEF_P1NC::contribuer_a_avec(const DoubleTab&, Matrice_Morse& matrice) const
{
  const Zone_VEF& zone_VEF = ref_cast(Zone_VEF,equation().zone_dis().valeur());
  const DoubleVect& volumes_entrelaces = zone_VEF.volumes_entrelaces();
  const DoubleVect& porosite_face = zone_VEF.porosite_face();
  //const DoubleTab& vit = la_vitesse->valeurs();
  const IntTab& face_voisins = zone_VEF.face_voisins();

  const DoubleVect& visco = mu_sur_K_.valeurs();

  int coeff_unif=0;
  double mu_sur_k=1e99;
  if (sub_type(Champ_Uniforme,mu_sur_K_.valeur()))
    {
      const Champ_Uniforme& ch = ref_cast(Champ_Uniforme,mu_sur_K_.valeur());
      mu_sur_k = ch(0,0);
      coeff_unif = 1;
    }



  int premiere_face_int=zone_VEF.premiere_face_int();
  int nb_faces =zone_VEF.nb_faces();
  for (int numfa=0; numfa<nb_faces; numfa++)
    {

      if (!coeff_unif)
        {
          int n0 = face_voisins(numfa,0);
          int n1 = face_voisins(numfa,1);
          if (numfa < premiere_face_int)
            {
              if (n0 != -1)
                mu_sur_k = visco[n0];
              else
                mu_sur_k = visco[n1];
            }
          else
            mu_sur_k = 0.5*(visco[n0]+visco[n1]);
        }
      for (int direction=0; direction<dimension; direction++)
        {
          //resu(numfa,direction) -= mu_sur_k*vit(numfa,direction)*volumes_entrelaces(numfa)*porosite_face(numfa);
          int n0bis=numfa*dimension+direction;
          matrice.coef(n0bis,n0bis)+=mu_sur_k*volumes_entrelaces(numfa)*porosite_face(numfa);
        }
    }

}

DoubleTab& Perte_Charge_Cas5_VEF_P1NC::calculer(DoubleTab& resu) const
{
  resu = 0;
  return ajouter(resu);
}


void Perte_Charge_Cas5_VEF_P1NC::completer()
{
  Source_base::completer();
}

