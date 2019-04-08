/****************************************************************************
* Copyright (c) 2017, CEA
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
// File      : Op_Diff_VEF_Tensor_base.cpp
// Directory : $DIFFUSION_ANISOTROPE_ROOT/src
//
/////////////////////////////////////////////////////////////////////////////

#include <Op_Diff_VEF_Tensor_base.h>
#include <Champ_P1NC.h>
#include <Champ_Q1NC.h>
#include <Ref_Champ_P1NC.h>
#include <Ref_Champ_Q1NC.h>
#include <Champ_Uniforme.h>
#include <Milieu_base.h>
#include <Probleme_base.h>
#include <Schema_Temps_base.h>
#include <Champ_Fonc_P0_base.h>
#include <Discretisation_base.h>
#include <Champ.h>
#include <Check_espace_virtuel.h>
#include <Debog.h>

#include <Domaine.h>
#include <MD_Vector.h>
#include <MD_Vector_composite.h>

Implemente_base( Op_Diff_VEF_Tensor_base, "Op_Diff_VEF_Tensor_base", Operateur_Diff_base ) ;

//// printOn
//

Sortie& Op_Diff_VEF_Tensor_base::printOn(Sortie& s ) const
{
  return s << que_suis_je() ;
}


//// readOn
//

Entree& Op_Diff_VEF_Tensor_base::readOn(Entree& s )
{
  return s ;
}


// Description : definit si on calcule div(phi nu grad Psi) ou div(nu grap Phi psi)
int Op_Diff_VEF_Tensor_base::phi_psi_diffuse(const Equation_base& eq) const
{
  if (eq.inconnue().le_nom()=="vitesse")
    return 1;
  return 0;
}

int Op_Diff_VEF_Tensor_base::impr(Sortie& os) const
{
  return Op_VEF_Face::impr(os, *this);
}

void Op_Diff_VEF_Tensor_base::associer(const Zone_dis& zone_dis,
                                       const Zone_Cl_dis& zone_cl_dis,
                                       const Champ_Inc& ch_transporte)
{

  const Zone_VEF& zvef = ref_cast(Zone_VEF,zone_dis.valeur());
  const Zone_Cl_VEF& zclvef = ref_cast(Zone_Cl_VEF,zone_cl_dis.valeur());

  if (sub_type(Champ_P1NC,ch_transporte.valeur()))
    {
      const Champ_P1NC& inco = ref_cast(Champ_P1NC,ch_transporte.valeur());
      REF(Champ_P1NC) inconnue;
      inconnue = inco;
    }
  if (sub_type(Champ_Q1NC,ch_transporte.valeur()))
    {
      const Champ_Q1NC& inco = ref_cast(Champ_Q1NC,ch_transporte.valeur());
      REF(Champ_Q1NC) inconnue;
      inconnue = inco;
    }

  la_zone_vef = zvef;
  la_zcl_vef = zclvef;


  const Zone_VEF& zone_VEF = la_zone_vef.valeur();
  int nb_dim = ch_transporte->valeurs().nb_dim();
  int nb_comp = 1;

  if(nb_dim==2)
    nb_comp = ch_transporte->valeurs().dimension(1);

  (ref_cast(DoubleTab,flux_bords_)).resize(zone_VEF.nb_faces_bord(),nb_comp);
  flux_bords_=0.;
}


double Op_Diff_VEF_Tensor_base::calculer_dt_stab() const
{
  //Cerr << "remplir_nu in the method: calculer_dt_stab" << finl;
  remplir_nu(nu_);
  // La diffusivite est constante dans le domaine donc
  //
  //          dt_diff = h*h/diffusivite

  double dt_stab = DMAXFLOAT;
  const Zone_VEF& zone_VEF = la_zone_vef.valeur();
  if (! has_champ_masse_volumique())
    {

      // Methode "standard" de calcul du pas de temps
      // Ce calcul est tres conservatif: si le max de la diffusivite
      // n'est pas atteint a l'endroit ou le min de delta_h_carre est atteint,
      // le pas de temps est sous-estime.
      const Champ_base& champ_diffusivite = diffusivite_pour_pas_de_temps();
      const DoubleVect&      valeurs_diffusivite = champ_diffusivite.valeurs();
      double alpha_max = local_max_vect(valeurs_diffusivite);

      if ((alpha_max == 0.)||(zone_VEF.nb_elem()==0))
        dt_stab = DMAXFLOAT;
      else
        {
          const double min_delta_h_carre = zone_VEF.carre_pas_du_maillage();
          dt_stab = min_delta_h_carre / (2. * dimension * alpha_max);
        }

    }
  else
    {
      const double           deux_dim      = 2. * Objet_U::dimension;
      const Champ_base& champ_diffu   = diffusivite();
      const DoubleTab&       valeurs_diffu = champ_diffu.valeurs();
      const Champ_base&      champ_rho     = get_champ_masse_volumique();
      const DoubleTab&       valeurs_rho   = champ_rho.valeurs();
      assert(sub_type(Champ_Fonc_P0_base, champ_rho));
      assert(sub_type(Champ_Fonc_P0_base, champ_diffu));
      const int nb_elem = zone_VEF.nb_elem();
      assert(valeurs_rho.size_array()==zone_VEF.nb_elem_tot());
      // Champ de masse volumique variable.
      for (int elem = 0; elem < nb_elem; elem++)
        {
          const double h_carre = zone_VEF.carre_pas_maille(elem);
          const double diffu   = valeurs_diffu(elem);
          const double rho     = valeurs_rho(elem);
          const double dt      = h_carre * rho / (deux_dim * (diffu+DMINFLOAT));
          if (dt_stab > dt)
            dt_stab = dt;
        }
    }
  dt_stab = Process::mp_min(dt_stab);
  return dt_stab;
}

// cf Op_Diff_VEF_Tensor_base::calculer_dt_stab() pour choix de calcul de dt_stab
void Op_Diff_VEF_Tensor_base::calculer_pour_post(Champ& espace_stockage,const Nom& option,int comp) const
{
  if (Motcle(option)=="stabilite")
    {
      DoubleTab& es_valeurs = espace_stockage->valeurs();

      if (la_zone_vef.non_nul())
        {
          remplir_nu(nu_);
          const Zone_VEF& zone_VEF = la_zone_vef.valeur();
          if (! has_champ_masse_volumique())
            {

              const Champ_base& champ_diffusivite = diffusivite_pour_pas_de_temps();
              const DoubleVect&      valeurs_diffusivite = champ_diffusivite.valeurs();
              double alpha_max = local_max_vect(valeurs_diffusivite);

              if ((alpha_max == 0.)||(zone_VEF.nb_elem()==0))
                es_valeurs = DMAXFLOAT;
              else
                {
                  const int nb_elem = zone_VEF.nb_elem();
                  for (int elem = 0; elem < nb_elem; elem++)
                    {
                      es_valeurs(elem) = zone_VEF.carre_pas_maille(elem) / (2. * dimension * alpha_max);
                    }
                }
            }
          else
            {
              const double           deux_dim      = 2. * Objet_U::dimension;
              const Champ_base& champ_diffu   = diffusivite();
              const DoubleTab&       valeurs_diffu = champ_diffu.valeurs();
              const Champ_base&      champ_rho     = get_champ_masse_volumique();
              const DoubleTab&       valeurs_rho   = champ_rho.valeurs();
              assert(sub_type(Champ_Fonc_P0_base, champ_rho));
              assert(sub_type(Champ_Fonc_P0_base, champ_diffu));
              const int nb_elem = zone_VEF.nb_elem();
              assert(valeurs_rho.size_array()==zone_VEF.nb_elem_tot());
              // Champ de masse volumique variable.
              for (int elem = 0; elem < nb_elem; elem++)
                {
                  const double h_carre = zone_VEF.carre_pas_maille(elem);
                  const double diffu   = valeurs_diffu(elem);
                  const double rho     = valeurs_rho(elem);
                  const double dt      = h_carre * rho / (deux_dim * diffu);
                  es_valeurs(elem) = dt;
                }
            }

          assert(mp_min_vect(es_valeurs)==calculer_dt_stab());
        }
    }
  else
    Operateur_Diff_base::calculer_pour_post(espace_stockage,option,comp);
}

Motcle Op_Diff_VEF_Tensor_base::get_localisation_pour_post(const Nom& option) const
{
  Motcle loc;
  if (Motcle(option)=="stabilite")
    loc = "elem";
  else
    return Operateur_Diff_base::get_localisation_pour_post(option);
  return loc;
}

// remplir les valeurs de diffusivité dans un tableau (soit 1-dim soit 3-dim)
void Op_Diff_VEF_Tensor_base::remplir_nu(DoubleTab& nu) const
{
  // original code
  const Zone_VEF& zone_VEF = la_zone_vef.valeur();
  // On dimensionne nu
  const DoubleTab& diffu=diffusivite().valeurs();

  if (!nu.get_md_vector().non_nul())
    {
      if (diffu.nb_dim() > 1 && diffu.dimension(1) > 1)
        {
          nu.resize(0, diffu.dimension(1));
        }
      zone_VEF.zone().creer_tableau_elements(nu, Array_base::NOCOPY_NOINIT);
    }
  if (! diffu.get_md_vector().non_nul())
    {
      // diffusvite uniforme
      const int n = nu.dimension_tot(0);
      const int nb_comp = nu.line_size();
      // Tableaux vus comme uni-dimenionnels:
      const DoubleVect& arr_diffu = diffu;
      DoubleVect& arr_nu = nu;
      for (int i = 0; i < n; i++)
        for (int j = 0; j < nb_comp; j++)
          arr_nu[i*nb_comp + j] = arr_diffu[j];
    }
  else
    {
      assert( nu.get_md_vector( ) == diffu.get_md_vector( ) );
      assert_espace_virtuel_vect(diffu);
      nu.inject_array(diffu);
    }
  // end original code

  /************************ dinhvan: changement pour tranformer le tableau 2D to un tableau 3D  ************************/
  /*
  const Zone_VEF& zone_VEF = la_zone_vef.valeur();
  // On dimensionne nu
  const DoubleTab& diffu=diffusivite().valeurs();
  // cas scalaire (tenseur 1x1): diffu est un tableau 1-dim
  // dinhvan 28/08/2018 jeu de donnees "lambda champ_fonc_xyz dom 1 1.0" --> diffusivite_ est un tableau de 2-dim [nb_elemx1], ce qui ne fonctionne pas avec la condition "if(diffu.nb_dim()==1)"
  // remplace par assert(diffu.nb_dim()==2) et if(diffu.dimension(0)==1)
  //if(diffu.nb_dim()==1)

  assert(diffu.nb_dim()==2);

  // cas scalaire
  if(diffu.dimension(1)==1)
  {
    if (!nu.get_md_vector().non_nul())
      {
        if (diffu.nb_dim() > 1 && diffu.dimension(1) > 1)
          nu.resize(0, diffu.dimension(1));
        zone_VEF.zone().creer_tableau_elements(nu, Array_base::NOCOPY_NOINIT);
      }
    if (! diffu.get_md_vector().non_nul())
      {
        // diffusvite uniforme
        const int n = nu.dimension_tot(0);
        const int nb_comp = nu.line_size();
        // Tableaux vus comme uni-dimenionnels:
        const DoubleVect& arr_diffu = diffu;
        DoubleVect& arr_nu = nu;
        for (int i = 0; i < n; i++)
          for (int j = 0; j < nb_comp; j++)
            arr_nu[i*nb_comp + j] = arr_diffu[j];
      }
    else
      {
        assert(nu.get_md_vector() == diffu.get_md_vector());
        assert_espace_virtuel_vect(diffu);
        nu.inject_array(diffu);
      }
  }
  // tenseur 2x2 ou 3x3 selon la dimension: diffu est un tableau de 2-dim
  // Attention: ne pas encore tenir compte le calcul parallele
  else
  {
    nu.resize(diffu.dimension(0),dimension,dimension);
    for (int i=0; i< nu.dimension(0); i++)
      {
        for(int j=0; j<dimension; j++)
          {
            for(int k=0; k<dimension; k++)
              {
                nu(i,j,k)=diffu(i,dimension*j+k);
              }
          }
      }
  }

  // check diffu et nu
  //  ArrOfDouble tmp(diffu);
  //  tmp.resize_array(9);
  //  Cerr << "tmp_diffu dim=" << tmp.size_array() << finl;
  //  Cerr << tmp << finl;
  //  Cerr << "remplir nu, nb_dim=" << nu.nb_dim() << ", dim= " << nu.dimension(0) << "x" << nu.dimension(1) << "x" << nu.dimension(2) << finl;
  */
  /************************ dinhvan: changement pour tranformer le tableau 2D to un tableau 3D  ************************/
}
