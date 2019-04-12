/****************************************************************************
* Copyright (c) 2015, CEA
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
// File:        Op_Diff_VEF_Face_Matricial.cpp
// Directory:   $TRUST_ROOT/src/VEF/Operateurs
// Version:     /main/55
//
//////////////////////////////////////////////////////////////////////////////

#include <Op_Diff_VEF_Face_Matricial.h>
#include <Champ_P1NC.h>
#include <Champ_Q1NC.h>
#include <Champ_Uniforme.h>
#include <Param.h>

#include <DoubleTrav.h>
#include <Porosites_champ.h>


#include <Periodique.h>
#include <Dirichlet_paroi_fixe.h>
#include <Dirichlet_paroi_defilante.h>
#include <Dirichlet_entree_fluide.h>
#include <Neumann_paroi.h>
#include <Echange_externe_impose.h>
#include <Neumann_homogene.h>
#include <Symetrie.h>
#include <Echange_global_impose.h>
#include <Temperature_imposee_paroi.h>
#include <Sortie_libre_pression_moyenne_imposee.h>
#include <Debog.h>
//

Implemente_instanciable(Op_Diff_VEF_Face_Matricial,"Op_Diff_VEFMATRICIAL_const_P1NC",Op_Diff_VEF_Face);
// XD diffusion_matricial diffusion_deriv matricial 1 matricial diffusion (only in VEF)

//// printOn
//

Sortie& Op_Diff_VEF_Face_Matricial::printOn(Sortie& s ) const
{
  return s << que_suis_je() ;
}

//// readOn
//

Entree& Op_Diff_VEF_Face_Matricial::readOn(Entree& s )
{
  Param param(que_suis_je());
  set_param(param);
  param.lire_avec_accolades_depuis(s);
  associer_diffusivite ( diffusivity_read_from_datafile_.valeur( ) );

  return s ;
}

void Op_Diff_VEF_Face_Matricial::set_param(Param& param)
{
  param.ajouter("Diffusivity",&diffusivity_read_from_datafile_); // XD_ADD_P field_base diffusivity coefficients

}


void Op_Diff_VEF_Face_Matricial::check_diffusivity( const int& nb_comp ) const
{
  const DoubleTab& diffusivity = diffusivite_.valeur().valeurs( );
  int nb_coeffs = diffusivity.dimension(1);
  int expected_nb_coeffs = nb_comp * nb_comp;
  assert( diffusivity.nb_dim( ) == 2 );
  //assert( diffusivity.dimension( 0 ) == 1 );
  //assert( diffusivity.dimension( 1 ) == nb_coeffs );
  if( nb_coeffs != expected_nb_coeffs )
    {
      Cerr << "Error in Op_Diff_VEF_Face_Matricial::check_diffusivity "<<finl;
      Cerr << "There are "<< nb_coeffs <<" coefficients instead of "<< expected_nb_coeffs << " as expected ( nb_comp : "<<nb_comp<<" )"<<finl;
      Cerr << "Aborting..."<<finl;
      Process::abort( );
    }
}

// Description:
// associe le champ de diffusivite
void Op_Diff_VEF_Face_Matricial::associer_diffusivite(const Champ_base& diffu)
{
  diffusivite_ = diffu;
}

const Champ_base& Op_Diff_VEF_Face_Matricial::diffusivite() const
{
  return diffusivite_.valeur();
}




DoubleTab& Op_Diff_VEF_Face_Matricial::ajouter(const DoubleTab& inconnue_org, DoubleTab& resu) const
{
  remplir_nu(nu_);
  Debog::verifier("Op_Diff_VEF_Face_Matricial inco", inconnue_org);
  const Zone_Cl_VEF& zone_Cl_VEF = la_zcl_vef.valeur();
  const Zone_VEF& zone_VEF = la_zone_vef.valeur();



  int nb_comp = 1;
  int nb_dim = resu.nb_dim();
  if(nb_dim==2)
    nb_comp=resu.dimension(1);

  check_diffusivity( nb_comp );

  // const DoubleTab& inconnue = inconnue_org;
  // const DoubleTab& nu = diffusivity_read_from_datafile_.valeurs( );
  DoubleTab nu;
  DoubleTab tab_inconnue;

  int marq=phi_psi_diffuse( equation() );
  const DoubleVect& porosite_face = zone_VEF.porosite_face();
  const DoubleVect& porosite_elem = zone_VEF.porosite_elem();
  // // soit on a div(phi nu grad inco)
  // // soit on a div(nu grad phi inco)
  // // cela depend si on diffuse phi_psi ou psi
  modif_par_porosite_si_flag(nu_,nu,!marq,porosite_elem);
  Debog::verifier("Op_Diff_VEF_Face_Matricial nu", nu);
  const DoubleTab& inconnue=modif_par_porosite_si_flag(inconnue_org,tab_inconnue,marq,porosite_face);


  const Champ_base& inco = equation().inconnue().valeur();
  const Nature_du_champ nature_champ = inco.nature_du_champ();

  if( nature_champ == multi_scalaire )
    {
      ajouter_cas_multi_scalaire( inconnue, resu, flux_bords_, nu , zone_Cl_VEF, zone_VEF, nb_comp );
    }
  else
    {
      Cerr <<"Error in Op_Diff_VEF_Face_Matricial::ajouter( )"<<finl;
      if (nature_champ==vectoriel)
        Cerr << "Vectorial fields are not expected ! "<<finl;
      if ( nature_champ == scalaire )
        Cerr << "Scalar fields are not expected ! "<<finl;
      Cerr << "Only scalars or multiscalars are expected."<<finl;
      Cerr << "Aborting..."<<finl;
      Process::abort( );
    }

  modifier_flux(*this);

  //resu.echange_espace_virtuel();
  Debog::verifier("Op_Diff_VEF_Face_Matricial resu ", resu);
  return resu;
}



void Op_Diff_VEF_Face_Matricial::ajouter_cas_multi_scalaire(const DoubleTab& inconnue,
                                                            DoubleTab& resu, DoubleTab& tab_flux_bords,
                                                            const DoubleTab& nu,
                                                            const Zone_Cl_VEF& zone_Cl_VEF,
                                                            const Zone_VEF& zone_VEF,
                                                            int nb_comp) const
{
  const IntTab& elemfaces = zone_VEF.elem_faces();
  const IntTab& face_voisins = zone_VEF.face_voisins();
  int i0,j,num_face;
  int nb_faces = zone_VEF.nb_faces();
  int elem0,elem2;
  int elem1;
  int nb_faces_elem = zone_VEF.zone().nb_faces_elem();
  double valA,flux0;
  DoubleVect n(Objet_U::dimension);
  DoubleTrav Tgrad(Objet_U::dimension,Objet_U::dimension);

  // dimensionning and initializing flow balance tabulars
  (ref_cast(DoubleTab,tab_flux_bords)).resize(zone_VEF.nb_faces_bord(),nb_comp);
  tab_flux_bords=0.;

  // assert( nb_comp>1 );
  int nb_bords=zone_VEF.nb_front_Cl();
  int ind_face;

  for (int n_bord=0; n_bord<nb_bords; n_bord++)
    {
      const Cond_lim& la_cl = zone_Cl_VEF.les_conditions_limites(n_bord);
      const Front_VF& le_bord = ref_cast(Front_VF,la_cl.frontiere_dis());

      int num1 = 0;
      int num2 = le_bord.nb_faces_tot();
      int nb_faces_bord_reel = le_bord.nb_faces();

      //periodic case
      if (sub_type(Periodique,la_cl.valeur()))
        {
          const Periodique& la_cl_perio = ref_cast(Periodique,la_cl.valeur());
          int fac_asso;


          for (ind_face=num1; ind_face<nb_faces_bord_reel; ind_face++)
            {
              fac_asso = la_cl_perio.face_associee(ind_face);
              fac_asso = le_bord.num_face(fac_asso);
              num_face = le_bord.num_face(ind_face);
              elem1 = face_voisins(num_face,0);
              for (i0=0; i0<nb_faces_elem; i0++)
                {
                  if ( ( (j= elemfaces(elem1,i0)) > num_face ) && (j != fac_asso ) )
                    {
                      for(int c1=0; c1<nb_comp; c1++)

                        for(int c2=0; c2<nb_comp; c2++)
                          {

                            int diffusivity_index = c1*nb_comp + c2;
                            valA = viscA(num_face,j,elem1,nu(elem1, diffusivity_index));
                            resu(num_face,c1)+=valA*inconnue(j,c2);
                            resu(num_face,c1)-=valA*inconnue(num_face,c2);
                            if(j<nb_faces) // real face
                              {
                                resu(j,c1)+=0.5*valA*inconnue(num_face,c2);
                                resu(j,c1)-=0.5*valA*inconnue(j,c2);
                              }
                          }

                    }

                }//loop over i0 (nb_faces) for elem1

              elem2 = face_voisins(num_face,1);
              for (i0=0; i0<nb_faces_elem; i0++)
                {
                  if ( ( (j= elemfaces(elem2,i0)) > num_face ) && (j!= fac_asso ) )
                    {
                      for(int c1=0; c1<nb_comp; c1++)

                        for(int c2=0; c2<nb_comp; c2++)
                          {

                            int diffusivity_index = c1*nb_comp + c2;
                            valA = viscA(num_face,j,elem2,nu(elem2, diffusivity_index));

                            resu(num_face,c1)+=valA*inconnue(j,c2);
                            resu(num_face,c1)-=valA*inconnue(num_face,c2);
                            if(j<nb_faces) // real face
                              {
                                resu(j,c1)+=0.5*valA*inconnue(num_face,c2);
                                resu(j,c1)-=0.5*valA*inconnue(j,c2);
                              }
                          }

                    }
                }//loop over i0 (nb_faces) for elem2
            }//loop over ind_face

        } // end of periodic case
      else
        {

          for (ind_face=num1; ind_face<num2; ind_face++)
            {
              num_face = le_bord.num_face(ind_face);
              elem1=face_voisins(num_face,0);

              // Loop over faces :
              for (int i=0; i<nb_faces_elem; i++)
                if (( (j= elemfaces(elem1,i)) > num_face ) || (ind_face>=nb_faces_bord_reel))
                  {
                    for(int c1=0; c1<nb_comp; c1++)


                      for(int c2=0; c2<nb_comp; c2++)
                        {
                          int diffusivity_index = c1*nb_comp + c2;
                          valA = viscA(num_face,j,elem1,nu(elem1, diffusivity_index));


                          if (ind_face<nb_faces_bord_reel)
                            {
                              double flux=valA*(inconnue(j,c2)-inconnue(num_face,c2));
                              resu(num_face,c1)+=flux;
                              tab_flux_bords(num_face,c1)+=flux;
                            }

                          if(j<nb_faces) // real faces
                            {
                              resu(j,c1)+=valA*inconnue(num_face,c2);
                              resu(j,c1)-=valA*inconnue(j,c2);
                            }
                        }
                  }
            }//end of ind_face loop

        }
    }//end for the boundaries loop

  // dealing now with internal faces


  for (num_face=zone_VEF.premiere_face_int(); num_face<nb_faces; num_face++)
    {
      for (int k=0; k<2; k++)
        {
          elem0 = face_voisins(num_face,k);
          for (i0=0; i0<nb_faces_elem; i0++)
            {
              if ( (j= elemfaces(elem0,i0)) > num_face )
                {

                  int el1,el2;
                  int contrib=1;
                  if(j>=nb_faces) // it is a virtual face
                    {
                      el1 = face_voisins(j,0);
                      el2 = face_voisins(j,1);
                      if((el1==-1)||(el2==-1))
                        contrib=0;
                    }
                  if(contrib)
                    {
                      for(int c1=0; c1<nb_comp; c1++)

                        for(int c2=0; c2<nb_comp; c2++)
                          {
                            int diffusivity_index = c1*nb_comp + c2;
                            valA = viscA(num_face,j,elem0,nu(elem0, diffusivity_index));
                            resu(num_face,c1)+=valA*inconnue(j,c2);
                            resu(num_face,c1)-=valA*inconnue(num_face,c2);
                            if(j<nb_faces) // dealing with real faces
                              {
                                resu(j,c1)+=valA*inconnue(num_face,c2);
                                resu(j,c1)-=valA*inconnue(j,c2);
                              }
                            else
                              {
                                //the face is virtual
                              }
                          }
                    }
                }
            }//loop i0 over elem faces
        }//loop over k
    }//loop over faces



  //Neumann :
  for (int n_bord=0; n_bord<nb_bords; n_bord++)
    {
      const Cond_lim& la_cl = zone_Cl_VEF.les_conditions_limites(n_bord);

      if (sub_type(Neumann_paroi,la_cl.valeur()))
        {
          const Neumann_paroi& la_cl_paroi = ref_cast(Neumann_paroi, la_cl.valeur());
          const Front_VF& le_bord = ref_cast(Front_VF,la_cl.frontiere_dis());
          int ndeb = le_bord.num_premiere_face();
          int nfin = ndeb + le_bord.nb_faces();
          for (int face=ndeb; face<nfin; face++)
            {
              for (int nc=0; nc<nb_comp; nc++)
                {
                  flux0=la_cl_paroi.flux_impose(face-ndeb,nc)*zone_VEF.surface(face);
                  resu(face,nc) += flux0;
                  tab_flux_bords(face,nc) = flux0;
                }
            }
        }
      else if (sub_type(Echange_externe_impose,la_cl.valeur()))
        {
          const Echange_externe_impose& la_cl_paroi = ref_cast(Echange_externe_impose, la_cl.valeur());
          const Front_VF& le_bord = ref_cast(Front_VF,la_cl.frontiere_dis());
          int ndeb = le_bord.num_premiere_face();
          int nfin = ndeb + le_bord.nb_faces();
          for (int face=ndeb; face<nfin; face++)
            {
              for (int nc=0; nc<nb_comp; nc++)
                {
                  flux0=la_cl_paroi.h_imp(face-ndeb,nc)*(la_cl_paroi.T_ext(face-ndeb,nc)-inconnue(face,nc))*zone_VEF.surface(face);
                  resu(face,nc) += flux0;
                  tab_flux_bords(face,nc) = flux0;
                }
            }
        }
      else if (sub_type(Neumann_homogene,la_cl.valeur())
               || sub_type(Symetrie,la_cl.valeur())
               || sub_type(Neumann_sortie_libre,la_cl.valeur()))
        {
          const Front_VF& le_bord = ref_cast(Front_VF,la_cl.frontiere_dis());
          int ndeb = le_bord.num_premiere_face();
          int nfin = ndeb + le_bord.nb_faces();
          for (int face=ndeb; face<nfin; face++)
            for (int nc=0; nc<nb_comp; nc++)
              tab_flux_bords(face,nc) = 0.;
        }
    }

}

DoubleTab& Op_Diff_VEF_Face_Matricial::calculer(const DoubleTab& inconnue, DoubleTab& resu) const
{
  resu = 0;
  return ajouter(inconnue,resu);
}



void Op_Diff_VEF_Face_Matricial::mettre_a_jour(double)
{
  remplir_nu(nu_);
}

void Op_Diff_VEF_Face_Matricial::contribuer_a_avec(const DoubleTab& transporte, Matrice_Morse& matrice) const
{
  remplir_nu(nu_);
  modifier_matrice_pour_periodique_avant_contribuer(matrice,equation());
  // On remplit le tableau nu car l'assemblage d'une
  // matrice avec ajouter_contribution peut se faire
  // avant le premier pas de temps

  const Zone_Cl_VEF& zone_Cl_VEF = la_zcl_vef.valeur();
  const Zone_VEF& zone_VEF = la_zone_vef.valeur();
  const IntTab& elem_faces = zone_VEF.elem_faces();
  const IntTab& face_voisins = zone_VEF.face_voisins();

  int n1 = zone_VEF.nb_faces();
  int nb_comp = 1;
  int nb_dim = transporte.nb_dim();

  DoubleTab nu;
  int marq=phi_psi_diffuse(equation());
  const DoubleVect& porosite_elem = zone_VEF.porosite_elem();

  // soit on a div(phi nu grad inco)
  // soit on a div(nu grad phi inco)
  // cela depend si on diffuse phi_psi ou psi
  modif_par_porosite_si_flag(nu_,nu,!marq,porosite_elem);
  DoubleVect porosite_eventuelle(zone_VEF.porosite_face());
  if (!marq)
    porosite_eventuelle=1;


  if(nb_dim==2)
    nb_comp=transporte.dimension(1);

  int i,j,num_face;
  int elem1,elem2;
  int nb_faces_elem = zone_VEF.zone().nb_faces_elem();
  double val;

  int nb_bords=zone_VEF.nb_front_Cl();
  for (int n_bord=0; n_bord<nb_bords; n_bord++)
    {
      const Cond_lim& la_cl = zone_Cl_VEF.les_conditions_limites(n_bord);
      const Front_VF& le_bord = ref_cast(Front_VF,la_cl.frontiere_dis());
      int num1 = le_bord.num_premiere_face();
      int num2 = num1 + le_bord.nb_faces();
      if (sub_type(Periodique,la_cl.valeur()))
        {
          const Periodique& la_cl_perio = ref_cast(Periodique,la_cl.valeur());
          int fac_asso;
          // on ne parcourt que la moitie des faces periodiques
          // on copiera a la fin le resultat dans la face associe..
          int num2b=num1+le_bord.nb_faces()/2;
          for (num_face=num1; num_face<num2b; num_face++)
            {
              elem1 = face_voisins(num_face,0);
              elem2 = face_voisins(num_face,1);
              fac_asso = la_cl_perio.face_associee(num_face-num1)+num1;
              for (i=0; i<nb_faces_elem; i++)
                {
                  if ( (j=elem_faces(elem1,i)) > num_face )
                    {
                      for(int c1=0; c1<nb_comp; c1++)
                        for(int c2=0; c2<nb_comp; c2++)
                          {

                            int diffusivity_index = c1*nb_comp + c2;
                            val = viscA(num_face,j,elem1,nu(elem1, diffusivity_index));
                            int n0=num_face*nb_comp+c1;
                            int j0=j*nb_comp+c1;
                            int n0b=num_face*nb_comp+c2;
                            int j0b=j*nb_comp+c2;

                            matrice(n0,n0b)+=val*porosite_eventuelle(num_face);
                            matrice(n0,j0b)-=val*porosite_eventuelle(j);
                            matrice(j0,n0b)-=val*porosite_eventuelle(num_face);
                            matrice(j0,j0b)+=val*porosite_eventuelle(j);


                          }
                    }
                  if (elem2!=-1)
                    if ( (j=elem_faces(elem2,i)) > num_face )
                      {
                        for (int nc=0; nc<nb_comp; nc++)
                          {
                            for(int c1=0; c1<nb_comp; c1++)
                              for(int c2=0; c2<nb_comp; c2++)
                                {

                                  int diffusivity_index = c1*nb_comp + c2;
                                  val = viscA(num_face,j,elem1,nu(elem1, diffusivity_index));
                                  int n0=num_face*nb_comp+c1;
                                  int j0=j*nb_comp+c1;
                                  int n0b=num_face*nb_comp+c2;
                                  int j0b=j*nb_comp+c2;

                                  matrice(n0,n0b)+=val*porosite_eventuelle(num_face);
                                  matrice(n0,j0b)-=val*porosite_eventuelle(j);
                                  matrice(j0,n0b)-=val*porosite_eventuelle(num_face);
                                  matrice(j0,j0b)+=val*porosite_eventuelle(j);

                                  abort();
                                  int  n0periob=fac_asso*nb_comp+c2;
                                  matrice(j0,n0periob)-=val*porosite_eventuelle(num_face);
                                }
                          }
                      }
                }
            }
        }
      else
        {
          for (num_face=num1; num_face<num2; num_face++)
            {
              elem1 = face_voisins(num_face,0);
              for (i=0; i<nb_faces_elem; i++)
                {
                  if ( (j= elem_faces(elem1,i)) > num_face )
                    {
                      for(int c1=0; c1<nb_comp; c1++)

                        for(int c2=0; c2<nb_comp; c2++)
                          {

                            int diffusivity_index = c1*nb_comp + c2;
                            val = viscA(num_face,j,elem1,nu(elem1, diffusivity_index));
                            int n0=num_face*nb_comp+c1;
                            int j0=j*nb_comp+c1;
                            int n0b=num_face*nb_comp+c2;
                            int j0b=j*nb_comp+c2;

                            matrice(n0,n0b)+=val*porosite_eventuelle(num_face);
                            matrice(n0,j0b)-=val*porosite_eventuelle(j);
                            matrice(j0,n0b)-=val*porosite_eventuelle(num_face);
                            matrice(j0,j0b)+=val*porosite_eventuelle(j);

                          }
                    }
                }
            }
        }
    }
  for (num_face=zone_VEF.premiere_face_int(); num_face<n1; num_face++)
    {
      elem1 = face_voisins(num_face,0);
      elem2 = face_voisins(num_face,1);

      for (i=0; i<nb_faces_elem; i++)
        {
          if ( (j=elem_faces(elem1,i)) > num_face )
            {
              for(int c1=0; c1<nb_comp; c1++)
                for(int c2=0; c2<nb_comp; c2++)
                  {
                    int diffusivity_index = c1*nb_comp + c2;
                    val = viscA(num_face,j,elem1,nu(elem1, diffusivity_index));
                    int n0=num_face*nb_comp+c1;
                    int j0=j*nb_comp+c1;
                    int n0b=num_face*nb_comp+c2;
                    int j0b=j*nb_comp+c2;

                    matrice(n0,n0b)+=val*porosite_eventuelle(num_face);
                    matrice(n0,j0b)-=val*porosite_eventuelle(j);
                    matrice(j0,n0b)-=val*porosite_eventuelle(num_face);
                    matrice(j0,j0b)+=val*porosite_eventuelle(j);


                  }
            }
          if (elem2!=-1)
            if ( (j=elem_faces(elem2,i)) > num_face )
              {
                for(int c1=0; c1<nb_comp; c1++)
                  for(int c2=0; c2<nb_comp; c2++)
                    {
                      int diffusivity_index = c1*nb_comp + c2;
                      val = viscA(num_face,j,elem2,nu(elem2, diffusivity_index));
                      int n0=num_face*nb_comp+c1;
                      int j0=j*nb_comp+c1;
                      int n0b=num_face*nb_comp+c2;
                      int j0b=j*nb_comp+c2;

                      matrice(n0,n0b)+=val*porosite_eventuelle(num_face);
                      matrice(n0,j0b)-=val*porosite_eventuelle(j);
                      matrice(j0,n0b)-=val*porosite_eventuelle(num_face);
                      matrice(j0,j0b)+=val*porosite_eventuelle(j);

                    }
              }
        }
    }
  modifier_matrice_pour_periodique_apres_contribuer(matrice,equation());
  //  matrice.imprimer_formatte(Cout);
}
