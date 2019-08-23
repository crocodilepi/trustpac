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
/////////////////////////////////////////////////////////////////////////////
//
// File      : Loi_Fermeture_PEMFC_base.cpp
// Directory : $PEMFC_ROOT/src
//
/////////////////////////////////////////////////////////////////////////////

#include <Loi_Fermeture_PEMFC_base.h>
#include <Param.h>
#include <Discretisation_base.h>
#include <Champ_P1NC.h>
#include <Milieu_base.h>
#include <Probleme_base.h>
#include <Zone_VF.h>
#include <Zone_VEF.h>
#include <Debog.h>
#include <Check_espace_virtuel.h>
#include <Fluide_Quasi_Compressible.h>
#include <Loi_Etat_Melange_GP_Fraction_Molaire.h>
#include <Neumann_paroi.h>
#include <Scalaire_impose_paroi.h>
#include <Interprete.h>
#include <Champ_Uniforme.h>
#include <Champ_front_calc.h>
#include <Champ_front_var.h>
#include <Champ_front.h>
#include <Domaine.h>

Implemente_base( Loi_Fermeture_PEMFC_base, "Loi_Fermeture_PEMFC_base", Loi_Fermeture_base ) ;

Sortie& Loi_Fermeture_PEMFC_base::printOn( Sortie& os ) const
{
  Loi_Fermeture_base::printOn( os );
  return os;
}

Entree& Loi_Fermeture_PEMFC_base::readOn( Entree& is )
{
  is_cas4_=0;
  Loi_Fermeture_base::readOn( is );
  return is;
}

void Loi_Fermeture_PEMFC_base::discretiser(const Discretisation_base& dis)

{
  Loi_Fermeture_base::discretiser(dis);
}

#include <Linear_algebra_tools_impl.h>

// evaluer la matrice de diffusion (cas2)
// quang: 16/06/19 modif sur la formulation qui peut fonctionner avec X2 = gaz H2 ou O2
double Loi_Fermeture_PEMFC_base::eval_matrice_diffusion_cas2(const double& cX2,const double& cvap,const double& cN2, DoubleTab& Diff,int elem)
{
  double T = T_(elem);

  double Rp = ch_Rp_.valeurs()(elem,0);
  double epsilon=ch_por_.valeurs()(elem,0);
  double tau= ch_tor_.valeurs()(elem,0);
  double K=ch_K_.valeurs()(elem,0);

  double  MX2=Mi_[0];
  double  Mvap=Mi_[1];
  double  MN2=Mi_[2];

  assert(MX2==32.e-3 || MX2==2.e-3);
  assert(Mvap==18.e-3);
  assert(MN2==28.e-3);

  const Fluide_Incompressible& le_fluide =ref_cast(Fluide_Incompressible,mon_probleme().milieu());
  const DoubleTab& mutab = le_fluide.viscosite_dynamique().valeurs();
  double mu = mutab(0,0);
  if(mutab.dimension(0)>1)
    mu = mutab(elem,0);

  double  Krg=1;						// permeabilite relative de gaz = const ??

  double  cg = cX2 + cN2 + cvap;
  double  Xvap =cvap/cg ;
  double  XX2=cX2/cg;
  double  XN2=cN2/cg;
  double  Pg=cg*R*T;

  // Darcy Binar diffusion coefficient [m2/s]
  double DX2N2 = 0., DX2vap = 0.;
  if(MX2==32.e-3) 	// gas O2
    {
      DX2N2 = f_D_O2_N2(T,Pg);
      DX2vap = f_D_O2_vap(T,Pg);
    }
  else if(MX2==2.e-3)    // gas H2
    {
      DX2N2 = f_D_H2_N2(T,Pg);
      DX2vap = f_D_H2_vap(T,Pg);
    }
  else
    {
      Cerr << "Unknown the gas" << MX2 << finl;
    }
  double DN2vap = f_D_N2_vap(T, Pg);
  double DN2X2=DX2N2;
  double DvapX2=DX2vap;
  double DvapN2=DN2vap;

  // Knudsen diffusion coefficient [m2/s]
  double  DKX2=2./3.*Rp*sqrt(8*R*T/(pi*MX2));
  double  DKvap=2./3.*Rp*sqrt(8*R*T/(pi*Mvap));
  double  DKN2=2./3.*Rp*sqrt(8*R*T/(pi*MN2));

  double  invDAX2vap=1./DX2vap+1./DKX2;
  double  invDAX2N2=1./DX2N2+1./DKX2;
  double  invDAvapX2=1./DvapX2+1./DKvap;
  double  invDAvapN2=1./DvapN2+1./DKvap;
  double  invDAN2X2=1./DN2X2+1./DKN2;
  double  invDAN2vap=1./DN2vap+1./DKN2;

  double  AK =0.75/Rp*sqrt(pi*R*T*0.5);
  double  ss=(XX2*sqrt(MX2)+Xvap*sqrt(Mvap)+XN2*sqrt(MN2));
  double  Ac=epsilon*mu/(cg*tau*tau*K*Krg*ss);
  double  AA=1./(1./Ac+1./AK);

  // Matrice
  double a11=-XN2*invDAX2N2-Xvap*invDAX2vap;
  double a12=XX2*invDAvapX2;
  double a13=XX2*invDAN2X2;

  double a21=Xvap*invDAX2vap;
  double a22=-XX2*invDAvapX2-XN2*invDAvapN2;
  double a23=Xvap*invDAN2vap;

  double a31= -AA*sqrt(MX2);
  double a32= -AA*sqrt(Mvap);
  double a33= -AA*sqrt(MN2);

  Matrice33 M(a11,a12,a13,
              a21,a22,a23,
              a31,a32,a33);
  Matrice33 invM;
  Matrice33::inverse(M,invM);

  Matrice33 coef(1-XX2,-XX2,-XX2,
                 -Xvap,1-Xvap,-Xvap,
                 R*T,R*T,R*T);								// VERIFIER

  for (int i=0; i<3; i++)
    for (int j=0; j<3; j++)
      {
        double p=0;
        for (int k=0; k<3; k++)
          p+=invM(i,k)*coef(k,j);
        p*=epsilon/(tau*tau);				// *epsilon_sur_tau2 OK
        Diff(elem,i*3+j)=-p;
      }

  return 1.;
}

double Loi_Fermeture_PEMFC_base::eval_matrice_diffusion_cas4(const double& cX2,const double& cvap,const double& cN2, DoubleTab& Diff,int elem)
{
  double T = T_(elem);
  //double eps_sur_tau2 = eps_sur_tau2_.valeurs()(elem,0);
  double eps_sur_tau2 = 1.;
  if(ch_por_.non_nul() && ch_tor_.non_nul())
    {
      double eps = ch_por_.valeurs()(0,0);
      if(!sub_type(Champ_Uniforme, ch_por_.valeur()))
        eps = ch_por_.valeurs()(elem,0);

      double tau = ch_tor_.valeurs()(0,0);
      if(!sub_type(Champ_Uniforme, ch_por_.valeur()))
        tau = ch_tor_.valeurs()(elem,0);

      eps_sur_tau2 *= eps;
      eps_sur_tau2 /= tau*tau;
    }
  else
    {
      Cerr << "Loi_Fermeture_PEMFC_Cas4 needs to input the tortuosity par keyword 'tor' and the porosity par keyword 'por'";
      abort();
    }

  const Fluide_Quasi_Compressible& le_fluideQC=ref_cast(Fluide_Quasi_Compressible,mon_probleme().milieu());
  const Loi_Etat_Melange_GP_Fraction_Molaire& loi_etat = ref_cast(Loi_Etat_Melange_GP_Fraction_Molaire,le_fluideQC.loi_etat().valeur());
  const ArrOfDouble& mi=loi_etat.get_masses_molaires();

  double  MX2=mi[0];
  assert(MX2==32.e-3 || MX2==2.e-3);
  double  cg = cX2 + cN2 + cvap;
  double  Xvap =cvap/cg ;
  double  XX2=cX2/cg;
  double  XN2=cN2/cg;
  double  Pg=cg*R*T;

  double DX2N2 = 0., DX2vap = 0.;
  if(MX2==32.e-3) 	// gas O2
    {
      DX2N2 = f_D_O2_N2(T,Pg);
      DX2vap = f_D_O2_vap(T,Pg);
    }
  else if(MX2==2.e-3)    // gas H2
    {
      DX2N2 = f_D_H2_N2(T,Pg);
      DX2vap = f_D_H2_vap(T,Pg);
    }
  else
    {
      Cerr << "Unknown the gas" << MX2 << finl;
    }
  double DN2vap = f_D_N2_vap(T, Pg);
  double DN2X2=DX2N2;
  double DvapX2=DX2vap;
  double DvapN2=DN2vap;

  double invDAX2vap=1./DX2vap;
  double invDAX2N2=1./DX2N2;
  double invDAvapX2=1./DvapX2;
  double invDAvapN2=1./DvapN2;
  double invDAN2X2=1./DN2X2;
  double invDAN2vap=1./DN2vap;

  // Matrice
  double a11=-XN2*invDAX2N2-Xvap*invDAX2vap;
  double a12=XX2*invDAvapX2;
  double a13=XX2*invDAN2X2;

  double a21=Xvap*invDAX2vap;
  double a22=-XX2*invDAvapX2-XN2*invDAvapN2;
  double a23=Xvap*invDAN2vap;

  Matrice33 invM;
  double a31= 1;
  double a32= 1;
  double a33=1;

  Matrice33 M(a11,a12,a13,
              a21,a22,a23,
              a31,a32,a33);
  Matrice33::inverse(M,invM);
  /*
  Matrice33 coef(eps_sur_tau2/(R*T),0,0,
                 0,eps_sur_tau2/(R*T),0,
                 0,0,0);
  */

  Matrice33 coef(eps_sur_tau2*cg,0,0,
                 0,eps_sur_tau2*cg,0,
                 0,0,0);

  for (int i=0; i<3; i++)
    for (int j=0; j<3; j++)
      {
        double p=0;
        for (int k=0; k<3; k++)
          p+=invM(i,k)*coef(k,j);
        Diff(elem,i*3+j)=-p;
      }

  return 1.;
}

// quang: cette methode a pour but de calculer la matrice de diffusition nu OK
void Loi_Fermeture_PEMFC_base::bidouille_nu(DoubleTab& nu,const DoubleTab&   inconnue_org,const Zone_VF& zone_VF)
{
  int nb_comp=inconnue_org.dimension(1);
  assert(nb_comp==3);
  ArrOfDouble c(3);

  const IntTab& elem_faces=zone_VF.elem_faces();
  int nbf=elem_faces.dimension(1);

  double invnbf=1./nbf;
  int nb_elem=zone_VF.zone().nb_elem();
  assert(nb_elem==nu.dimension(0));
  assert(zone_VF.nb_faces()==inconnue_org.dimension(0));
  assert(elem_faces.dimension(0)==nb_elem);

  if(is_cas4_) 		// inconnu = fraction_molaire Xi
    {
      const DoubleTab& cg=mon_probleme().get_champ("cg").valeurs();
      assert(cg.dimension(0)==inconnue_org.dimension(0));
      for (int elem=0; elem<nb_elem; elem++)
        {
          c=0;
          for (int nc=0; nc<nb_comp; nc++)
            for (int f=0; f<nbf; f++)
              {
                int face=elem_faces(elem,f);
                c(nc)+=inconnue_org(face,nc)*cg(face);		// X = c/cg -> c = X.cg
              }
          c*=invnbf;										// moyen de concentration = sum(c_face)/nb_faces ?? -> interpolation P1NC vers P0 ?
          eval_matrice_diffusion_cas4(c[0],c[1],c[2],  nu, elem);
        }
      nu.echange_espace_virtuel();
    }
  else			// cas2 inconnu = concentration
    {
      for (int elem=0; elem<nb_elem; elem++)
        {
          c=0;
          for (int nc=0; nc<nb_comp; nc++)
            for (int f=0; f<nbf; f++)
              {
                int face=elem_faces(elem,f);
                c(nc)+=inconnue_org(face,nc);
              }
          c*=invnbf;			// moyen de concentration = sum(c_face)/nb_faces ?? -> interpolation P1NC vers P0 ?
          eval_matrice_diffusion_cas2(c[0],c[1],c[2],  nu, elem);
        }
      nu.echange_espace_virtuel();
    }


  // for only debug
  //double maxc=local_max_vect(nu);
  //if (local_min_vect(nu)<0)
  // {
  //Cerr<<"oooooo !!!!!!!!!!!! "<< local_min_vect(nu) <<finl;
  //abort();
  //}
  //Cerr<<"iii "<<nu(0,0)<<finl;
  //Cerr<<"max nu "<<maxc<<finl;
}

// copy of discretiation_tools
void Loi_Fermeture_PEMFC_base::discretisation_tools_cells_to_faces(const Champ_base& He,  Champ_base& Hf)
{
  DoubleTab& tabHf=Hf.valeurs();
  const DoubleTab& tabHe=He.valeurs();
  Debog::verifier("element_face entreee",tabHe);
  assert_espace_virtuel_vect(tabHe);
  const Zone_dis_base& zone_dis_base=He.zone_dis_base();

  const Zone_VF& zone_vf= ref_cast(Zone_VF,zone_dis_base);
  // en realite on fait P1B vers face
  //assert(tabHe.dimension_tot(0)==zone_dis_base.nb_elem_tot());
  assert(tabHf.dimension_tot(0)==zone_vf.nb_faces_tot());

  const IntTab& elem_faces=zone_vf.elem_faces();
  const DoubleVect& volumes=zone_vf.volumes();
  const DoubleVect& volumes_entrelaces=zone_vf.volumes_entrelaces();



  tabHf=0;
  int nb_face_elem=elem_faces.dimension(1);
  int nb_elem_tot=zone_dis_base.nb_elem_tot();
  int nb_dim=tabHe.nb_dim();


  double coeffb=nb_face_elem;
  double coeffi=coeffb;
  if (zone_vf.que_suis_je()=="Zone_VDF")
    {
      coeffb=1;
      coeffi=2;
    }

  if (nb_dim==1)
    {
      for (int ele=0; ele<nb_elem_tot; ele++)
        {
          for (int s=0; s<nb_face_elem; s++)
            {

              tabHf(elem_faces(ele,s))+=tabHe(ele)*volumes(ele);
              //	      Cerr<<elem_faces(ele,s)<<" "<<tabHe(ele)*volumes(ele)<<" "<<tabHf(elem_faces(ele,s))<<finl;
            }
        }
      for (int f=0; f<zone_vf.premiere_face_int(); f++)
        tabHf(f)/=volumes_entrelaces(f)*coeffb;
      for (int f=zone_vf.premiere_face_int(); f<zone_vf.nb_faces(); f++)
        tabHf(f)/=volumes_entrelaces(f)*coeffi;

    }
  else
    {

      if (tabHf.nb_dim()==1)
        {
          assert(coeffi==2);
          assert(coeffb==1);
          // VDF
          //abort();
          for (int ele=0; ele<nb_elem_tot; ele++)
            {
              for (int s=0; s<nb_face_elem; s++)
                {
                  int face=elem_faces(ele,s);
                  //for (int comp=0;comp<nb_comp;comp++)
                  int comp=zone_vf.orientation()[face];
                  tabHf(face)+=tabHe(ele,comp)*volumes(ele);
                }
            }
          for (int f=0; f<zone_vf.premiere_face_int(); f++)
            tabHf(f)/=volumes_entrelaces(f)*coeffb;
          for (int f=zone_vf.premiere_face_int(); f<zone_vf.nb_faces(); f++)
            tabHf(f)/=volumes_entrelaces(f)*coeffi;
        }
      else
        {
          //abort();
          int nb_comp=tabHf.dimension(1);
          for (int ele=0; ele<nb_elem_tot; ele++)
            {
              for (int s=0; s<nb_face_elem; s++)
                {
                  int face=elem_faces(ele,s);
                  for (int comp=0; comp<nb_comp; comp++)
                    tabHf(face,comp)+=tabHe(ele,comp)*volumes(ele);
                }
            }

          for (int f=0; f<zone_vf.premiere_face_int(); f++)
            for (int comp=0; comp<nb_comp; comp++)
              tabHf(f,comp)/=volumes_entrelaces(f)*coeffb;
          for (int f=zone_vf.premiere_face_int(); f<zone_vf.nb_faces(); f++)
            for (int comp=0; comp<nb_comp; comp++)
              tabHf(f,comp)/=volumes_entrelaces(f)*coeffi;

        }
    }
  tabHf.echange_espace_virtuel();
//  Cerr<<min_array(tabHe)<<" elem  "<<max_array(tabHe)<<finl;
//  Cerr<<min_array(tabHf)<<" face  "<<max_array(tabHf)<<finl;

  Debog::verifier("element_face sortie",tabHf);
}


// mettre a jour Ni et ud (Stefan Maxwell dans CH) cas4 + concentration "ci"
void Loi_Fermeture_PEMFC_base::calculer_Ni_ud(const double& temps)
{
//  Cerr << "Loi_Fermeture_PEMFC_base::calculer_Ni_ud Cas4 au temps=" << temps << finl;
//  Cerr << "DEBUG valeur champ_Inc avant MAJ" << finl;
//  Cerr << "DEBUG ci_CH avant MAJ temps_courant=" << ci_ch_.temps() << " "<< ci_ch_.valeurs() << finl;
//  Cerr << "DEBUG Ns_CH avant MAJ temps_courant=" << Ns_ch_.temps() << " "<< Ns_ch_.valeurs() << finl;

  const  DoubleTab& Xi= equation().inconnue().valeurs();		// note fraction molaire Xi (cas4)
  int nb_comp=Xi.dimension(1);									// must be equal to 3

  assert(nb_comp == 3);

  // calcul du gradient grad(Xi)
  DoubleTab grad(0, nb_comp, Objet_U::dimension);
  equation().zone_dis().zone().creer_tableau_elements(grad);

  // flux diffusive Ni_d = -alpha.grad(Xi)
  DoubleTab& Ni=Ni_.valeurs();

  // vitesse diffusion (P0) ud = sum(Mi Ni_d)
  DoubleTab& ude=ud_.valeurs();

  const Champ_base& ch_rho_face=mon_probleme().get_champ("masse_volumique");
  const DoubleTab& rho_face=ch_rho_face.valeurs();						// note rhog
  assert(rho_face.get_md_vector()==Um_.valeurs().get_md_vector());

  const Champ_base& vitesse=mon_probleme().get_champ("vitesse");		// note vitesse massique 'ug'
  const DoubleTab& cg=mon_probleme().get_champ("cg").valeurs();			// note concentration totale 'cg'

  const Fluide_Quasi_Compressible& le_fluideQC=ref_cast(Fluide_Quasi_Compressible,mon_probleme().milieu());
  const Loi_Etat_Melange_GP_Fraction_Molaire& loi_etat = ref_cast(Loi_Etat_Melange_GP_Fraction_Molaire,le_fluideQC.loi_etat().valeur());
  const ArrOfDouble& mi=loi_etat.get_masses_molaires();					// note Mi

  const DoubleTab& nu=diffu_.valeurs();									// note Dx

  Debog::verifier("loi nu",nu);
  Debog::verifier("loi inco",Xi);

  // debug masse volumique -> not changed ??? masse_volumique needs to be updated rhog = cg.sum(Xi.Mi)
  double max_masse = mp_max_abs_vect(rho_face);
  double min_masse = mp_min_abs_vect(rho_face);
//  for (int face = 1; face < rho_face.dimension(0); ++face)
//    {
//      double drho = rho_face(face) - rho_face(face-1);
//      if(abs(drho)>max_masse) max_masse = abs(drho);
//    }
  Cerr << "Max masse [kg/m3] " << max_masse << finl;
  Cerr << "Min masse [kg/m3] " << min_masse << finl;
  // end debug masse volumique

  // debug flux_bords
  DoubleTab& flux_diff = equation().operateur(0).l_op_base().flux_bords();
  DoubleTab& flux_conv = equation().operateur(1).l_op_base().flux_bords();
  Cerr << "DEBUG: Flux diffusion : nb_faces=" << flux_diff.dimension(0) << ",nb_comp= " << flux_diff.dimension(1) << finl;
  if(flux_conv.dimension(0) != 0) Cerr << "DEBUG: Flux convection: nb_faces=" << flux_conv.dimension(0) << ",nb_comp= " << flux_conv.dimension(1) << finl;
  // end debug flux_bords

  // mettre a jour ci = cg*Xi
  DoubleTab& ci = ci_ch_.valeurs();
  for (int face = 0; face < ci.dimension(0); ++face)
    {
      for (int j = 0; j < nb_comp; ++j)
        {
          ci(face,j) = cg(face)*Xi(face,j);
        }
    }
  ci_ch_.mettre_a_jour(temps);

  /*
  Champ_P1NC::calcul_gradient(inco,grad,ref_cast(Zone_Cl_VEF,equation().zone_Cl_dis().valeur()));		// gradXi
  Debog::verifier("loi grad",grad);
  grad.echange_espace_virtuel();
  Debog::verifier("loi grad",grad);
  for (int elem=0; elem<grad.dimension(0); elem++)
    {
      for (int j=0; j<dimension; j++)
        {
          ude(elem,j)=0;
          for (int i=0; i<nb_comp; i++)
            {
              double p=0;
              for (int k=0; k<nb_comp; k++)
                p+=nu(elem,i*nb_comp+k)*grad(elem,k,j);

              Ni(elem,i*dimension+j)=-p;				// Ni = -nu.gradXi

              ude(elem,j)-=p*mi[i] ;					// ud = sum(-nu.gradXi.Mi) = sum(Ni.Mi)
            }
          // ude(elem,j)/=rho_elem(elem);
        }

      // DVQ comment: this for checking sum(flux_passed_elem) = 0 for every element
      for (int j=0; j<dimension; j++)
        {
          double pt=0;
          for (int i=0; i<nb_comp; i++)
            pt+=Ni(elem,i*dimension+j);
          if (!est_egal(pt,0))
            {
              Cerr<<"Ni oooo "<< pt<<" "<<Ni(elem,0*dimension+j)<<" "<<Ni(elem,1*dimension+j)<<" "<<Ni(elem,2*dimension+j)<< finl;
              //exit();
            }

        }
    }
  ude.echange_espace_virtuel();
  Ni.echange_espace_virtuel();
  */
  //Champ_P1NC::calcul_gradient(ci,grad,ref_cast(Zone_Cl_VEF,equation().zone_Cl_dis().valeur()));		// gradci

  // debug gradXi
  Champ_P1NC::calcul_gradient(Xi,grad,ref_cast(Zone_Cl_VEF,equation().zone_Cl_dis().valeur()));		// gradXi
  Debog::verifier("loi grad",grad);
  grad.echange_espace_virtuel();
  Debog::verifier("loi grad",grad);

  // test methode est_egal
  for (int elem=0; elem<grad.dimension(0); elem++)
    {
      for (int j=0; j<dimension; j++)
        {
          ude(elem,j)=0;
          for (int i=0; i<nb_comp; i++)
            {
              double p=0;
              for (int k=0; k<nb_comp; k++)
                p+=nu(elem,i*nb_comp+k)*grad(elem,k,j);					// Dx.grad(Xi)

              Ni(elem,i*dimension+j)=-p;				// Ni = -Dx.gradXi
              ude(elem,j)-=p*mi[i] ;					// ud = sum(-nu.gradci.Mi) = sum(Ni.Mi)

              // DEBUG Ni ud
//              if(elem == 10) Cerr << "DEBUG Ni ud: elem=" << elem << ", compo i=" << i << ", dimension j=" << j << ", Ni(i,j)=" << Ni(elem,i*dimension+j) << ", ud(j)=" << ude(elem,j) << finl;
              // End DEBUG Ni ud
            }
        }

      // DVQ comment: this for checking sum(Ni) = 0 par composant
      for (int j=0; j<dimension; j++)
        {
          double pt=0;
          for (int i=0; i<nb_comp; i++)
            pt+=Ni(elem,i*dimension+j);

          if (!est_egal(pt,0))
            {
              Cerr<<"Ni oooo "<< pt<<" "<<Ni(elem,0*dimension+j)<<" "<<Ni(elem,1*dimension+j)<<" "<<Ni(elem,2*dimension+j)<< finl;
              exit();
            }
        }
    }
  ude.echange_espace_virtuel();
  Ni.echange_espace_virtuel();

  // MAJ Um = ug - ud
  discretisation_tools_cells_to_faces(ud_,Um_);						// P0 -> P1NC
  Um_.valeurs()*=-1;												// Um = -ud = -sum(Mi.Ni)
  tab_divide_any_shape(Um_.valeurs(),rho_face);						// Um = -ud/rhog = -sum(Mi.Ni)/rhog
  Um_.valeurs()+=vitesse.valeurs();									// Um = +ug-ud/rhog = +ug-sum(Mi.Ni)/rhog

  // assert( nb_comp>1 );
  const Zone_VF& zone_VF = ref_cast(Zone_VF,equation().zone_dis().valeur());
  int nb_bords= equation().zone_Cl_dis().nb_cond_lim();

  int ok=0;
  for (int i_bord=0; i_bord<nb_bords; i_bord++)
    {
      const Cond_lim& la_cl_sur_X = equation().zone_Cl_dis().les_conditions_limites(i_bord);

      // quang: note : condition limite ug_CH = ug_GDL se traduit par Um = -sum(flux_surfacique_passant/cg)*normale_face
      if (sub_type(Neumann_paroi,la_cl_sur_X.valeur()))
        {
          const Neumann_paroi& la_cl_flux=ref_cast(Neumann_paroi,la_cl_sur_X.valeur());
          const Nom& nom_bord =la_cl_sur_X.frontiere_dis().le_nom();
          Cerr <<"on calcule VGDL sur " <<nom_bord<<finl;

          // attention: On recupere la condition limite de vitesse ici
          Cond_lim_base& la_cl = ref_equation_.valeur().probleme().equation(0).zone_Cl_dis().valeur().condition_limite_de_la_frontiere(nom_bord);
          const Front_VF& le_bord = ref_cast(Front_VF,la_cl.frontiere_dis());

          // quang: ajout initialiser(temps)
          //if(temps==0) la_cl.initialiser(temps); // TRUST erreur: distributed array cannot be resized

          // attention: valeur de vitesse imposee
          DoubleTab& val_ug= la_cl.champ_front().valeur().valeurs_au_temps(temps);

          //          const Champ_front& ch_front_flux = la_cl_flux.champ_front();
          //          if(sub_type(Champ_front_calc, ch_front_flux.valeur()))
          //            {
          //              Cerr << "Debug: un champ_front_calc needs to be initialized!!" << finl;
          //              Champ_front_calc ch_fr_cl = ref_cast(Champ_front_calc, ch_front_flux.valeur());
          //              val = ch_fr_cl.valeurs_au_temps(temps);
          //            }

          // debug vitesse imposee
          assert(val_ug.dimension(0)==le_bord.nb_faces());
          Cerr << "DEBUG vitesse ug_CH impose au temps=" << temps << " sur le bord bord=" << le_bord.le_nom() << finl;
          for (int ind_face=0; ind_face < le_bord.nb_faces(); ind_face++)
            {
              int face = le_bord.num_face(ind_face);
              if(dimension==2)
                Cerr << "DEBUG ug_CH: , face=" << face << ", ind_face_bord=" << ind_face << ", surface=" << zone_VF.face_surfaces(face) << ", ug_CH(face)=" << val_ug(ind_face,0) << " " << val_ug(ind_face,1) << finl;
              else
                Cerr << "DEBUG ug_CH: , face=" << face << ", ind_face_bord=" << ind_face << ", surface=" << zone_VF.face_surfaces(face) << ", ug_CH(face)=" << val_ug(ind_face,0) << " " << val_ug(ind_face,1) << " " << val_ug(ind_face,2) << finl;
            }
          // end vitesse imposee


          ok=1;
          int num1 = 0;
          int num2 = le_bord.nb_faces_tot();
          //int nb_faces_bord_reel = le_bord.nb_faces();
          Cerr << "DEBUG flux surfacique diffusive impose au temps=" << temps << " sur le bord bord=" << le_bord.le_nom() << finl;
          const DoubleTab& val_flux= la_cl_flux.champ_front().valeur().valeurs_au_temps(temps);

          for (int ind_face=num1; ind_face<num2; ind_face++)
            {
              int num_face = le_bord.num_face(ind_face);
              double flux=0;
              for (int nc=0; nc<nb_comp; nc++)
                //flux+=la_cl_flux.flux_impose(ind_face,nc);					// sum (Ns)
                flux+=val_flux(ind_face,nc);
              flux/=cg(num_face)*zone_VF.face_surfaces(num_face);			// sum(Ns)/cg
              for (int dir=0 ; dir<dimension; dir++)
                {
                  double vgdl=-flux*zone_VF.face_normales(num_face,dir); // - car normale sortante
                  //val(ind_face,dir)=(vitesse.valeurs()(num_face,dir)-(Um_(num_face,dir)-vgdl));		// ??? DVQ: pas compris, pourquoi changer la valeur de cond lim
                  //Um_(num_face,dir)=vgdl;
                  val_ug(ind_face,dir)=vgdl;														// ug = sum(Ns)/cg
                  Um_(num_face,dir)=vgdl - (vitesse.valeurs()(num_face,dir)-Um_(num_face,dir));	// Um = vgdl - ud
                }
              // debug flux impose
              //Cerr << "DEBUG Paroi_flux_impose (CH): , face=" << num_face << ", ind_face_bord=" << ind_face << ", surface=" << zone_VF.face_surfaces(num_face) << ", flux_impose(face)=" << la_cl_flux.flux_impose(ind_face,0) << " " << la_cl_flux.flux_impose(ind_face,1) << " " << la_cl_flux.flux_impose(ind_face,2) << finl;
              Cerr << "DEBUG Paroi_flux_impose (CH): , face=" << num_face << ", ind_face_bord=" << ind_face << ", surface=" << zone_VF.face_surfaces(num_face) << ", flux_impose(face)=" << val_flux(ind_face,0) << " " << val_flux(ind_face,1) << " " << val_flux(ind_face,2) << finl;
              // debug flux impose
            }
          val_ug.echange_espace_virtuel();
        }
    }

  if (!ok) abort();

  Debog::verifier("loi Ni",Ni_);
  Debog::verifier("loi ud",ud_);
  Debog::verifier("loi Um",Um_);
  Ni_.changer_temps(temps);
  ud_.changer_temps(temps);
  Um_.changer_temps(temps);

  // mettre a jour flux surfacique Ns = - integral((Nd_i + ciUm),norma) // - car normale sortante ? -> a verifier
  DoubleTab& Ns = Ns_ch_.valeurs();
  DoubleTab& Um = Um_.valeurs();

  const Zone_VEF& zone_VEF = ref_cast(Zone_VEF,equation().zone_dis().valeur());
  const IntTab& face_voisins = zone_VEF.face_voisins();
  const DoubleTab& face_norm = zone_VEF.face_normales();		// attention face_norm_i = |S_i|.n_i

  // Debug ci sur frontiere
  Cerr << "DEBUG ci dom=" << zone_VEF.domaine_dis().domaine().le_nom() << ", nombre de bords=" << zone_VEF.nb_front_Cl() << finl;
  for (int i = 0; i < zone_VEF.nb_front_Cl(); ++i)
    {
      const Front_VF& bord = zone_VEF.front_VF(i);
      if(bord.le_nom()=="I_CH_GDL")
        {
          Cerr << "DEBUG ci au temps=" << temps << ", i=" << i << ", bord=" << bord.le_nom() << finl;
          int nfaces = bord.nb_faces();
          for (int ind_face=0; ind_face < nfaces; ind_face++)
            {
              int face = bord.num_face(ind_face);
              Cerr << "DEBUG ci: , face=" << face << ", ind_face_bord=" << ind_face << ", surface=" << zone_VEF.face_surfaces(face) << ", ci(face)=" << ci(face,0) << " " << ci(face,1) << " " << ci(face,2) << finl;
            }
        }
    }
  // End Debug ci sur frontiere

  // debug check face_normale=norma.Surface -> done OK
//  double max_norma = 0;
//  for (int face = 0; face < ci.dimension(0); ++face)
//    {
//      double norm2 = face_norm(face,0)*face_norm(face,0)+face_norm(face,1)*face_norm(face,1);
//      if(dimension == 3)
//        norm2 += face_norm(face,2)*face_norm(face,2);
//      double norm = sqrt(norm2);
//      if(max_norma < norm) max_norma = norm;
//      if(abs(norm-1.)>1e-10)
//        {
//          Cerr << "DEBUG: face_norm n'est pas un vecteur unitaire" << ", face_ID "<< face <<  ", norm " << norm << ", norm/surface "<< norm/zone_VEF.face_surfaces(face) << finl;
//        }
//    }
//  Cerr << "DEBUG: max face_norm" << max_norma << finl;
  // end debug face_normale

  for (int face = 0; face < ci.dimension(0); ++face)
    {
      int n0 = face_voisins(face,0);
      int n1 = face_voisins(face,1);
      double face_surf = zone_VEF.face_surfaces(face);			// surface |S_i|
//      if(dimension==2)
//        Cerr << "DEBUG Ns_CH: face=" << face << ", face_norm=" << face_norm(face,0) << ", " << face_norm(face,1) << ", n0=" << n0 << ", n1=" << n1;
//      else
//        Cerr << "DEBUG Ns_CH: face=" << face << ", face_norm=" << face_norm(face,0) << ", " << face_norm(face,1) << ", " << face_norm(face,2) << ", n0=" << n0 << ", n1=" << n1;
//
//      Cerr << ", Ni=[";
      for (int i = 0; i < nb_comp; ++i)
        {
          double res = 0;
          double cif = ci(face,i);
//          Cerr << "[";
          for (int j = 0; j < dimension; ++j)
            {
              double Niface = 0;
              if(n0 != -1 && n1 != -1)
                {
                  // Ni moyen de deux elements
                  Niface = 0.5 * (Ni(n0,i*dimension+j)+Ni(n1,i*dimension+j));
                }
              else if(n0 != -1)
                {
                  // Ni de element nO
                  Niface = Ni(n0,i*dimension+j);
                }
              else
                {
                  // Ni de element n1
                  Niface = Ni(n1,i*dimension+j);
                }
//              Cerr << Niface << " ";
              double Ni_plus_ciUm = Niface+cif*Um(face,j);				// Ni +ciUm
              //res -= Ni_plus_ciUm*face_norm(face,j);					// produit scalaire -(Ni+ciUm)*norma*Surface
              res += Ni_plus_ciUm*face_norm(face,j);					// produit scalaire +(Ni+ciUm)*norma*Surface

              // Debug Ns -> done OK
              //if(face == 10)
              //{
              //Cerr << "DEBUG Ns_CH: face=" << face << ", face_norm(j)=" << face_norm(face,j) << ", n0=" << n0 << ", n1=" << n1 << ", compo i=" << i << ", dimension j=" << j << ", Ni(i,j)=" << Niface << ", ci=" << cif << ", Um(j)=" << Um(face,j) << ", Ni+ciUm=" << Ni_plus_ciUm <<", (Ni+ciUm).n.S="<< res << finl;
              //}
              // End Debug Ns

              // debug Ni
              //res += Niface*face_norm(face,j);					// produit scalaire +Ni*norma*Surface
              // End debug Ni
            }
//          Cerr << "]";
          Ns(face,i) = res / face_surf;
          // Debug Ns -> done OK
          //if(face == 10) Cerr << "DEBUG Ns_CH: face=" << face << ", surface=" << face_surf << ", i=" << i << ", Ns(i)=" << Ns(face,i) << finl;
          // End Debug Ns
        }
//      Cerr << "], ci=" << ci(face,0) << " " << ci(face,1) << " " << ci(face,2);
//      if(dimension==2)
//        Cerr << ", Um=" << Um(face,0) << " valeur" << Um(face,1);
//      else
//        Cerr << ", Um=" << Um(face,0) << " " << Um(face,1) << " " << Um(face,2);
//      Cerr << ", Ns=" << Ns(face,0) << " " << Ns(face,1) << " " << Ns(face,2) << finl;
    }

  Ns_ch_.mettre_a_jour(temps);

  // Debug Ns flux_bord -> done OK (compare with flux_bords post-traitement)
  int nb_front = zone_VEF.nb_front_Cl();
  Cerr << "DEBUG flux_bord dom=" << zone_VEF.domaine_dis().domaine().le_nom() << ", nombre de bords=" << nb_front << finl;
  for (int i = 0; i < nb_front; ++i)
    {
      const Front_VF& bord = zone_VEF.front_VF(i);
      Cerr << "DEBUG flux_bord au temps=" << temps << ", i=" << i << ", bord=" << bord.le_nom() << finl;
      if(bord.le_nom()=="I_CH_GDL")
        {
          int nfaces = bord.nb_faces();
          for (int ind_face=0; ind_face < nfaces; ind_face++)
            {
              int face = bord.num_face(ind_face);
              //if(dimension==2)
              Cerr << "DEBUG flux_bord: face=" << face << ", face_norm=" << face_norm(face,0) << " " << face_norm(face,1) <<  ", ind_face=" << ind_face << ", surface=" << zone_VEF.face_surfaces(face) << ", Ns(face)=" << Ns(face,0) << " " << Ns(face,1) << " " << Ns(face,2) << " " << finl;
              //else

            }
        }

    }

  // End Debug Ns flux_bord
}

// mettre a jour Ni ug (Stefan Maxwell) dans GDL (cas2) + fraction molaire Xi + flux surfacique de diffusion Ns
void Loi_Fermeture_PEMFC_base::calculer_Ni_ug(const double& temps)
{
//  Cerr << "Loi_Fermeture_PEMFC_base::calculer_Ni_ug Cas2 au temps=" << temps << finl;
//  Cerr << "DEBUG valeur champ_Inc avant MAJ" << finl;
//  Cerr << "DEBUG ug_GDL avant MAJ temps_courant=" << ug_gdl_.temps() << " " << ug_gdl_.valeurs() << finl;
//  Cerr << "DEBUG cg_GDL avant MAJ temps_courant=" << cg_gdl_.temps() << " " << cg_gdl_.valeurs() << finl;
//  Cerr << "DEBUG Xi_GDL avant MAJ temps_courant=" << Xi_gdl_.temps() << " " << Xi_gdl_.valeurs() << finl;
//  Cerr << "DEBUG Ns_GDL avant MAJ temps_courant=" << Ns_gdl_.temps() << " " << Ns_gdl_.valeurs() << finl;
  const  DoubleTab& ci= equation().inconnue().valeurs();			// note inco = concentration ci (cas2)
  int nb_comp=ci.dimension(1);

  // calcul du gradient
  DoubleTab grad(0, nb_comp, Objet_U::dimension);
  equation().zone_dis().zone().creer_tableau_elements(grad);

  DoubleTab& Ni=Ni_.valeurs();
  DoubleTab& ude=ud_.valeurs();

  const DoubleTab& nu=diffu_.valeurs();								// matrice de diffusion Dc

  Debog::verifier("loi nu",nu);
  Debog::verifier("loi inco",ci);

  Champ_P1NC::calcul_gradient(ci,grad,ref_cast(Zone_Cl_VEF,equation().zone_Cl_dis().valeur()));		// gradci
  Debog::verifier("loi grad",grad);
  grad.echange_espace_virtuel();
  Debog::verifier("loi grad",grad);

  for (int elem=0; elem<grad.dimension(0); elem++)
    {
      for (int j=0; j<dimension; j++)
        {
          ude(elem,j)=0;
          for (int i=0; i<nb_comp; i++)
            {
              double p=0;
              for (int k=0; k<nb_comp; k++)
                {
                  double nik = nu(elem,i*nb_comp+k);
                  double gradkj=grad(elem,k,j);
                  double nugrad=nik*gradkj;
                  //p+=nu(elem,i*nb_comp+k)*grad(elem,k,j);
                  p+= nugrad;
                }
              Ni(elem,i*dimension+j)=-p;				// Ni = -nu.grad(ci)
              ude(elem,j)-=p;							// ud = sum(-nu.gradci)

              // DEBUG Ni ud
              //if(elem == 10) Cerr << "DEBUG Ni ud (GDL): elem=" << elem << ", compo i=" << i << ", dimension j=" << j << ", Ni(i,j)=" << Ni(elem,i*dimension+j) << ", ud(j)=" << ude(elem,j) << finl;
              // End DEBUG Ni ud
            }

        }
    }
  Ni.echange_espace_virtuel();
  Ni_.mettre_a_jour(temps);
  ud_.mettre_a_jour(temps);

  // mettre a jour cg = sum(ci)
  DoubleTab& cg = cg_gdl_.valeurs();
  for (int face = 0; face < cg.dimension(0); ++face)
    {
      double cgf = 0;
      for (int j = 0; j < nb_comp; ++j)
        {
          cgf += ci(face,j);
        }
      cg(face) = cgf;
    }
  cg_gdl_.mettre_a_jour(temps);

  // here we calculate ug = sum(Ni)/cg
  discretisation_tools_cells_to_faces(ud_,ug_gdl_);			// ud=sum(-nu.gradci) P0 -> P1NC
  for(int face=0; face < cg.dimension(0); face++)
    {
      for (int i = 0; i < dimension; ++i)
        {
          ug_gdl_(face,i) /= cg(face);								// ug = sum(-nu.gradci)/cg
        }
    }
  //tab_divide_any_shape(ug_gdl_.valeurs(),cg_gdl_.valeurs());	// ug = sum(-nu.gradci)/cg
  ug_gdl_.mettre_a_jour(temps);

  // MAJ fraction molaire Xi = ci/cg
  DoubleTab& Xi = Xi_gdl_.valeurs();
  for (int face = 0; face < Xi.dimension(0); ++face)
    {
      for (int j = 0; j < nb_comp; ++j)
        {
          Xi(face,j) = ci(face,j)/cg(face);
        }
      // Debug Xi
      //if(face==10)
      //Cerr << "DEBUG Xi_GDL: face=" << face << ", cg=" << cg(face) << ", ci=" << ci(face,0) << " " << ci(face,1) << " " << ci(face,2) << ", Xi=" << Xi(face,0) << " " << Xi(face,1) << " " << Xi(face,2) << finl;
      // End Debug Xi
    }
  Xi_gdl_.mettre_a_jour(temps);

  // here we calculate the integral surfacique de flux de concentration Ns
  const Zone_VEF& zone_VEF = ref_cast(Zone_VEF,equation().zone_dis().valeur());
  const IntTab& face_voisins = zone_VEF.face_voisins();
  const DoubleTab& face_norm = zone_VEF.face_normales();

  // Debug ci impose sur le bord (GDL) au temps present -> pareil a ci sur le bord (CH) au temps precedent
  for (int i = 0; i < zone_VEF.nb_front_Cl(); ++i)
    {
      const Front_VF& bord = zone_VEF.front_VF(i);
      const Cond_lim& la_cl = equation().zone_Cl_dis().les_conditions_limites(i);
      if(bord.le_nom()=="I_CH_GDL")
        {
          //assert(sub_type(Scalaire_impose_paroi, la_cl.valeur()));
          if(sub_type(Scalaire_impose_paroi, la_cl.valeur()))
            {

              const Scalaire_impose_paroi la_cl_ci = ref_cast(Scalaire_impose_paroi, la_cl.valeur());
              // valeur ci imposee
              const DoubleTab& val= la_cl_ci.champ_front().valeur().valeurs();
              Cerr << "DEBUG ci impose sur le bord bord=" << bord.le_nom() << "au temps defaut=" << la_cl_ci.champ_front().valeur().get_temps_defaut() << finl;
              int nfaces = bord.nb_faces();
              for (int ind_face=0; ind_face < nfaces; ind_face++)
                {
                  int face = bord.num_face(ind_face);
                  Cerr << "DEBUG ci: , face=" << face << ", ind_face_bord=" << ind_face << ", surface=" << zone_VEF.face_surfaces(face) << ", ci_impose(face)=" << val(ind_face,0) << " " << val(ind_face,1) << " " << val(ind_face,2) << finl;
                }
            }

        }
    }
  // End ci impose sur le bord

  // Debug Ug sur le bord -> done ug = sum(Ni)/cg
  for (int i = 0; i < zone_VEF.nb_front_Cl(); ++i)
    {
      const Front_VF& bord = zone_VEF.front_VF(i);
      if(bord.le_nom()=="I_CH_GDL")
        {
          Cerr << "DEBUG ug_GDL sur le bord bord=" << bord.le_nom() << "au temps=" << temps << finl;
          int nfaces = bord.nb_faces();
          for (int ind_face=0; ind_face < nfaces; ind_face++)
            {
              int face = bord.num_face(ind_face);
              if(dimension == 2)
                Cerr << "DEBUG ug_GDL: face=" << face << ", ug=" << ug_gdl_(face,0) << " " << ug_gdl_(face,1) << finl;
              else
                Cerr << "DEBUG ug_GDL: face=" << face << ", ug=" << ug_gdl_(face,0) << " " << ug_gdl_(face,1) << " " << ug_gdl_(face,2) << finl;
            }
        }
    }
  // End Debug Ug

  // MAJ Ns flux surfacique de diffusion Ns = (integral(-nu.grad(ci)).n.dS) / S -> checked with flux_bords OK
  DoubleTab& Ns = Ns_gdl_.valeurs();
  for (int face = 0; face < Ns.dimension(0); ++face)
    {
      double face_surf = zone_VEF.face_surfaces(face);
      int n0 = face_voisins(face,0);
      int n1 = face_voisins(face,1);

      for (int i = 0; i < nb_comp; ++i)
        {
          double res = 0;
          for (int j = 0; j < dimension; ++j)
            {
              double Niface = 0;
              if(n0 != -1 && n1 != -1)
                {
                  // Ni moyen de deux elements
                  Niface = 0.5 * (Ni(n0,i*dimension+j)+Ni(n1,i*dimension+j));
                }
              else if(n0 != -1)
                {
                  // Ni de element nO
                  Niface = Ni(n0,i*dimension+j);
                }
              else
                {
                  // Ni de element n1
                  Niface = Ni(n1,i*dimension+j);
                }
              //res -= Niface*face_norm(face,j);				// Ns = -Ni*norma*Surface // - car normale sortante -> a verifier
              res += Niface*face_norm(face,j);				// Ns = +Ni*norma*Surface
            }
          Ns(face,i) = res / face_surf;						// Ns = -Ni*norma
        }
    }
  Ns_gdl_.mettre_a_jour(temps);

  // Debug Ns flux_bord -> done OK (compare with flux_bords post-traitement)
  int nb_front = zone_VEF.nb_front_Cl();
  Cerr << "DEBUG flux_bord dom=" << zone_VEF.domaine_dis().domaine().le_nom() << ", nombre de bords=" << nb_front << finl;
  for (int i = 0; i < nb_front; ++i)
    {
      const Front_VF& bord = zone_VEF.front_VF(i);
      Cerr << "DEBUG flux_bord au temps=" << temps << ", i=" << i << ", bord=" << bord.le_nom() << finl;
      if(bord.le_nom()=="I_CH_GDL")
        {
          int nfaces = bord.nb_faces();
          for (int ind_face=0; ind_face < nfaces; ind_face++)
            {
              int face = bord.num_face(ind_face);
              Cerr << "DEBUG flux_bord: , bord=" << bord.le_nom() << ", face=" << face << ", ind_face_bord=" << ind_face << ", surface=" << zone_VEF.face_surfaces(face) << ", Ns(face)=" << Ns(face,0) << " " << Ns(face,1) << " " << Ns(face,2) << " " << finl;
            }
        }

    }
  // End Debug Ns flux_bord
}

void Loi_Fermeture_PEMFC_base::mettre_a_jour(double temps)
{
  if(is_cas4_)
    Cerr << "Loi_Fermeture_PEMFC_base::mettre_a_jour Cas4 temps=" << temps << finl;
  else
    Cerr << "Loi_Fermeture_PEMFC_base::mettre_a_jour Cas2 temps=" << temps << finl;
  // update champ T
  if(ch_T_.non_nul())
    {
      ch_T_.valeur().mettre_a_jour(temps);
      const Zone_VF& zone_VF = ref_cast(Zone_VF,equation().zone_dis().valeur());
      const DoubleTab xp = zone_VF.xp();
      ch_T_.valeur().valeur_aux(xp, T_);
    }

  const Zone_VF& zone_VF = ref_cast(Zone_VF,equation().zone_dis().valeur());
  const DoubleTab& inconnue_org=equation().inconnue().valeurs();				// inco = concentration ci (cas2) fraction_molaire Xi (cas4)

  DoubleTab& nu=diffu_.valeurs();
  diffu_->changer_temps(temps);

  bidouille_nu(nu,inconnue_org,zone_VF);

  if (is_cas4_)
    calculer_Ni_ud(temps);
  else
    calculer_Ni_ug(temps); // TO-DO
}


void Loi_Fermeture_PEMFC_base::completer()
{
  Loi_Fermeture_base::completer();
  equation().zone_dis().zone().creer_tableau_elements(T_);
  if(nom_pb_T_ != "??")
    {
      Probleme_base& pb_T = ref_cast(Probleme_base,interprete().objet(nom_pb_T_));
      ch_T_ = pb_T.get_champ("temperature");
    }
  else
    {
      T_ = T_0;
    }
}



