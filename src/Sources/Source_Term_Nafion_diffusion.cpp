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
// File      : Source_Term_Nafion_diffusion.cpp
// Directory : $PEMFC_ROOT/src/Sources
//
/////////////////////////////////////////////////////////////////////////////

#include <Source_Term_Nafion_diffusion.h>
#include <Probleme_base.h>
#include <Conduction.h>
#include <Param.h>
#include <Interprete.h>
#include <Matrice_Morse.h>
#include <Champ_base.h>
#include <Champ_Inc.h>
#include <Domaine.h>
#include <Zone_VF.h>

Implemente_instanciable_sans_constructeur( Source_Term_Nafion_diffusion, "Source_Term_Nafion_diffusion_VEF_P1NC", Source_base ) ;

Source_Term_Nafion_diffusion::Source_Term_Nafion_diffusion(void)
{
  //T_0_ = 353.15;
  C_SO3_ = 2036.;
}

Sortie& Source_Term_Nafion_diffusion::printOn( Sortie& os ) const
{
  Source_base::printOn( os );
  return os;
}

Entree& Source_Term_Nafion_diffusion::readOn( Entree& is )
{
  //Source_base::readOn( is );
  Cerr << " Source_Term_Nafion_diffusion::readOn " << finl  ;
  Param param(que_suis_je());
  set_param(param);
  param.lire_avec_accolades(is);

  Motcles nom_especes_compris_(5);
  nom_especes_compris_[0] = "H2";
  nom_especes_compris_[1] = "O2";
  nom_especes_compris_[2] = "H2O";
  nom_especes_compris_[3] = "N2";
  nom_especes_compris_[4] = "vap";
  if(nom_especes_compris_.search(nom_espece_) == -1)
    {
      Cerr <<" unknown species in the list "<<finl;
      Process::exit();
    }
  //thickness_ionomer_ = (1-por_naf_)*eps_naf_ / gamma_CL_;
  dom_ = equation().probleme().domaine();
  dom_.valeur().creer_tableau_elements(C_);
  dom_.valeur().creer_tableau_elements(diffu_);
  dom_.valeur().creer_tableau_elements(ci_);
  dom_.valeur().creer_tableau_elements(T_);

  if(nom_espece_ == "H2")
    {
      // check
      assert(nom_ssz_CLa_ != "??");
      assert((nom_pb_ci_anode_ != "??" && nom_champ_ci_anode_ != "??") || concentration_.non_nul());
    }
  else if (nom_espece_ == "O2")
    {
      assert(nom_ssz_CLc_ != "??");
      assert((nom_pb_ci_cathode_ != "??" && nom_champ_ci_cathode_ != "??") || concentration_.non_nul());
    }
  else
    {
      assert(nom_ssz_CLa_ != "??" && nom_ssz_CLc_ != "??");
      assert((nom_pb_ci_cathode_ != "??" && nom_champ_ci_cathode_ != "??" && nom_pb_ci_anode_ != "??" && nom_champ_ci_anode_ != "??")|| concentration_.non_nul());
    }
  if(nom_ssz_CLa_ != "??")
    CL_a_ = dom_.valeur().ss_zone(nom_ssz_CLa_);
  if(nom_ssz_CLc_ != "??")
    CL_c_ = dom_.valeur().ss_zone(nom_ssz_CLc_);

  // discretiser le champ source -> pour post-traiter (diffusion des multi-especes dans le milieu poreux)
  discretiser(equation().discretisation());

  return is;
}

void Source_Term_Nafion_diffusion::associer_pb(const Probleme_base& pb)
{
  Cerr << " Source_Term_Nafion_diffusion::associer_pb " << finl ;
  assert(pb.que_suis_je() == "Pb_Conduction");
}

void Source_Term_Nafion_diffusion::associer_zones(const Zone_dis& zone_dis, const Zone_Cl_dis& zcl_dis)
{
  la_zone_ = ref_cast(Zone_VF,zone_dis.valeur());
}

void Source_Term_Nafion_diffusion::set_param(Param& param)
{
  param.ajouter("nom_espece",&nom_espece_,Param::REQUIRED); // XD_ADD_P chaine in list 'O2' 'H2' 'H2O' 'vap' 'N2'
  param.ajouter("nom_domaine",&nom_domaine_,Param::REQUIRED); // XD_ADD_P chaine domaine's name
  param.ajouter("nom_CLa",&nom_ssz_CLa_,Param::OPTIONAL); // XD_ADD_P chaine sub-area where the source exists, required for H2, N2, H2O
  param.ajouter("nom_CLc",&nom_ssz_CLc_,Param::OPTIONAL); // XD_ADD_P chaine sub-area where the source exists, required for O2, N2, H2O
  param.ajouter("por_naf", &por_naf_, Param::REQUIRED);	// XD_ADD_P Champ_Don porosity
  param.ajouter("eps_naf", &eps_naf_, Param::REQUIRED); // XD_ADD_P Champ_Don ionomer proportion
  param.ajouter("gamma_CL", &gamma_CL_, Param::REQUIRED); // XD_ADD_P Champ_Don specific surface of Platine atom m2/m3
  param.ajouter("C_SO3", &C_SO3_, Param::OPTIONAL);		// XD_ADD_P double default value of concentration mol/m3, required if H2O
  param.ajouter("nom_pb_C", &nom_pb_C_, Param::REQUIRED);	// XD_ADD_P chaine name of problem Co
  param.ajouter("nom_champ_C", &nom_champ_C_, Param::REQUIRED);	// XD_ADD_P chaine default 'temperature'
  param.ajouter("nom_champ_D", &nom_champ_D_, Param::REQUIRED); // XD_ADD_P chaine default 'diffusion_nafion'

  param.ajouter("nom_pb_T", &nom_pb_T_, Param::OPTIONAL);	// XD_ADD_P chaine name of problem of temperature
  param.ajouter("nom_champ_T", &nom_champ_T_, Param::OPTIONAL);	// XD_ADD_P chaine default 'temperature'
  param.ajouter("temperature", &temperature_, Param::OPTIONAL); // XD_ADD_P Champ_Don if the field of temperature is given

  param.ajouter("nom_pb_ci_cathode", &nom_pb_ci_cathode_, Param::OPTIONAL);// XD_ADD_P chaine name of problem Stefan Maxwell for cathode,required for O2, N2, H2O
  param.ajouter("nom_champ_ci_cathode", &nom_champ_ci_cathode_, Param::OPTIONAL); // XD_ADD_P chaine default 'concentration'
  param.ajouter("nom_pb_ci_anode", &nom_pb_ci_anode_, Param::OPTIONAL);// XD_ADD_P chaine name of problem Stefan Maxwell for anode, required for H2, N2, H2O
  param.ajouter("nom_champ_ci_anode", &nom_champ_ci_anode_, Param::OPTIONAL);// XD_ADD_P chaine default 'concentration'
  param.ajouter("concentration", &concentration_, Param::OPTIONAL);// XD_ADD_P chaine given a field of 'concentration'
}

void Source_Term_Nafion_diffusion::discretiser(const Discretisation_base& dis)
{
  dis.discretiser_champ("champ_elem",equation().zone_dis().valeur(),"S","unit", 1 , 0. , ch_S_);
  champs_compris_.ajoute_champ(ch_S_);

  //dis.discretiser_champ("temperature",equation().zone_dis().valeur(),"Ceq","unit", 1 , 0. , ch_Ceq_);
  dis.discretiser_champ("champ_elem",equation().zone_dis().valeur(),"Ceq","unit", 1 , 0. , ch_Ceq_);
  champs_compris_.ajoute_champ(ch_Ceq_);
}

// Sa = D.gamma/ep_naf*(Ceq - C)
double Source_Term_Nafion_diffusion::f_Sa(double diffu, double Ci, double ci, double T, double por, double eps, double gamma) const
{
  // expression_source
  double Ceq = 0.;
  double H = 0.;
  double Psat = 0.;
  double a_H2O = 0.;
  //double a_H2O_lim = 0.;
  double lambda_eq = 0.;
  if (nom_espece_ == "H2O" || nom_espece_ == "vap")
    {
      Psat = f_Psat(T);
      a_H2O = ci * R * T / Psat;
      //a_H2O_lim = max(a_lim, a_H2O);
      //lambda_eq = f_lambda(a_H2O_lim);
      lambda_eq = f_lambda(a_H2O);
      Ceq = lambda_eq * C_SO3_;
    }
  else
    {
      if (nom_espece_ == "H2")
        {
          H = f_Henry_H2(T);
        }
      else if (nom_espece_ == "O2")
        {
          H = f_Henry_O2(T);
        }
      else if (nom_espece_ == "N2")
        {
          H = f_Henry_N2(T);
        }
      else
        {
          Cerr <<" unknown species in the list "<<finl;
          Process::exit();
        }
      Ceq = ci * R * T * H;
    }

  double e_naf = (1-por)*eps / gamma;
  double val = diffu * gamma / e_naf * (Ceq - Ci);
  return val;
}

double Source_Term_Nafion_diffusion::f_Ceq(double Ci, double ci, double T) const
{
  double Ceq = 0.;
  double H = 0.;
  double Psat = 0.;
  double a_H2O = 0.;
  //double a_H2O_lim = 0.;
  double lambda_eq = 0.;
  if (nom_espece_ == "H2O" || nom_espece_ == "vap")
    {
      Psat = f_Psat(T);
      a_H2O = ci * R * T / Psat;
      //a_H2O_lim = max(a_lim, a_H2O);
      //lambda_eq = f_lambda(a_H2O_lim);
      lambda_eq = f_lambda(a_H2O);
      Ceq = lambda_eq * C_SO3_;
    }
  else
    {
      if (nom_espece_ == "H2")
        {
          H = f_Henry_H2(T);
        }
      else if (nom_espece_ == "O2")
        {
          H = f_Henry_O2(T);
        }
      else if (nom_espece_ == "N2")
        {
          H = f_Henry_N2(T);
        }
      else
        {
          Cerr <<" unknown species in the list "<<finl;
          Process::exit();
        }
      Ceq = ci * R * T * H;
    }

  return Ceq;
}

// derivee de Sa par rapport C = dSa/dc = -Dgamma/ep_naf
double Source_Term_Nafion_diffusion::f_dSadC(double diffu, double por, double eps, double gamma) const
{
  // expression_derivee_par_rapport_inconnue
  double e_naf = (1-por)*eps/gamma;
  return - diffu * gamma / e_naf;
}

DoubleTab& Source_Term_Nafion_diffusion::calculer(DoubleTab& resu) const
{
  resu = 0;
  ajouter(resu);
  return resu;
}

void Source_Term_Nafion_diffusion::mettre_a_jour(double temps)
{
  Cerr << "Source_Term_Nafion_diffusion::mettre_a_jour " << equation().probleme().le_nom() << finl;

  // update the involving fields
  ch_C_.valeur().mettre_a_jour(temps);
  ch_D_i_naf_.valeur().mettre_a_jour(temps);
  if(ch_ci_cathode_.non_nul())
    ch_ci_cathode_.valeur().mettre_a_jour(temps);
  if(ch_ci_anode_.non_nul())
    ch_ci_anode_.valeur().mettre_a_jour(temps);
  if(ch_T_.non_nul())
    ch_T_.valeur().mettre_a_jour(temps);
  if(temperature_.non_nul())
    temperature_.mettre_a_jour(temps);
  if(concentration_.non_nul())
    concentration_.mettre_a_jour(temps);

  const DoubleTab& xp=la_zone_.valeur().xp();
  ch_C_.valeur().valeur_aux(xp, C_);
  C_.echange_espace_virtuel();
  ch_D_i_naf_.valeur().valeur_aux(xp, diffu_);
  diffu_.echange_espace_virtuel();

  //assert(xp.dimension_tot(0) == ci_.dimension_tot(0));
  //Cerr << "dimension tot xp " << xp.dimension_tot(0) << finl;

  if(ch_T_.non_nul())
    {
      ch_T_.valeur().valeur_aux(xp, T_);
      T_.echange_espace_virtuel();
    }

  if(ch_ci_cathode_.non_nul() && !ch_ci_anode_.non_nul())
    {
      // case O2
      assert(nom_espece_ == "O2");
      ch_ci_cathode_.valeur().valeur_aux_compo(xp, ci_, 0);										// ncomp = 0 pour O2
      ci_.echange_espace_virtuel();
    }
  else if(ch_ci_anode_.non_nul() && !ch_ci_cathode_.non_nul())
    {
      // case H2
      assert(nom_espece_=="H2");
      ch_ci_anode_.valeur().valeur_aux_compo(xp, ci_, 0);										// ncomp = 0 pour O2
      ci_.echange_espace_virtuel();
    }
  else if(ch_ci_anode_.non_nul() && ch_ci_cathode_.non_nul())
    {
      // case N2, H2O
      assert(nom_espece_=="N2"||nom_espece_=="H2O"||nom_espece_=="vap");

      DoubleTab val_ci_cathode, val_ci_anode;
      la_zone_.valeur().zone().creer_tableau_elements(val_ci_cathode);
      la_zone_.valeur().zone().creer_tableau_elements(val_ci_anode);
      if(nom_espece_ == "N2")
        {
          ch_ci_cathode_.valeur().valeur_aux_compo(xp, val_ci_cathode, 2);							// ncomp = 1 pour H20, 2 pour N2
          ch_ci_anode_.valeur().valeur_aux_compo(xp, val_ci_anode, 2);
        }
      else
        {
          ch_ci_cathode_.valeur().valeur_aux_compo(xp, val_ci_cathode, 1);							// ncomp = 1 pour H20, 2 pour N2
          ch_ci_anode_.valeur().valeur_aux_compo(xp, val_ci_anode, 1);
        }
      val_ci_cathode.echange_espace_virtuel();
      val_ci_anode.echange_espace_virtuel();

      for (int poly = 0; poly < CL_a_.valeur().nb_elem_tot(); poly++)
        {
          int elem = CL_a_.valeur()(poly);
          ci_(elem) = val_ci_anode(elem);
        }
      for (int poly = 0; poly < CL_c_.valeur().nb_elem_tot(); poly++)
        {
          int elem = CL_c_.valeur()(poly);
          ci_(elem) = val_ci_cathode(elem);
        }
    }
  if(temperature_.non_nul())
    {
      temperature_.valeur().valeur_aux(xp, T_);
      T_.echange_espace_virtuel();
    }
  if(concentration_.non_nul())
    {
      concentration_.valeur().valeur_aux(xp, ci_);
      ci_.echange_espace_virtuel();
    }

  // mettre a jour ch_S (P0)
  DoubleTab& val_S = ch_S_.valeurs();
  DoubleTab& por = por_naf_.valeurs();
  DoubleTab& eps = eps_naf_.valeurs();
  DoubleTab& gamma = gamma_CL_.valeurs();

  if(CL_a_.non_nul())
    {
      for (int poly = 0; poly < CL_a_.valeur().nb_elem_tot(); poly++)
        {
          int elem = CL_a_.valeur()(poly);
          val_S(elem) = f_Sa(diffu_(elem), C_(elem), ci_(elem), T_(elem), por(elem,0), eps(elem, 0), gamma(elem,0));
        }
    }

  if(CL_c_.non_nul())
    {
      for (int poly = 0; poly < CL_c_.valeur().nb_elem_tot(); poly++)
        {
          int elem = CL_c_.valeur()(poly);
          val_S(elem) = f_Sa(diffu_(elem), C_(elem), ci_(elem), T_(elem), por(elem,0), eps(elem, 0), gamma(elem,0));

        }
    }
  ch_S_.mettre_a_jour(temps);

  /*
    // mettre a jour ch_Ceq (P1NC)
    DoubleTab& val_Ceq = ch_Ceq_.valeurs();
    const DoubleTab& xv=la_zone_.valeur().xv();		// bary centre des faces
    DoubleTab Tf;
    la_zone_.valeur().creer_tableau_faces(Tf);
    ch_T_.valeur().valeur_aux(xv, Tf);
    Tf.echange_espace_virtuel();

    DoubleTab& Cf = ch_C_.valeur().valeurs();
    DoubleTab cif;
    la_zone_.valeur().creer_tableau_faces(cif);
    if(ch_ci_cathode_.non_nul() && !ch_ci_anode_.non_nul())
      {
        // case O2
        ch_ci_cathode_.valeur().valeur_aux_compo(xv, cif, 0);										// ncomp = 0 pour O2
      }
    else if(ch_ci_anode_.non_nul() && !ch_ci_cathode_.non_nul())
      {
        // case H2
        ch_ci_anode_.valeur().valeur_aux_compo(xv, cif, 0);										// ncomp = 0 pour H2
      }
    else
      {
        // case N2, H20
        DoubleTab val_ci_cathode, val_ci_anode;
        la_zone_.valeur().creer_tableau_faces(val_ci_cathode);
        la_zone_.valeur().creer_tableau_faces(val_ci_anode);
        if(nom_espece_ == "N2")
          {
            ch_ci_cathode_.valeur().valeur_aux_compo(xv, val_ci_cathode, 2);							// ncomp = 1 pour H20, 2 pour N2
            ch_ci_anode_.valeur().valeur_aux_compo(xv, val_ci_anode, 2);
          }
        else
          {
            ch_ci_cathode_.valeur().valeur_aux_compo(xv, val_ci_cathode, 1);							// ncomp = 1 pour H20, 2 pour N2
            ch_ci_anode_.valeur().valeur_aux_compo(xv, val_ci_anode, 1);
          }
        val_ci_anode.echange_espace_virtuel();
        val_ci_cathode.echange_espace_virtuel();

        for (int face = 0; face < la_zone_.valeur().nb_faces(); ++face)
          {
            cif(face) = val_ci_anode(face) + val_ci_cathode(face);
          }
        //cif = val_ci_anode;
        //cif += val_ci_cathode;
      }


    for (int face = 0; face < la_zone_.valeur().nb_faces(); ++face)
      {
        val_Ceq(face) = f_Ceq(Cf(face), cif(face), Tf(face));						// Verifier
      }
   */

  // mettre a jour ch_Ceq (P0)
  DoubleTab& val_Ceq = ch_Ceq_.valeurs();
  if(CL_a_.non_nul())
    {
      for (int poly = 0; poly < CL_a_.valeur().nb_elem_tot(); poly++)
        {
          int elem = CL_a_.valeur()(poly);
          val_Ceq(elem) = f_Ceq(C_(elem), ci_(elem), T_(elem));
        }
    }
  if(CL_c_.non_nul())
    {
      for (int poly = 0; poly < CL_c_.valeur().nb_elem_tot(); poly++)
        {
          int elem = CL_c_.valeur()(poly);
          val_Ceq(elem) = f_Ceq(C_(elem), ci_(elem), T_(elem));
        }
    }
  ch_Ceq_.mettre_a_jour(temps);

  Cerr << "concentration min max " << mp_min_vect(ci_) << " " << mp_max_vect(ci_) << finl;
  Cerr << "diffu min max " << mp_min_vect(diffu_) << " " << mp_max_vect(diffu_) << finl;
  Cerr << "source Sa min max " << mp_min_vect(val_S) << " " << mp_max_vect(val_S) << finl;
  Cerr << "Ceq min max " << mp_min_vect(val_Ceq) << " " << mp_max_vect(val_Ceq) << finl;
}

void Source_Term_Nafion_diffusion::completer()
{
  Source_base::completer();
  // get the reference to the coupling fields
  Probleme_base& pb_C = ref_cast(Probleme_base,interprete().objet(nom_pb_C_));
  ch_C_ = pb_C.get_champ(nom_champ_C_);
  ch_D_i_naf_ = pb_C.get_champ(nom_champ_D_);

  if(nom_pb_T_ != "??")
    {
      Probleme_base& pb_T = ref_cast(Probleme_base,interprete().objet(nom_pb_T_));
      ch_T_ = pb_T.get_champ(nom_champ_T_);
    }

  if(nom_pb_ci_cathode_ != "??")
    {
      Probleme_base& pb_ci_cathode = ref_cast(Probleme_base,interprete().objet(nom_pb_ci_cathode_));
      ch_ci_cathode_ = pb_ci_cathode.get_champ(nom_champ_ci_cathode_);
    }
  if(nom_pb_ci_anode_ != "??")
    {
      Probleme_base& pb_ci_anode = ref_cast(Probleme_base,interprete().objet(nom_pb_ci_anode_));
      ch_ci_anode_ = pb_ci_anode.get_champ(nom_champ_ci_anode_);
    }
}


DoubleTab& Source_Term_Nafion_diffusion::ajouter(DoubleTab& resu) const
{

  assert(resu.dimension(0)==la_zone_.valeur().nb_faces());
  const DoubleTab& por = por_naf_.valeurs();
  const DoubleTab& eps = eps_naf_.valeurs();
  const DoubleTab& gamma = gamma_CL_.valeurs();

  DoubleVect vol = la_zone_.valeur().volumes();

  if(CL_a_.non_nul())
    {
      for (int poly = 0; poly < CL_a_.valeur().nb_elem_tot(); poly++)
        {
          int elem = CL_a_.valeur()(poly);
          double Sa = f_Sa(diffu_(elem), C_(elem), ci_(elem), T_(elem), por(elem,0), eps(elem, 0), gamma(elem,0));
          double inv_rhoCp = 1./((1-por(elem,0))*eps(elem,0));

          int nb_face_elem = la_zone_.valeur().zone().nb_faces_elem(0);
          for (int f = 0; f < nb_face_elem; ++f)
            {
              int face = la_zone_.valeur().elem_faces(elem, f);
              resu(face) +=  Sa * inv_rhoCp * vol(elem) / nb_face_elem ;
            }
        }
    }
  if(CL_c_.non_nul())
    {
      {
        for (int poly = 0; poly < CL_c_.valeur().nb_elem_tot(); poly++)
          {
            int elem = CL_c_.valeur()(poly);
            double Sa = f_Sa(diffu_(elem), C_(elem), ci_(elem), T_(elem), por(elem,0), eps(elem, 0), gamma(elem,0));
            double inv_rhoCp = 1./((1-por(elem,0))*eps(elem,0));
            int nb_face_elem = la_zone_.valeur().zone().nb_faces_elem(0);
            for (int f = 0; f < nb_face_elem; ++f)
              {
                int face = la_zone_.valeur().elem_faces(elem, f);
                resu(face) +=  Sa * inv_rhoCp * vol(elem) / nb_face_elem ;
              }
          }
      }
    }

  return resu;
}

void Source_Term_Nafion_diffusion::contribuer_a_avec(const DoubleTab& inco, Matrice_Morse& mat) const
{

  DoubleVect vol = la_zone_.valeur().volumes();
  const DoubleTab& por = por_naf_.valeurs();
  const DoubleTab& eps = eps_naf_.valeurs();
  const DoubleTab& gamma = gamma_CL_.valeurs();

  if(CL_a_.non_nul())
    {
      for (int poly = 0; poly < CL_a_.valeur().nb_elem_tot(); poly++)
        {
          int elem = CL_a_.valeur()(poly);
          double terme = - f_dSadC(diffu_(elem), por(elem,0), eps(elem,0), gamma(elem,0));	// quang: 23/05/19 correction: -dSdC
          double coef = 1./((1-por(elem,0))*eps(elem,0));
          int nb_face_elem = la_zone_.valeur().zone().nb_faces_elem(0);

          for (int f = 0; f < nb_face_elem; ++f)
            {
              int face = la_zone_.valeur().elem_faces(elem, f);
              mat.coef(face,face) +=  terme * coef * vol(elem) / nb_face_elem ;
            }
        }
    }
  if(CL_c_.non_nul())
    {
      {
        for (int poly = 0; poly < CL_c_.valeur().nb_elem_tot(); poly++)
          {
            int elem = CL_c_.valeur()(poly);
            double terme = - f_dSadC(diffu_(elem), por(elem,0), eps(elem,0), gamma(elem,0));		// quang: 23/05/19 correction: -dSdC
            double coef = 1./((1-por(elem,0))*eps(elem,0));
            int nb_face_elem = la_zone_.valeur().zone().nb_faces_elem(0);

            for (int f = 0; f < nb_face_elem; ++f)
              {
                int face = la_zone_.valeur().elem_faces(elem, f);
                mat.coef(face,face) +=  terme * coef * vol(elem) / nb_face_elem ;
              }
          }
      }
    }
}
