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
// File      : Young_Todd_Source_multi.cpp
// Directory : $PEMFC_ROOT/src/Sources
//
/////////////////////////////////////////////////////////////////////////////

#include <Young_Todd_Source_multi.h>
#include <Probleme_base.h>
#include <Conduction.h>
#include <Param.h>
#include <Interprete.h>
#include <Matrice_Morse.h>
#include <Champ_base.h>
#include <Domaine.h>
#include <Zone_VF.h>

Implemente_instanciable( Young_Todd_Source_multi, "Young_Todd_Source_multi_VEF_P1NC", Source_base ) ;

Sortie& Young_Todd_Source_multi::printOn( Sortie& os ) const
{
  Source_base::printOn( os );
  return os;
}

Entree& Young_Todd_Source_multi::readOn( Entree& is )
{
  Source_base::readOn( is );
  Cerr << " Young_Todd_Source_multi::readOn " << finl  ;
  Param param(que_suis_je());
  set_param(param);
  param.lire_avec_accolades(is);

  dom_ = ref_cast(Domaine, Interprete::objet(nom_dom_));
  if(nom_ssz_CLa_ != "??")
    {
      assert(nom_ssz_CLc_ == "??");
      CL_a_ = dom_.valeur().ss_zone(nom_ssz_CLa_);
    }
  if(nom_ssz_CLc_ != "??")
    {
      assert(nom_ssz_CLa_ == "??");
      CL_c_ = dom_.valeur().ss_zone(nom_ssz_CLc_);
    }
  return is;
}

// ajouter S_t
DoubleTab& Young_Todd_Source_multi::ajouter(DoubleTab& resu) const
{
  assert(resu.dimension(0)==la_zone_.valeur().nb_faces());
  int nb_comp = resu.dimension(1);
  DoubleVect vol = la_zone_.valeur().volumes();

  if(CL_a_.non_nul())
    {
      for (int poly = 0; poly < CL_a_.valeur().nb_elem_tot(); poly++)
        {
          int elem = CL_a_.valeur()(poly);
          int nb_face_elem = la_zone_.valeur().zone().nb_faces_elem(0);

          for (int f = 0; f < nb_face_elem; ++f)
            {
              int face = la_zone_.valeur().elem_faces(elem, f);
              for (int j = 0; j < nb_comp; ++j)
                {
                  resu(face,j) +=  S_(elem,j) * vol(elem) / nb_face_elem ;
                }
            }
        }
    }
  if(CL_c_.non_nul())
    {
      {
        for (int poly = 0; poly < CL_c_.valeur().nb_elem_tot(); poly++)
          {
            int elem = CL_c_.valeur()(poly);
            int nb_face_elem = la_zone_.valeur().zone().nb_faces_elem(0);
            for (int f = 0; f < nb_face_elem; ++f)
              {
                int face = la_zone_.valeur().elem_faces(elem, f);
                for (int j = 0; j < nb_comp; ++j)
                  {
                    resu(face,j) +=  S_(elem,j) * vol(elem) / nb_face_elem ;
                  }
              }
          }
      }
    }

  return resu;
}


DoubleTab& Young_Todd_Source_multi::calculer(DoubleTab& resu) const
{
  resu = 0;
  ajouter(resu);
  return resu;
}

void Young_Todd_Source_multi::mettre_a_jour(double temps)
{
  Cerr << "Young_Todd_Source_multi::mettre_a_jour " << equation().probleme().le_nom() << finl;

  // update the involving fields -> to-do: nothing
  ch_ci_.valeur().mettre_a_jour(temps);
  ch_T_.valeur().mettre_a_jour(temps);
  ch_C_X2_.valeur().mettre_a_jour(temps);
  ch_D_X2_.valeur().mettre_a_jour(temps);
  ch_C_vap_.valeur().mettre_a_jour(temps);
  ch_D_vap_.valeur().mettre_a_jour(temps);
  ch_C_N2_.valeur().mettre_a_jour(temps);
  ch_D_N2_.valeur().mettre_a_jour(temps);

  const DoubleTab& xp=la_zone_.valeur().xp();		// centre gravite des elements
  ch_C_X2_.valeur().valeur_aux(xp, C_X2);
  ch_D_X2_.valeur().valeur_aux(xp, D_X2);
  ch_C_vap_.valeur().valeur_aux(xp, C_vap);
  ch_D_vap_.valeur().valeur_aux(xp, D_vap);
  ch_C_N2_.valeur().valeur_aux(xp, C_N2);
  ch_D_N2_.valeur().valeur_aux(xp, D_N2);
  ch_T_.valeur().valeur_aux(xp, T_);
  ch_ci_.valeur().valeur_aux(xp, ci_);
  C_X2.echange_espace_virtuel();
  D_X2.echange_espace_virtuel();
  C_vap.echange_espace_virtuel();
  D_vap.echange_espace_virtuel();
  C_N2.echange_espace_virtuel();
  D_N2.echange_espace_virtuel();
  T_.echange_espace_virtuel();
  ci_.echange_espace_virtuel();

  // mettre a jour les tableaux de S et dSdc
  const DoubleTab& por = por_naf_.valeurs();
  const DoubleTab& eps = eps_naf_.valeurs();
  const DoubleTab& gamma = gamma_CL_.valeurs();

  if(CL_a_.non_nul())
    {
      for (int poly = 0; poly < CL_a_.valeur().nb_elem_tot(); poly++)
        {
          int elem = CL_a_.valeur()(poly);
          double S_X2 = f_S_H2(D_X2(elem), ci_(elem,0), C_X2(elem), T_(elem), por(elem,0), eps(elem, 0), gamma(elem,0));
          double dSdc_X2 = f_dSdc_H2(D_X2(elem), ci_(elem,0), C_X2(elem), T_(elem), por(elem,0), eps(elem, 0), gamma(elem,0));

          double S_vap = f_S_vap(D_vap(elem), ci_(elem,1), C_vap(elem), T_(elem), por(elem,0), eps(elem, 0), gamma(elem,0));
          double dSdc_vap = f_dSdc_vap(D_vap(elem), ci_(elem,1), C_vap(elem), T_(elem), por(elem,0), eps(elem, 0), gamma(elem,0));

          double S_N2 = f_S_N2(D_N2(elem), ci_(elem,2), C_N2(elem), T_(elem), por(elem,0), eps(elem, 0), gamma(elem,0));
          double dSdc_N2 = f_dSdc_N2(D_N2(elem), ci_(elem,2), C_N2(elem), T_(elem), por(elem,0), eps(elem, 0), gamma(elem,0));

          S_(elem,0) = S_X2;
          S_(elem,1) = S_vap;
          S_(elem,2) = S_N2;
          dSdc_(elem,0) = dSdc_X2;
          dSdc_(elem,1) = dSdc_vap;
          dSdc_(elem,2) = dSdc_N2;
        }
    }
  if(CL_c_.non_nul())
    {
      {
        for (int poly = 0; poly < CL_c_.valeur().nb_elem_tot(); poly++)
          {
            int elem = CL_c_.valeur()(poly);
            double S_X2 = f_S_O2(D_X2(elem), ci_(elem,0), C_X2(elem), T_(elem), por(elem,0), eps(elem, 0), gamma(elem,0));
            double dSdc_X2 = f_dSdc_O2(D_X2(elem), ci_(elem,0), C_X2(elem), T_(elem), por(elem,0), eps(elem, 0), gamma(elem,0));

            double S_vap = f_S_vap(D_vap(elem), ci_(elem,1), C_vap(elem), T_(elem), por(elem,0), eps(elem, 0), gamma(elem,0));
            double dSdc_vap = f_dSdc_vap(D_vap(elem), ci_(elem,1), C_vap(elem), T_(elem), por(elem,0), eps(elem, 0), gamma(elem,0));

            double S_N2 = f_S_N2(D_N2(elem), ci_(elem,2), C_N2(elem), T_(elem), por(elem,0), eps(elem, 0), gamma(elem,0));
            double dSdc_N2 = f_dSdc_N2(D_N2(elem), ci_(elem,2), C_N2(elem), T_(elem), por(elem,0), eps(elem, 0), gamma(elem,0));

            S_(elem,0) = S_X2;
            S_(elem,1) = S_vap;
            S_(elem,2) = S_N2;
            dSdc_(elem,0) = dSdc_X2;
            dSdc_(elem,1) = dSdc_vap;
            dSdc_(elem,2) = dSdc_N2;
          }
      }
    }
  Cerr << "C_X2 min max " << mp_min_vect(C_X2) << " " << mp_max_vect(C_X2) << finl;
  Cerr << "D_X2 min max " << mp_min_vect(D_X2) << " " << mp_max_vect(D_X2) << finl;
  Cerr << "C_vap min max " << mp_min_vect(C_vap) << " " << mp_max_vect(C_vap) << finl;
  Cerr << "D_vap min max " << mp_min_vect(D_vap) << " " << mp_max_vect(D_vap) << finl;
  Cerr << "C_N2 min max " << mp_min_vect(C_N2) << " " << mp_max_vect(C_N2) << finl;
  Cerr << "D_N2 min max " << mp_min_vect(D_N2) << " " << mp_max_vect(D_N2) << finl;
  Cerr << "T min max " << mp_min_vect(T_) << " " << mp_max_vect(T_) << finl;
  Cerr << "ci min max " << mp_min_vect(ci_) << " " << mp_max_vect(ci_) << finl;
  Cerr << "S min max " << mp_min_vect(S_) << " " << mp_max_vect(S_) << finl;
  Cerr << "dSdc min max " << mp_min_vect(dSdc_) << " " << mp_max_vect(dSdc_) << finl;
}

// ajouter -dS/dc
void Young_Todd_Source_multi::contribuer_a_avec(const DoubleTab& inco, Matrice_Morse& mat) const
{

  DoubleVect vol = la_zone_.valeur().volumes();

  int nb_comp = inco.dimension(1);

  if(CL_a_.non_nul())
    {
      for (int poly = 0; poly < CL_a_.valeur().nb_elem_tot(); poly++)
        {
          int elem = CL_a_.valeur()(poly);
          int nb_face_elem = la_zone_.valeur().zone().nb_faces_elem(0);


          for (int f = 0; f < nb_face_elem; ++f)
            {
              int face = la_zone_.valeur().elem_faces(elem, f);
              for (int j = 0; j < nb_comp; ++j)
                {
                  int idface = face * nb_comp + j;
                  mat.coef(idface,idface) +=  - dSdc_(elem,j) * vol(elem) / nb_face_elem ;
                }

            }
        }
    }
  if(CL_c_.non_nul())
    {
      {
        for (int poly = 0; poly < CL_c_.valeur().nb_elem_tot(); poly++)
          {
            int elem = CL_c_.valeur()(poly);
            int nb_face_elem = la_zone_.valeur().zone().nb_faces_elem(0);

            for (int f = 0; f < nb_face_elem; ++f)
              {
                int face = la_zone_.valeur().elem_faces(elem, f);
                for (int j = 0; j < nb_comp; ++j)
                  {
                    int idface = face * nb_comp + j;
                    mat.coef(idface,idface) +=  - dSdc_(elem,j) * vol(elem) / nb_face_elem ;
                  }

              }
          }
      }
    }
}

void Young_Todd_Source_multi::associer_zones(const Zone_dis& zone_dis, const Zone_Cl_dis& zcl_dis)
{
  la_zone_ = ref_cast(Zone_VF,zone_dis.valeur());
}

void Young_Todd_Source_multi::associer_pb(const Probleme_base& pb)
{
  // nothing to be associated
}

void Young_Todd_Source_multi::completer()
{
  Source_base::completer();
  // get the reference to the coupling fields

  Probleme_base& pb_T = ref_cast(Probleme_base,interprete().objet(nom_pb_T_));
  ch_T_ = pb_T.get_champ(nom_champ_T_);

  Probleme_base& pb_ci = ref_cast(Probleme_base,interprete().objet(nom_pb_ci_));
  ch_ci_ = pb_ci.get_champ(nom_champ_ci_);

  Probleme_base& pb_C_X2 = ref_cast(Probleme_base,interprete().objet(nom_pb_C_X2_));
  ch_C_X2_ = pb_C_X2.get_champ(nom_champ_C_X2_);
  ch_D_X2_ = pb_C_X2.get_champ(nom_champ_D_X2_);

  Probleme_base& pb_C_vap = ref_cast(Probleme_base,interprete().objet(nom_pb_C_vap_));
  ch_C_vap_ = pb_C_vap.get_champ(nom_champ_C_vap_);
  ch_D_vap_ = pb_C_vap.get_champ(nom_champ_D_vap_);

  Probleme_base& pb_C_N2 = ref_cast(Probleme_base,interprete().objet(nom_pb_C_N2_));
  ch_C_N2_ = pb_C_N2.get_champ(nom_champ_C_N2_);
  ch_D_N2_ = pb_C_N2.get_champ(nom_champ_D_N2_);

  int nb_comp = ch_ci_.valeur().valeurs().dimension(1);
  assert(nb_comp == 3);
  ci_.resize(0, nb_comp);
  S_.resize(0, nb_comp);
  dSdc_.resize(0, nb_comp);
  la_zone_.valeur().zone().creer_tableau_elements(C_X2);
  la_zone_.valeur().zone().creer_tableau_elements(C_vap);
  la_zone_.valeur().zone().creer_tableau_elements(C_N2);
  la_zone_.valeur().zone().creer_tableau_elements(D_X2);
  la_zone_.valeur().zone().creer_tableau_elements(D_vap);
  la_zone_.valeur().zone().creer_tableau_elements(D_N2);
  la_zone_.valeur().zone().creer_tableau_elements(ci_);
  la_zone_.valeur().zone().creer_tableau_elements(T_);
  la_zone_.valeur().zone().creer_tableau_elements(S_);
  la_zone_.valeur().zone().creer_tableau_elements(dSdc_);
}

void Young_Todd_Source_multi::set_param(Param& param)
{
  param.ajouter("nom_domaine",&nom_dom_,Param::REQUIRED);
  param.ajouter("nom_ssz_CLa",&nom_ssz_CLa_,Param::OPTIONAL);	// requis pour anode
  param.ajouter("nom_ssz_CLc",&nom_ssz_CLc_,Param::OPTIONAL);	// requis pour cathode
  param.ajouter("por_naf", &por_naf_, Param::REQUIRED);
  param.ajouter("eps_naf", &eps_naf_, Param::REQUIRED);
  param.ajouter("gamma_CL", &gamma_CL_, Param::REQUIRED);
  param.ajouter("nom_pb_C_ci", &nom_pb_ci_, Param::REQUIRED);
  param.ajouter("nom_champ_ci", &nom_champ_ci_, Param::REQUIRED);
  param.ajouter("nom_pb_C_X2", &nom_pb_C_X2_, Param::REQUIRED);
  param.ajouter("nom_champ_C_X2", &nom_champ_C_X2_, Param::REQUIRED);
  param.ajouter("nom_champ_D_X2", &nom_champ_D_X2_, Param::REQUIRED);
  param.ajouter("nom_pb_C_vap", &nom_pb_C_vap_, Param::REQUIRED);
  param.ajouter("nom_champ_C_vap", &nom_champ_C_vap_, Param::REQUIRED);
  param.ajouter("nom_champ_D_vap", &nom_champ_D_vap_, Param::REQUIRED);
  param.ajouter("nom_pb_C_N2", &nom_pb_C_N2_, Param::REQUIRED);
  param.ajouter("nom_champ_C_N2", &nom_champ_C_N2_, Param::REQUIRED);
  param.ajouter("nom_champ_D_N2", &nom_champ_D_N2_, Param::REQUIRED);
  param.ajouter("nom_pb_T", &nom_pb_T_, Param::REQUIRED);
  param.ajouter("nom_champ_T", &nom_champ_T_, Param::REQUIRED);
}

// S_X2 = -D.gamma_CL/ep_naf(cRTH - C)
double Young_Todd_Source_multi::f_S_H2(double diffu, double ci, double C, double T, double por, double eps, double gamma) const
{
  double H = f_Henry_H2(T);
  double Ceq = ci * R * T * H;
  double e_naf = (1-por)*eps / gamma;
  double S = - diffu * gamma / e_naf * (Ceq - C);
  return S;
}

// dSdc_X2 = -D.gamma_CL/ep_naf.RTH
double Young_Todd_Source_multi::f_dSdc_H2(double diffu, double ci, double C, double T, double por, double eps, double gamma) const
{
  double H = f_Henry_H2(T);
  double RTH = R * T * H;
  double e_naf = (1-por)*eps / gamma;
  double dSdc = - diffu * gamma / e_naf * RTH;
  return dSdc;
}

// S_O2 = -D.gamma_CL/ep_naf(cRTH - C)
double Young_Todd_Source_multi::f_S_O2(double diffu, double ci, double C, double T, double por, double eps, double gamma) const
{
  double H = f_Henry_O2(T);
  double Ceq = ci * R * T * H;
  double e_naf = (1-por)*eps / gamma;
  double S = - diffu * gamma / e_naf * (Ceq - C);
  return S;
}

// dSdc_O2 = -D.gamma_CL/ep_naf.RTH
double Young_Todd_Source_multi::f_dSdc_O2(double diffu, double ci, double C, double T, double por, double eps, double gamma) const
{
  double H = f_Henry_O2(T);
  double RTH = R * T * H;
  double e_naf = (1-por)*eps / gamma;
  double dSdc = - diffu * gamma / e_naf * RTH;
  return dSdc;
}

// S_N2 = -D.gamma_CL/ep_naf(cRTH - C)
double Young_Todd_Source_multi::f_S_N2(double diffu, double ci, double C, double T, double por, double eps, double gamma) const
{
  double H = f_Henry_N2(T);
  double Ceq = ci * R * T * H;
  double e_naf = (1-por)*eps / gamma;
  double S = - diffu * gamma / e_naf * (Ceq - C);
  return S;
}

// dSdc_N2 = -D.gamma_CL/ep_naf.RTH
double Young_Todd_Source_multi::f_dSdc_N2(double diffu, double ci, double C, double T, double por, double eps, double gamma) const
{
  double H = f_Henry_N2(T);
  double RTH = R * T * H;
  double e_naf = (1-por)*eps / gamma;
  double dSdc = - diffu * gamma / e_naf * RTH;
  return dSdc;
}

// S_vap = -Dw.gamma_CL/ep_naf(Ceq - C) avec Ceq = CSO3.f_lambda(a) et a=cRT/Psat
double Young_Todd_Source_multi::f_S_vap(double diffu, double ci, double C, double T, double por, double eps, double gamma) const
{
  double Psat = f_Psat(T);
  double a_H2O = ci * R * T / Psat;
  double lambda_eq = f_lambda(a_H2O);
  double Ceq = lambda_eq * C_SO3;
  double e_naf = (1-por)*eps / gamma;
  double S = - diffu * gamma / e_naf * (Ceq - C);
  return S;
}

// dSdc_vap = -Dw.gamma_CL/ep_naf.f_derivee_lambda(a).CSO3.R.T/Psat
double Young_Todd_Source_multi::f_dSdc_vap(double diffu, double ci, double C, double T, double por, double eps, double gamma) const
{
  double Psat = f_Psat(T);
  double a_H2O = ci * R * T / Psat;
  double d_lambda_eq = f_derivee_lambda(a_H2O);
  double e_naf = (1-por)*eps / gamma;
  double dSdc = - diffu * gamma / e_naf * d_lambda_eq * C_SO3 * R * T / Psat;
  return dSdc;
}
