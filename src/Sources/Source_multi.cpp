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
// File      : Source_multi.cpp
// Directory : $PEMFC_ROOT/src/Sources
//
/////////////////////////////////////////////////////////////////////////////

#include <Source_multi.h>
#include <Probleme_base.h>
#include <Conduction.h>
#include <Param.h>
#include <Interprete.h>
#include <Matrice_Morse.h>
#include <Champ_base.h>
#include <Domaine.h>
#include <Zone_VF.h>

Implemente_instanciable( Source_multi, "Source_multi_VEF_P1NC", Source_base ) ;

Sortie& Source_multi::printOn( Sortie& os ) const
{
  Source_base::printOn( os );
  return os;
}

Entree& Source_multi::readOn( Entree& is )
{
  Source_base::readOn( is );
  Cerr << " Source_multi::readOn " << finl  ;
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

// ajouter le terme2 = S_t
DoubleTab& Source_multi::ajouter(DoubleTab& resu) const
{
  assert(resu.dimension(0)==la_zone_.valeur().nb_faces());

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
              resu(face,0) +=  S_X2_(elem) * vol(elem) / nb_face_elem ;
              resu(face,1) +=  S_vap_(elem) * vol(elem) / nb_face_elem ;
              resu(face,2) +=  S_N2_(elem) * vol(elem) / nb_face_elem ;
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
                resu(face,0) +=  S_X2_(elem) * vol(elem) / nb_face_elem ;
                resu(face,1) +=  S_vap_(elem) * vol(elem) / nb_face_elem ;
                resu(face,2) +=  S_N2_(elem) * vol(elem) / nb_face_elem ;
              }
          }
      }
    }

  return resu;
}

DoubleTab& Source_multi::calculer(DoubleTab& resu) const
{
  resu = 0;
  ajouter(resu);
  return resu;
}

void Source_multi::mettre_a_jour(double temps)
{
  Cerr << "Source_multi::mettre_a_jour " << equation().probleme().le_nom() << finl;

  // update the involving fields -> to-do: nothing
  ch_ci_.valeur().mettre_a_jour(temps);

  const DoubleTab& xp=la_zone_.valeur().xp();	// centre gravite des elements
  ch_ci_.valeur().valeur_aux_compo(xp, c_X2_, 0);
  ch_ci_.valeur().valeur_aux_compo(xp, c_vap_, 1);
  ch_ci_.valeur().valeur_aux_compo(xp, c_N2_, 2);

  // mettre a jour les tableaux de S et dSdc
  const DoubleTab& por = por_naf_.valeurs();
  const DoubleTab& eps = eps_naf_.valeurs();
  const DoubleTab& gamma = gamma_CL_.valeurs();
  const DoubleTab& C_X2 = ch_C_X2_.valeurs();
  const DoubleTab& D_X2 = ch_D_X2_.valeurs();
  const DoubleTab& C_vap = ch_C_vap_.valeurs();
  const DoubleTab& D_vap = ch_D_vap_.valeurs();
  const DoubleTab& C_N2 = ch_C_N2_.valeurs();
  const DoubleTab& D_N2 = ch_D_N2_.valeurs();
  const DoubleTab& T = ch_T_.valeurs();

  if(CL_a_.non_nul())		// case anode H2 vap N2
    {
      for (int poly = 0; poly < CL_a_.valeur().nb_elem_tot(); poly++)
        {
          int elem = CL_a_.valeur()(poly);
          S_X2_(elem) = f_S_H2(D_X2(elem,0), c_X2_(elem), C_X2(elem,0), T(elem,0), por(elem,0), eps(elem, 0), gamma(elem,0));
          dSdc_X2_(elem) = f_dSdc_H2(D_X2(elem,0), c_X2_(elem), C_X2(elem,0), T(elem,0), por(elem,0), eps(elem, 0), gamma(elem,0));
          S_vap_(elem) = f_S_vap(D_vap(elem,0), c_vap_(elem), C_vap(elem,0), T(elem,0), por(elem,0), eps(elem, 0), gamma(elem,0));
          dSdc_vap_(elem) = f_dSdc_vap(D_vap(elem,0), c_vap_(elem), C_vap(elem,0), T(elem,0), por(elem,0), eps(elem, 0), gamma(elem,0));
          S_N2_(elem) = f_S_N2(D_N2(elem,0), c_N2_(elem), C_N2(elem,0), T(elem,0), por(elem,0), eps(elem, 0), gamma(elem,0));
          dSdc_N2_(elem) = f_dSdc_N2(D_N2(elem,0), c_N2_(elem), C_N2(elem,0), T(elem,0), por(elem,0), eps(elem, 0), gamma(elem,0));
        }
    }
  if(CL_c_.non_nul())
    {
      {
        for (int poly = 0; poly < CL_c_.valeur().nb_elem_tot(); poly++)
          {
            int elem = CL_c_.valeur()(poly);
            S_X2_(elem) = f_S_O2(D_X2(elem,0), c_X2_(elem), C_X2(elem,0), T(elem,0), por(elem,0), eps(elem, 0), gamma(elem,0));
            dSdc_X2_(elem) = f_dSdc_O2(D_X2(elem,0), c_X2_(elem), C_X2(elem,0), T(elem,0), por(elem,0), eps(elem, 0), gamma(elem,0));
            S_vap_(elem) = f_S_vap(D_vap(elem,0), c_vap_(elem), C_vap(elem,0), T(elem,0), por(elem,0), eps(elem, 0), gamma(elem,0));
            dSdc_vap_(elem) = f_dSdc_vap(D_vap(elem,0), c_vap_(elem), C_vap(elem,0), T(elem,0), por(elem,0), eps(elem, 0), gamma(elem,0));
            S_N2_(elem) = f_S_N2(D_N2(elem,0), c_N2_(elem), C_N2(elem,0), T(elem,0), por(elem,0), eps(elem, 0), gamma(elem,0));
            dSdc_N2_(elem) = f_dSdc_N2(D_N2(elem,0), c_N2_(elem), C_N2(elem,0), T(elem,0), por(elem,0), eps(elem, 0), gamma(elem,0));
          }
      }
    }
  Cerr << "S_X2 min max " << mp_min_vect(S_X2_) << " " << mp_max_vect(S_X2_) << finl;
  Cerr << "S_vap min max " << mp_min_vect(S_vap_) << " " << mp_max_vect(S_vap_) << finl;
  Cerr << "S_N2 min max " << mp_min_vect(S_N2_) << " " << mp_max_vect(S_N2_) << finl;
  Cerr << "dSdc_X2 min max " << mp_min_vect(dSdc_X2_) << " " << mp_max_vect(dSdc_X2_) << finl;
  Cerr << "dSdc_vap min max " << mp_min_vect(dSdc_vap_) << " " << mp_max_vect(dSdc_vap_) << finl;
  Cerr << "dSdc_N2 min max " << mp_min_vect(dSdc_N2_) << " " << mp_max_vect(dSdc_N2_) << finl;
}

// ajouter terme1 = -dS/dc
void Source_multi::contribuer_a_avec(const DoubleTab& inco, Matrice_Morse& mat) const
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
              mat.coef(face * nb_comp,face * nb_comp) +=  - dSdc_X2_(elem) * vol(elem) / nb_face_elem ;
              mat.coef(face*nb_comp+1,face*nb_comp+1) +=  - dSdc_vap_(elem) * vol(elem) / nb_face_elem ;
              mat.coef(face*nb_comp+2,face*nb_comp+2) +=  - dSdc_N2_(elem) * vol(elem) / nb_face_elem ;
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
                mat.coef(face * nb_comp,face * nb_comp) +=  - dSdc_X2_(elem) * vol(elem) / nb_face_elem ;
                mat.coef(face*nb_comp+1,face*nb_comp+1) +=  - dSdc_vap_(elem) * vol(elem) / nb_face_elem ;
                mat.coef(face*nb_comp+2,face*nb_comp+2) +=  - dSdc_N2_(elem) * vol(elem) / nb_face_elem ;
              }
          }
      }
    }
}

void Source_multi::associer_zones(const Zone_dis& zone_dis, const Zone_Cl_dis& zcl_dis)
{
  la_zone_ = ref_cast(Zone_VF,zone_dis.valeur());
}

void Source_multi::associer_pb(const Probleme_base& pb)
{
  // nothing to be associated
}

void Source_multi::completer()
{
  Source_base::completer();
  // get the reference to the coupling fields
  ch_ci_ = equation().inconnue();

  assert(ch_ci_.valeur().valeurs().dimension(1) == 3);
  la_zone_.valeur().zone().creer_tableau_elements(c_X2_);
  la_zone_.valeur().zone().creer_tableau_elements(c_vap_);
  la_zone_.valeur().zone().creer_tableau_elements(c_N2_);
  la_zone_.valeur().zone().creer_tableau_elements(S_X2_);
  la_zone_.valeur().zone().creer_tableau_elements(S_vap_);
  la_zone_.valeur().zone().creer_tableau_elements(S_N2_);
  la_zone_.valeur().zone().creer_tableau_elements(dSdc_X2_);
  la_zone_.valeur().zone().creer_tableau_elements(dSdc_vap_);
  la_zone_.valeur().zone().creer_tableau_elements(dSdc_N2_);
}

void Source_multi::set_param(Param& param)
{
  param.ajouter("nom_domaine",&nom_dom_,Param::REQUIRED);
  param.ajouter("nom_ssz_CLa",&nom_ssz_CLa_,Param::OPTIONAL);	// requis pour anode
  param.ajouter("nom_ssz_CLc",&nom_ssz_CLc_,Param::OPTIONAL);	// requis pour cathode
  param.ajouter("por_naf", &por_naf_, Param::REQUIRED);
  param.ajouter("eps_naf", &eps_naf_, Param::REQUIRED);
  param.ajouter("gamma_CL", &gamma_CL_, Param::REQUIRED);

  param.ajouter("ch_C_X2", &ch_C_X2_, Param::REQUIRED);
  param.ajouter("ch_D_X2", &ch_D_X2_, Param::REQUIRED);
  param.ajouter("ch_C_vap", &ch_C_vap_, Param::REQUIRED);
  param.ajouter("ch_D_vap", &ch_D_vap_, Param::REQUIRED);
  param.ajouter("ch_C_N2", &ch_C_N2_, Param::REQUIRED);
  param.ajouter("ch_D_N2", &ch_D_N2_, Param::REQUIRED);
  param.ajouter("ch_T", &ch_T_, Param::REQUIRED);
}

// S_X2 = -D.gamma_CL/ep_naf(cRTH - C)
double Source_multi::f_S_H2(double diffu, double ci, double C, double T, double por, double eps, double gamma) const
{
  double H = f_Henry_H2(T);
  double Ceq = ci * R * T * H;
  double e_naf = (1-por)*eps / gamma;
  double S = - diffu * gamma / e_naf * (Ceq - C);
  return S;
}

// dSdc_X2 = -D.gamma_CL/ep_naf.RTH
double Source_multi::f_dSdc_H2(double diffu, double ci, double C, double T, double por, double eps, double gamma) const
{
  double H = f_Henry_H2(T);
  double RTH = R * T * H;
  double e_naf = (1-por)*eps / gamma;
  double dSdc = - diffu * gamma / e_naf * RTH;
  return dSdc;
}

// S_O2 = -D.gamma_CL/ep_naf(cRTH - C)
double Source_multi::f_S_O2(double diffu, double ci, double C, double T, double por, double eps, double gamma) const
{
  double H = f_Henry_O2(T);
  double Ceq = ci * R * T * H;
  double e_naf = (1-por)*eps / gamma;
  double S = - diffu * gamma / e_naf * (Ceq - C);
  return S;
}

// dSdc_O2 = -D.gamma_CL/ep_naf.RTH
double Source_multi::f_dSdc_O2(double diffu, double ci, double C, double T, double por, double eps, double gamma) const
{
  double H = f_Henry_O2(T);
  double RTH = R * T * H;
  double e_naf = (1-por)*eps / gamma;
  double dSdc = - diffu * gamma / e_naf * RTH;
  return dSdc;
}

// S_N2 = -D.gamma_CL/ep_naf(cRTH - C)
double Source_multi::f_S_N2(double diffu, double ci, double C, double T, double por, double eps, double gamma) const
{
  double H = f_Henry_N2(T);
  double Ceq = ci * R * T * H;
  double e_naf = (1-por)*eps / gamma;
  double S = - diffu * gamma / e_naf * (Ceq - C);
  return S;
}

// dSdc_N2 = -D.gamma_CL/ep_naf.RTH
double Source_multi::f_dSdc_N2(double diffu, double ci, double C, double T, double por, double eps, double gamma) const
{
  double H = f_Henry_N2(T);
  double RTH = R * T * H;
  double e_naf = (1-por)*eps / gamma;
  double dSdc = - diffu * gamma / e_naf * RTH;
  return dSdc;
}

// S_vap = -Dw.gamma_CL/ep_naf(Ceq - C) avec Ceq = CSO3.f_lambda(a) et a=cRT/Psat
double Source_multi::f_S_vap(double diffu, double ci, double C, double T, double por, double eps, double gamma) const
{
  double Psat = f_Psat(T);
  double a_H2O = ci * R * T / Psat;
  double lambda_eq = f_lambda(a_H2O);
  double Ceq = lambda_eq * C_SO3;
  double e_naf = (1-por)*eps / gamma;
  double S = - diffu * gamma / e_naf * (Ceq - C);
  return S;
}

// dSdc_vap = -Dw.gamma_CL/ep_naf.f_derivee_lambda(a).RT/Psat
double Source_multi::f_dSdc_vap(double diffu, double ci, double C, double T, double por, double eps, double gamma) const
{
  double Psat = f_Psat(T);
  double a_H2O = ci * R * T / Psat;
  double d_lambda_eq = f_derivee_lambda(a_H2O);
  double e_naf = (1-por)*eps / gamma;
  double dSdc = - diffu * gamma / e_naf * C_SO3 * d_lambda_eq * R * T / Psat;

  return dSdc;
}
