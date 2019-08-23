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
// File      : Source_Term_pemfc_base.cpp
// Directory : $PEMFC_ROOT/src/Sources
//
/////////////////////////////////////////////////////////////////////////////

#include <Source_Term_pemfc_base.h>
#include <Probleme_base.h>
#include <Conduction.h>
#include <Param.h>
#include <Interprete.h>
#include <Matrice_Morse.h>
#include <Champ_base.h>
#include <Champ_Inc.h>
#include <Domaine.h>


Implemente_base( Source_Term_pemfc_base, "Source_Term_pemfc_base", Source_base ) ;

Sortie& Source_Term_pemfc_base::printOn( Sortie& os ) const
{
  Source_base::printOn( os );
  return os;
}

Entree& Source_Term_pemfc_base::readOn( Entree& is )
{
  //Source_base::readOn( is );
  Cerr << " Source_Term_pemfc_base::readOn " << finl  ;
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

  if(nom_espece_ == "H2")
    {
      // check
      assert(nom_ssz_CLa_ != "??");
      assert(nom_pb_ci_anode_ != "??");
      assert(nom_champ_ci_anode_ != "??");
    }
  else if (nom_espece_ == "O2")
    {
      assert(nom_ssz_CLc_ != "??");
      assert(nom_pb_ci_cathode_ != "??");
      assert(nom_champ_ci_cathode_ != "??");
    }
  else
    {
      assert(nom_ssz_CLa_ != "??");
      assert(nom_ssz_CLc_ != "??");
      assert(nom_pb_ci_cathode_ != "??");
      assert(nom_champ_ci_cathode_ != "??");
      assert(nom_pb_ci_anode_ != "??");
      assert(nom_champ_ci_anode_ != "??");
//      if (nom_espece_ == "H2O" || nom_espece_ == "vap")
//        {
//          assert(nom_pb_phi_ != "??");
//          assert(nom_champ_ir_ != "??");
//          assert(nom_champ_ip_ != "??");
//          assert(nom_op_diff_ != "??");
//        }
    }

  if(nom_ssz_CLa_ != "??")
    CL_a_ = dom_.valeur().ss_zone(nom_ssz_CLa_);
  if(nom_ssz_CLc_ != "??")
    CL_c_ = dom_.valeur().ss_zone(nom_ssz_CLc_);

  // discretiser le champ source -> pour post-traiter (diffusion des multi-especes dans le milieu poreux)
  discretiser(equation().discretisation());

  return is;
}


void Source_Term_pemfc_base::set_param(Param& param)
{
  param.ajouter("nom_espece",&nom_espece_,Param::REQUIRED);
  param.ajouter("nom_domaine",&nom_domaine_,Param::REQUIRED);		// pas necessaire ? dom_ = equation().problem().domaine() ?
  param.ajouter("nom_CLa",&nom_ssz_CLa_,Param::OPTIONAL);					// requis pour H2, N2, H20
  param.ajouter("nom_CLc",&nom_ssz_CLc_,Param::OPTIONAL);					// requis pour O2, N2, H20
  param.ajouter("por_naf", &por_naf_, Param::REQUIRED);
  param.ajouter("eps_naf", &eps_naf_, Param::REQUIRED);
  param.ajouter("gamma_CL", &gamma_CL_, Param::REQUIRED);
  param.ajouter("C_SO3", &C_SO3_, Param::REQUIRED);
  param.ajouter("nom_champ_D", &nom_champ_D_, Param::REQUIRED);			// diffusion_nafion requis si couplage
  param.ajouter("nom_pb_T", &nom_pb_T_, Param::REQUIRED);
  param.ajouter("nom_champ_T", &nom_champ_T_, Param::REQUIRED);
  param.ajouter("nom_pb_ci_cathode", &nom_pb_ci_cathode_, Param::OPTIONAL);
  param.ajouter("nom_champ_ci_cathode", &nom_champ_ci_cathode_, Param::OPTIONAL);
  param.ajouter("nom_pb_ci_anode", &nom_pb_ci_anode_, Param::OPTIONAL);
  param.ajouter("nom_champ_ci_anode", &nom_champ_ci_anode_, Param::OPTIONAL);
  //param.ajouter("nom_pb_phi", &nom_pb_phi_, Param::OPTIONAL);
  //param.ajouter("nom_champ_ir", &nom_champ_ir_, Param::OPTIONAL);
  //param.ajouter("nom_champ_ip", &nom_champ_ip_, Param::OPTIONAL);
  //param.ajouter("nom_op_diff", &nom_op_diff_, Param::OPTIONAL);
}

void Source_Term_pemfc_base::discretiser(const Discretisation_base& dis)
{
  dis.discretiser_champ("temperature",equation().zone_dis().valeur(),"S","unit", 1 , 0. , ch_S_);
  champs_compris_.ajoute_champ(ch_S_);
}

double Source_Term_pemfc_base::eval_f(double diffu, double Ci, double ci, double T) const
{
  double Ceq;
  if (nom_espece_ == "H2O" || nom_espece_ == "vap")
    {
      double Psat = f_Psat(T);
      double a_H20 = ci * R * T / Psat;
      double lambda_eq = f_lambda(a_H20);
      Ceq = lambda_eq * C_SO3_;
    }
  else
    {
      double H = 0.;
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
      Ceq = ci * R * T * H;
    }
  double e_naf = (1-por_naf_)*eps_naf_ / gamma_CL_;
  return diffu * gamma_CL_ / e_naf * (Ceq - Ci);
}

double Source_Term_pemfc_base::eval_derivee_f(double diffu) const
{
  // expression_derivee_par_rapport_inconnue
  double e_naf = (1-por_naf_)*eps_naf_ / gamma_CL_;
  return (- diffu * gamma_CL_ / e_naf);
}

DoubleTab& Source_Term_pemfc_base::calculer(DoubleTab& resu) const
{
  resu = 0;
  ajouter(resu);
  return resu;
}

void Source_Term_pemfc_base::mettre_a_jour(double temps)
{
  // update the involving fields -> to-do: nothing
  if(ch_T_.non_nul())
    ch_T_.valeur().mettre_a_jour(temps);
  if(ch_ci_cathode_.non_nul())
    ch_ci_cathode_.valeur().mettre_a_jour(temps);
  if(ch_ci_anode_.non_nul())
    ch_ci_anode_.valeur().mettre_a_jour(temps);
//  if(ch_ir_.non_nul())
//    ch_ir_.valeur().mettre_a_jour(temps);
//  if(ch_ip_.non_nul())
//    ch_ip_.valeur().mettre_a_jour(temps);
//  if(ch_op_.non_nul())
//    ch_op_.valeur().mettre_a_jour(temps);
}

void Source_Term_pemfc_base::completer()
{
  Source_base::completer();
  // get the reference to the coupling fields
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
//  if(nom_pb_phi_ != "??" )
//    {
//      Probleme_base& pb_phi = ref_cast(Probleme_base,interprete().objet(nom_pb_phi_));
//      ch_ir_ = pb_phi.get_champ(nom_champ_ir_);
//      ch_ir_ = pb_phi.get_champ(nom_champ_ir_);
//      ch_op_ = pb_phi.get_champ(nom_op_diff_);
//    }
  ch_C_ = equation().inconnue();
  ch_D_i_naf_ = equation().probleme().get_champ(nom_champ_D_);
}
