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
// File      : Loi_Fermeture_Diffusion_Nafion.cpp
// Directory : $PEMFC_ROOT/src/Loi_Fermeture
//
/////////////////////////////////////////////////////////////////////////////

#include <Loi_Fermeture_Diffusion_Nafion.h>
#include <Pb_Conduction.h>
#include <Param.h>
#include <Discretisation_base.h>
#include <Champ_base.h>
#include <Champ_Generique_base.h>

Implemente_instanciable( Loi_Fermeture_Diffusion_Nafion, "Loi_Fermeture_Diffusion_Nafion", Loi_Fermeture_base ) ;

Sortie& Loi_Fermeture_Diffusion_Nafion::printOn( Sortie& os ) const
{
  Loi_Fermeture_base::printOn( os );
  return os;
}

Entree& Loi_Fermeture_Diffusion_Nafion::readOn( Entree& is )
{
  Loi_Fermeture_base::readOn( is );

  // check param input
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
  return is;
}

void Loi_Fermeture_Diffusion_Nafion::associer_pb_base(const Probleme_base& pb)
{
  assert(sub_type(Pb_Conduction, pb));
  Loi_Fermeture_base::associer_pb_base(pb);
}

void Loi_Fermeture_Diffusion_Nafion::discretiser(const Discretisation_base& dis)
{
  Loi_Fermeture_base::discretiser(dis);
  ref_equation_=mon_probleme().get_equation_by_name("Conduction");
  dis.discretiser_champ("champ_elem",equation().zone_dis().valeur(),"diffusion_nafion","unit", 1 ,0. , diffu_);
  champs_compris_.ajoute_champ(diffu_);
}

void Loi_Fermeture_Diffusion_Nafion::completer()
{
  Loi_Fermeture_base::completer();

  if(nom_champ_T_ != "??")
    {
      Champ sto;
      const Champ_base& ch_T = equation().probleme().get_champ_post(nom_champ_T_).get_champ(sto);
      assert(ch_T.que_suis_je().find("P0") !=-1);
      T_.ref(ch_T.valeurs());
    }
  else
    {
      equation().zone_dis().zone().creer_tableau_elements(T_);
      T_ = default_value_T_;
    }
  if(nom_champ_C_ != "??")
    {
      Champ sto;
      const Champ_base& ch_C = equation().probleme().get_champ_post(nom_champ_C_).get_champ(sto);
      assert(ch_C.que_suis_je().find("P0") !=-1);
      C_.ref(ch_C.valeurs());
    }
  else
    {
      equation().zone_dis().zone().creer_tableau_elements(T_);
      C_ = CSO3_;
    }
}

void Loi_Fermeture_Diffusion_Nafion::preparer_calcul()
{
  Loi_Fermeture_base::preparer_calcul();
}

void Loi_Fermeture_Diffusion_Nafion::mettre_a_jour(double temps)
{
  // mettre a jour le champ diffusivite
  DoubleTab& diffu = diffu_.valeur().valeurs();
  int nb_elem = diffu.dimension(0);
  for (int elem = 0; elem < nb_elem; ++elem)
    {
      diffu(elem) = eval(T_(elem), C_(elem));
    }
  Cerr << "Loi_Fermeture_Diffusion_Nafion::mettre_a_jour" << finl;
  Cerr << "diffu min max " << mp_min_vect(diffu)<< " " <<mp_max_vect(diffu) << finl;
}

void Loi_Fermeture_Diffusion_Nafion::set_param(Param& param)
{
  param.ajouter("nom_espece",&nom_espece_,Param::REQUIRED);
  param.ajouter("default_value_T", &default_value_T_, Param::REQUIRED);
  param.ajouter("nom_champ_T",&nom_champ_T_,Param::OPTIONAL);
  param.ajouter("nom_champ_C",&nom_champ_C_,Param::OPTIONAL);
  param.ajouter("CSO3", &CSO3_, Param::OPTIONAL);
}

double Loi_Fermeture_Diffusion_Nafion::eval(double T, double C)
{
  if (nom_espece_ == "H2")
    {
      return 4.1e-7*exp((-2602.)/T);
    }
  else if (nom_espece_ == "02")
    {
      return 4.24e-6*exp(-2246. / T);
    }
  else if (nom_espece_ == "N2")
    {
      return 3.1e-7 * exp(-2768. /T);
    }
  else if (nom_espece_ == "H2O" || nom_espece_ == "vap")
    {
      double lambda = C / CSO3_;
      return (6.707e-8*lambda + 6.387e-7)*exp(-2416. / T);
    }
  // pour compilos.
  return 0.;
}
