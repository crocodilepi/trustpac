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
// File      : Source_Term_Nafion_Reaction.cpp
// Directory : $PEMFC_ROOT/src/Sources
//
/////////////////////////////////////////////////////////////////////////////

#include <Source_Term_Nafion_Reaction.h>
#include <Probleme_base.h>
#include <Equation_base.h>
#include <Conduction.h>
#include <Champ_P1NC.h>
#include <Zone_VEF.h>
#include <Zone_Cl_VEF.h>
#include <Domaine.h>
#include <Champ.h>
#include <Champ_Generique_base.h>

Implemente_instanciable( Source_Term_Nafion_Reaction, "Source_Term_Nafion_Reaction_VEF_P1NC", Source_base ) ;

Sortie& Source_Term_Nafion_Reaction::printOn( Sortie& os ) const
{
  Source_base::printOn( os );
  return os;
}

Entree& Source_Term_Nafion_Reaction::readOn( Entree& is )
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
  dom_ = equation().probleme().domaine();
  ssz_ = dom_.valeur().ss_zone(nom_ssz_);
  return is;
  return is;
}

void Source_Term_Nafion_Reaction::associer_pb(const Probleme_base& pb)
{
  Cerr << " Source_Term_Nafion_Reaction::associer_pb " << finl ;
  assert(pb.que_suis_je() == "Pb_Conduction");
  int ok = 0;
  const Equation_base& eqn = pb.equation(0);
  if  (eqn.que_suis_je() == "Conduction")
    {
      associer_zones(eqn.zone_dis(),eqn.zone_Cl_dis());
      ok = 1;
    }
  if (!ok)
    {
      Cerr << "Erreur TRUST dans Source_Term_Nafion_Reaction::associer_pb()" << finl;
      Cerr << "On ne trouve pas d'equation de conduction dans le probleme" << finl;
      exit();
    }
}

void Source_Term_Nafion_Reaction::associer_zones(const Zone_dis& zone_dis, const Zone_Cl_dis& zcl_dis)
{
  // cast VEF
  la_zone_VEF = ref_cast(Zone_VEF,zone_dis.valeur());
  la_zcl_VEF = ref_cast(Zone_Cl_VEF,zcl_dis.valeur());
}

void Source_Term_Nafion_Reaction::set_param(Param& param)
{
  param.ajouter("nom_espece",&nom_espece_,Param::REQUIRED);
  param.ajouter("nom_domaine",&nom_domaine_,Param::REQUIRED);		// pas necessaire ? dom_ = equation().problem().domaine() ?
  param.ajouter("nom_ssz",&nom_ssz_,Param::REQUIRED);
  param.ajouter("nom_champ_jr",&nom_champ_jr_,Param::REQUIRED);
  param.ajouter("nom_champ_jp",&nom_champ_jp_,Param::REQUIRED);
}

void Source_Term_Nafion_Reaction::completer()
{
  Champ stoJr;
  const Champ_base& ch_jr = equation().probleme().get_champ_post(nom_champ_jr_).get_champ(stoJr);
  assert(ch_jr.que_suis_je().find("P1NC") !=-1);
  jr_.ref(ch_jr.valeurs());

  Champ stoJp;
  const Champ_base& ch_jp = equation().probleme().get_champ_post(nom_champ_jp_).get_champ(stoJp);
  assert(ch_jp.que_suis_je().find("P1NC") !=-1);
  jp_.ref(ch_jp.valeurs());
}

void Source_Term_Nafion_Reaction::mettre_a_jour(double temps)
{
  // update the involving fields -> to-do: nothing
}

DoubleTab& Source_Term_Nafion_Reaction::calculer(DoubleTab& resu) const
{
  resu = 0;
  ajouter(resu);
  return resu;
}

void Source_Term_Nafion_Reaction::remplir_volumes()
{
  volumes_.ref(la_zone_VEF.valeur().volumes_entrelaces());
}

DoubleTab& Source_Term_Nafion_Reaction::ajouter(DoubleTab& resu) const
{

  // to-do

  return resu;
}

double Source_Term_Nafion_Reaction::eval_f(double jr, double jp) const
{
  double F = 96500;
  if (nom_espece_ == "H2")
    {
      return - jr / (2. * F);
    }
  else if (nom_espece_ == "02")
    {
      return (jr + jp) / (4. * F);
    }
  else if (nom_espece_ == "N2")
    {
      return 0.;
    }
  else if (nom_espece_ == "H2O" || nom_espece_ == "vap")
    {
      return -(jr + jp)/(2. * F);
    }
  return 0.;
}
