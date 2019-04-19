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
// File      : Source_Term_pemfc_VDF_P0_VDF.cpp
// Directory : $PEMFC_ROOT/src/Sources
//
/////////////////////////////////////////////////////////////////////////////
#include <Source_Term_pemfc_VDF_P0_VDF.h>
#include <Zone_VDF.h>
#include <Zone_Cl_VDF.h>
#include <Equation_base.h>
#include <Champ_Generique_base.h>

Implemente_instanciable( Source_Term_pemfc_VDF_P0_VDF, "Source_Term_pemfc_VDF_P0_VDF", Source_Term_pemfc_base );

Sortie& Source_Term_pemfc_VDF_P0_VDF::printOn(Sortie& os) const
{
  Source_Term_pemfc_base::printOn(os);
  return os;
}

Entree& Source_Term_pemfc_VDF_P0_VDF::readOn(Entree& is)
{
  Source_Term_pemfc_base::readOn(is);
  return is;
}

void Source_Term_pemfc_VDF_P0_VDF::associer_zones(const Zone_dis& zone_dis, const Zone_Cl_dis& zcl_dis)
{
  la_zone_VDF = ref_cast(Zone_VDF,zone_dis.valeur());
  la_zcl_VDF = ref_cast(Zone_Cl_VDF,zcl_dis.valeur());
  remplir_volumes();
}

void Source_Term_pemfc_VDF_P0_VDF::associer_pb(const Probleme_base& pb)
{
  Cerr << " Source_Term_pemfc_VDF_P0_VDF::associer_pb " << finl ;
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
      Cerr << "Erreur TRUST dans Source_Term_pemfc_VDF_P0_VDF::associer_pb()" << finl;
      Cerr << "On ne trouve pas d'equation de conduction dans le probleme" << finl;
      exit();
    }

}

DoubleTab& Source_Term_pemfc_VDF_P0_VDF::ajouter(DoubleTab& resu) const
{
  assert(resu.dimension(0)==volumes_.size());
  assert(resu.dimension(0)==T_.size());
  assert(resu.dimension(0)==C_.size());
  assert(resu.dimension(0)==diffu_.size());
  assert(resu.dimension(0)==c_.size());

  double inv_rhoCp = 1./((1-epsilon_)*epsilon_ionomer_);
  int poly;
  for (poly = 0; poly < ssz_.valeur().nb_elem_tot(); poly++)
    {
      int elem = ssz_.valeur()(poly);
      resu(elem) = eval_f(diffu_(elem), C_(elem), c_(elem), T_(elem)) * inv_rhoCp;
    }
  return resu;
}

void Source_Term_pemfc_VDF_P0_VDF::contribuer_a_avec(const DoubleTab& inco, Matrice_Morse& mat) const
{
  assert(inco.dimension(0)==volumes_.size());
  assert(inco.dimension(0)==diffu_.size());

  double inv_rhoCp = 1./((1-epsilon_)*epsilon_ionomer_);
  int poly;
  for (poly = 0; poly < ssz_.valeur().nb_elem_tot(); poly++)
    {
      int elem = ssz_.valeur()(poly);
      mat.coef(elem,elem) += volumes_(elem) * eval_derivee_f(diffu_(elem))*inv_rhoCp;
    }
}

void Source_Term_pemfc_VDF_P0_VDF::completer()
{
  Source_Term_pemfc_base::completer();

  // recuperer les references vers les champs
  Champ stoT;
  const Champ_base& ch_T = equation().probleme().get_champ_post(nom_champ_T_).get_champ(stoT);
  assert(ch_T.que_suis_je().find("P0") !=-1);
  T_.ref(ch_T.valeurs());

  Champ stoc;
  const Champ_base& ch_c = equation().probleme().get_champ_post(nom_champ_c_).get_champ(stoc);
  assert(ch_c.que_suis_je().find("P0") !=-1);
  c_.ref(ch_c.valeurs());

  Champ_base& ch_C = equation().inconnue();
  if(ch_C.que_suis_je().find("P0") !=-1)
    {
      C_.ref(ch_C.valeurs());
    }
  else
    {
      Champ stoC;
      const Champ_base& ch_C_post = equation().probleme().get_champ_post(nom_champ_C_).get_champ(stoC);
      assert(ch_C_post.que_suis_je().find("P0") !=-1);
      C_.ref(ch_C_post.valeurs());
    }

  const Champ_base& ch_D = equation().probleme().get_champ("diffusion_nafion");
  if(ch_D.que_suis_je().find("P0") !=-1)
    {
      diffu_.ref(ch_D.valeurs());
    }
  else
    {
      Champ stoD;
      const Champ_base& ch_D_post = equation().probleme().get_champ_post(nom_champ_D_).get_champ(stoD);
      assert(ch_D_post.que_suis_je().find("P0") !=-1);
      diffu_.ref(ch_D_post.valeurs());
    }
}

void Source_Term_pemfc_VDF_P0_VDF::remplir_volumes()
{
  volumes_.ref(ref_cast(Zone_VF,equation().zone_dis().valeur()).volumes_entrelaces());
}
