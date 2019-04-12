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
// File      : Source_Term_pemfc_VEF_P1NC.cpp
// Directory : $PEMFC_ROOT/src/Sources
//
/////////////////////////////////////////////////////////////////////////////
#include <Source_Term_pemfc_VEF_P1NC.h>
#include <Zone_VEF.h>
#include <Zone_Cl_VEF.h>
#include <Probleme_base.h>
#include <Equation_base.h>
#include <Conduction.h>
#include <Champ_P1NC.h>

Implemente_instanciable( Source_Term_pemfc_VEF_P1NC, "Source_Term_pemfc_VEF_P1NC", Source_Term_pemfc_base );

Sortie& Source_Term_pemfc_VEF_P1NC::printOn(Sortie& os) const
{
  Source_Term_pemfc_base::printOn(os);
  return os;
}

Entree& Source_Term_pemfc_VEF_P1NC::readOn(Entree& is)
{
  Source_Term_pemfc_base::readOn(is);
  return is;
}

void Source_Term_pemfc_VEF_P1NC::associer_zones(const Zone_dis& zone_dis, const Zone_Cl_dis& zcl_dis)
{
  la_zone_VEF = ref_cast(Zone_VEF,zone_dis.valeur());
  la_zcl_VEF = ref_cast(Zone_Cl_VEF,zcl_dis.valeur());
}

void Source_Term_pemfc_VEF_P1NC::associer_pb(const Probleme_base& pb)
{
  Cerr << " Source_Term_pemfc_VEF_P1NC::associer_pb " << finl ;
  int ok = 0;
  const Equation_base& eqn = pb.equation(0);
  if  (sub_type(Conduction,eqn))
    {
      C_ = ref_cast(Champ_P1NC,eqn.inconnue().valeur());
      Da_ = ref_cast(Champ_P1NC,pb.get_champ("conductivite"));
      associer_zones(eqn.zone_dis(),eqn.zone_Cl_dis());
      ok = 1;
    }
  if (!ok)
    {
      Cerr << "Erreur TRUST dans Source_Term_pemfc_VEF_P1NC::associer_pb()" << finl;
      Cerr << "On ne trouve pas d'equation de conduction dans le probleme" << finl;
      exit();
    }

}

void Source_Term_pemfc_VEF_P1NC::remplir_volumes()
{
  volumes_.ref(la_zone_VEF.valeur().volumes_entrelaces());
}
