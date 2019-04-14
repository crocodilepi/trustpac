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
}

void Source_Term_pemfc_VDF_P0_VDF::associer_pb(const Probleme_base& pb)
{

}

void Source_Term_pemfc_VDF_P0_VDF::remplir_volumes()
{
  volumes_.ref(ref_cast(Zone_VF,equation().zone_dis().valeur()).volumes_entrelaces());
}

