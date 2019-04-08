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
// File      : Loi_Fermeture_PEMFC_Cas2.cpp
// Directory : $PEMFC_ROOT/src
//
/////////////////////////////////////////////////////////////////////////////

#include <Loi_Fermeture_PEMFC_Cas2.h>
#include <Probleme_base.h>
#include <Discretisation_base.h>
#include <Equation_base.h>

Implemente_instanciable( Loi_Fermeture_PEMFC_Cas2, "Loi_Fermeture_PEMFC_Cas2", Loi_Fermeture_PEMFC_base ) ;
// XD loi_fermeture_pemfc_cas2 loi_fermeture_base loi_fermeture_pemfc_cas2 -1 Loi for test only

Sortie& Loi_Fermeture_PEMFC_Cas2::printOn( Sortie& os ) const
{
  Loi_Fermeture_PEMFC_base::printOn( os );
  return os;
}

Entree& Loi_Fermeture_PEMFC_Cas2::readOn( Entree& is )
{
  Loi_Fermeture_PEMFC_base::readOn( is );
  return is;
}

void Loi_Fermeture_PEMFC_Cas2::discretiser(const Discretisation_base& dis)
{
  Loi_Fermeture_PEMFC_base::discretiser(dis);
  ref_equation_=mon_probleme().get_equation_by_name("Convection_diffusion_Concentration");

  dis.discretiser_champ("champ_elem",equation().zone_dis().valeur(),"pemfc_cas2","unit", equation().inconnue().valeur().nb_comp()*equation().inconnue().valeur().nb_comp(),0.,diffu_);
  champs_compris_.ajoute_champ(diffu_);
  diffu_.valeur().fixer_nature_du_champ(multi_scalaire);

// READON optionel
//  status_ = READON_FAIT;
}
