/****************************************************************************
* Copyright (c) 2015, CEA
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
// File      : Op_Diff_VEF_Face_PEMFC.cpp
// Directory : $PEMFC_ROOT/src/Cas2
//
/////////////////////////////////////////////////////////////////////////////

#include <Op_Diff_VEF_Face_PEMFC.h>
#include <Param.h>
#include <EChaine.h>
#include <Discretisation_base.h>
#include <Champ_P1NC.h>
#include <Milieu_base.h>
#include <Probleme_base.h>

Implemente_instanciable( Op_Diff_VEF_Face_PEMFC, "Op_Diff_VEFPEMFC_const_P1NC", Op_Diff_VEF_Face_Matricial ) ;
// XD diffusion_pemfc diffusion_deriv pemfc 1 diffusion (only in VEF)

Sortie& Op_Diff_VEF_Face_PEMFC::printOn( Sortie& os ) const
{
  Op_Diff_VEF_Face_Matricial::printOn( os );
  return os;
}

Entree& Op_Diff_VEF_Face_PEMFC::readOn( Entree& is )
{
  Param param(que_suis_je());
  set_param(param);
  param.lire_avec_accolades_depuis(is);

  const Champ_base& diffu=equation().probleme().get_champ(diffu_name_);
  associer_diffusivite ( diffu );

  return is;
}
void Op_Diff_VEF_Face_PEMFC::set_param(Param& param)
{
  param.ajouter("diffusivity_fieldname",&diffu_name_,Param::REQUIRED); // XD_ADD_P chaine name of the diffusity field
}




void Op_Diff_VEF_Face_PEMFC::completer()
{
  Op_Diff_VEF_Face_Matricial::completer();

  // throw;
}
