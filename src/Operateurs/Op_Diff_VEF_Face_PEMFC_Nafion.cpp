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
// File      : Op_Diff_VEF_Face_PEMFC_Nafion.cpp
// Directory : $PEMFC_ROOT/src/Operateurs
//
/////////////////////////////////////////////////////////////////////////////

#include <Op_Diff_VEF_Face_PEMFC_Nafion.h>
#include <Champ_base.h>
#include <Probleme_base.h>
#include <Champ_Uniforme.h>

Implemente_instanciable( Op_Diff_VEF_Face_PEMFC_Nafion, "Op_Diff_VEFNAFION_const_P1NC", Op_Diff_VEF_Face ) ;

Sortie& Op_Diff_VEF_Face_PEMFC_Nafion::printOn( Sortie& os ) const
{
  Op_Diff_VEF_Face::printOn( os );
  return os;
}

Entree& Op_Diff_VEF_Face_PEMFC_Nafion::readOn( Entree& is )
{
  Param param(que_suis_je());
  set_param(param);
  param.lire_avec_accolades_depuis(is);

  const Champ_base& diffu=equation().probleme().get_champ(diffu_name_);
  associer_diffusivite ( diffu );

  return is;
}

void Op_Diff_VEF_Face_PEMFC_Nafion::set_param(Param& param)
{
  param.ajouter("diffusivity_fieldname",&diffu_name_,Param::REQUIRED); // XD_ADD_P chaine name of the diffusity field
  param.ajouter("por_naf", &por_naf_, Param::REQUIRED);
  param.ajouter("eps_naf", &eps_naf_, Param::REQUIRED);
}

void Op_Diff_VEF_Face_PEMFC_Nafion::completer()
{
  Op_Diff_VEF_Face::completer();
}

void Op_Diff_VEF_Face_PEMFC_Nafion::remplir_nu(DoubleTab& nu) const
{
  Op_Diff_VEF_base::remplir_nu(nu);
  const DoubleTab& por = por_naf_.valeurs();
  const DoubleTab& eps = eps_naf_.valeurs();

  if(sub_type(Champ_Uniforme, por_naf_))
    {
      nu /= 1.-por(0,0);
    }
  else
    {
      DoubleTab subpor(por);
      subpor *= -1.;
      subpor += 1.;		// 1-por_naf
      tab_divide_any_shape(nu,subpor);
    }

  if(sub_type(Champ_Uniforme, eps_naf_))
    {
      nu /= eps(0,0);
    }
  else
    {
      tab_divide_any_shape(nu,eps);
    }

  Cerr << "Op_Diff_VEF_Face_PEMFC_Nafion::remplir_nu" << finl;
  Cerr << "nu min max " << mp_min_vect(nu) << " " << mp_max_vect(nu) << finl;
}

void Op_Diff_VEF_Face_PEMFC_Nafion::mettre_a_jour(double temps)
{
  remplir_nu(nu_);
}
