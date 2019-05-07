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
// File      : Op_Diff_VEF_Face_PEMFC_Nafion.h
// Directory : $PEMFC_ROOT/src/Operateurs
//
/////////////////////////////////////////////////////////////////////////////

#ifndef Op_Diff_VEF_Face_PEMFC_Nafion_included
#define Op_Diff_VEF_Face_PEMFC_Nafion_included

#include <Op_Diff_VEF_Face.h>
#include <Param.h>
#include <Champ_Don.h>

/////////////////////////////////////////////////////////////////////////////
//
// .DESCRIPTION : class Op_Diff_VEF_Face_PEMFC_Nafion
//  Cette classe a pour but de calculer l'operateur de diffusion de type: div(D.grad(C) avec
//  un champ de coefficient de diffusion particulier D = D_eff / ((1-por_naf)eps_naf)
//  D_eff champ de coefficient de diffusion (couple avec un loi_fermeture_diffusion_nafion via le champ D_i_eff)
//  por_naf porosite de Nafion
//  eps_naf ionomer proportiton de Nafion
// <Description of class Op_Diff_VEF_Face_PEMFC_Nafion>
//
/////////////////////////////////////////////////////////////////////////////

class Op_Diff_VEF_Face_PEMFC_Nafion : public Op_Diff_VEF_Face
{

  Declare_instanciable( Op_Diff_VEF_Face_PEMFC_Nafion ) ;

public :
  void set_param(Param& param);
  void completer();
  void remplir_nu(DoubleTab&) const;
  void mettre_a_jour(double);
protected :
  Nom diffu_name_;
  Champ_Don por_naf_;			// porosite de Nafion
  Champ_Don eps_naf_;			// ionomer proportionnel de Nafion
};

#endif /* Op_Diff_VEF_Face_PEMFC_Nafion_included */
