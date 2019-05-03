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
// File      : Op_Diff_VEF_Face_transport_ie.h
// Directory : $PEMFC_ROOT/src/Operateurs
//
/////////////////////////////////////////////////////////////////////////////

#ifndef Op_Diff_VEF_Face_transport_ie_included
#define Op_Diff_VEF_Face_transport_ie_included

#include <Op_Diff_VEF_Face.h>
#include <Param.h>
#include <Champ_Don.h>

/////////////////////////////////////////////////////////////////////////////
//
// .DESCRIPTION : class Op_Diff_VEF_Face_transport_ie
//
// <Description of class Op_Diff_VEF_Face_transport_ie>
//
/////////////////////////////////////////////////////////////////////////////

class Op_Diff_VEF_Face_transport_ie : public Op_Diff_VEF_Face
{

  Declare_instanciable( Op_Diff_VEF_Face_transport_ie ) ;

public :
  void set_param(Param& param);
  void completer();
  void remplir_nu(DoubleTab&) const;
  void mettre_a_jour(double temps);
protected :
  Nom diffu_name_;
  Champ_Don Cdl_;
};

#endif /* Op_Diff_VEF_Face_transport_ie_included */
