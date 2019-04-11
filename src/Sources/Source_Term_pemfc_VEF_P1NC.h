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
// File      : Source_Term_pemfc_VEF_P1NC.h
// Directory : $PEMFC_ROOT/src/Sources
//
/////////////////////////////////////////////////////////////////////////////

#ifndef Source_Term_pemfc_VEF_P1NC_included
#define Source_Term_pemfc_VEF_P1NC_included

#include <Source_Term_pemfc_base.h>

/////////////////////////////////////////////////////////////////////////////
//
// .DESCRIPTION : class Source_Term_pemfc_VEF_P1NC
//
// <Description of class Source_Term_pemfc_VEF_P1NC>
//
/////////////////////////////////////////////////////////////////////////////

class Source_Term_pemfc_VEF_P1NC : public Source_Term_pemfc_base
{

  Declare_instanciable( Source_Term_pemfc_VEF_P1NC ) ;

public :

protected :
  void associer_zones(const Zone_dis& ,const Zone_Cl_dis& );
  void associer_pb(const Probleme_base& );
  void remplir_volumes();

  REF(Zone_VEF) la_zone_VEF;
  REF(Zone_Cl_VEF) la_zcl_VEF;
};

#endif /* Source_Term_pemfc_VEF_P1NC_included */
