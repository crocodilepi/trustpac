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
// File      : Loi_Fermeture_PEMFC_base.h
// Directory : $PEMFC_ROOT/src
//
/////////////////////////////////////////////////////////////////////////////

#ifndef Loi_Fermeture_PEMFC_base_included
#define Loi_Fermeture_PEMFC_base_included

#include <Loi_Fermeture_base.h>
#include <Champ_Fonc.h>
#include <Champ_Don.h>
#include <Champ_Inc.h>
#include <Ref_Equation_base.h>

/////////////////////////////////////////////////////////////////////////////
//
// .DESCRIPTION : class Loi_Fermeture_PEMFC_base
//
// <Description of class Loi_Fermeture_PEMFC_base>
//
/////////////////////////////////////////////////////////////////////////////

class Loi_Fermeture_PEMFC_base : public Loi_Fermeture_base
{

  Declare_base( Loi_Fermeture_PEMFC_base ) ;

public :

  void discretiser(const Discretisation_base& );
  void mettre_a_jour(double);
  // void set_param(Param& param);
  // void completer();
  void calculer_Ni_ud(const double& temps) ;
  const Equation_base& equation() const
  {
    return ref_equation_.valeur();
  } ;
protected :
  int is_cas4_;
  Champ_Fonc Ni_,ud_,diffu_;
  Champ_Don eps_sur_tau2_;
  Champ_Inc Um_;
  REF(Equation_base) ref_equation_;

  ArrOfDouble VGDL_;

};

#endif /* Loi_Fermeture_PEMFC_base_included */
