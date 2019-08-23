/****************************************************************************
* Copyright (c) 2019, CEA
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
// File      : Champ_front_calc_couplage.h
// Directory : $PEMFC_ROOT/src/Champs
//
/////////////////////////////////////////////////////////////////////////////

#ifndef Champ_front_calc_couplage_included
#define Champ_front_calc_couplage_included

#include <Champ_front_calc.h>
#include <Ref_Loi_Fermeture_base.h>

/////////////////////////////////////////////////////////////////////////////
//
// .DESCRIPTION : class Champ_front_calc_couplage
//
// <Description of class Champ_front_calc_couplage>
//
/////////////////////////////////////////////////////////////////////////////

class Champ_front_calc_couplage : public Champ_front_calc
{

  Declare_instanciable( Champ_front_calc_couplage ) ;

public :
  void mettre_a_jour(double temps);
  int initialiser(double, const Champ_Inc_base&);
  DoubleTab& valeurs_au_temps(double temps);
  const DoubleTab& valeurs_au_temps(double temps) const;
protected :

  Nom nom_loi_;
  REF(Loi_Fermeture_base) loi_;


};

#endif /* Champ_front_calc_couplage_included */
