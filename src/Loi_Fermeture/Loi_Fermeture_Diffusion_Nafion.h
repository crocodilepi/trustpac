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
// File      : Loi_Fermeture_Diffusion_Nafion.h
// Directory : $PEMFC_ROOT/src/Loi_Fermeture
//
/////////////////////////////////////////////////////////////////////////////

#ifndef Loi_Fermeture_Diffusion_Nafion_included
#define Loi_Fermeture_Diffusion_Nafion_included

#include <Loi_Fermeture_base.h>
#include <Ref_Champ_base.h>
#include <Champ_Fonc.h>
#include <Ref_Equation_base.h>
#include <DoubleTab.h>
#include <Param.h>

/////////////////////////////////////////////////////////////////////////////
//
// .DESCRIPTION : class Loi_Fermeture_Diffusion_Nafion
//
// <Description of class Loi_Fermeture_Diffusion_Nafion>
//
/////////////////////////////////////////////////////////////////////////////

class Loi_Fermeture_Diffusion_Nafion : public Loi_Fermeture_base
{

  Declare_instanciable( Loi_Fermeture_Diffusion_Nafion ) ;

public :
  // associer avec le probleme de diffusion
  void associer_pb_base(const Probleme_base&);
  // discretiser le champ diffusivite
  void discretiser(const Discretisation_base& );
  // recuperer tous Reference
  void completer();
  // mettre a jour une fois avant du calcul iteratif
  void preparer_calcul();
  // mettre a jour le tableau de valeurs
  void mettre_a_jour(double temps);
  //
  void set_param(Param& param);
protected :
  Champ_Fonc diffu_;						// champ P0
  const Equation_base& equation() const
  {
    return ref_equation_.valeur();
  } ;
  REF(Equation_base) ref_equation_;

  // Motcles nom_especes_compris_ = {"H2", "O2", "H2O", "N2", "vap"};

  Nom nom_espece_, nom_champ_T_, nom_champ_C_;	// pour readOn
  double default_value_T_ = 353.15;				// dans le cas T constant
  double CSO3_ = 2036.;							// au cas ou "H2O" ou "vap"

  //REF(Champ_base) T_;			// champ reference pour temperature (optionnel) -> champ P0 avec meme support que diffu_
  //REF(Champ_base) C_;			// champ reference pour concentration dissolue (optionnel) -> champ P0 avec meme support que diffu_
  DoubleTab T_, C_;				// tableau des valeurs du champ T et C (P0)
  double eval(double T, double C);
};

#endif /* Loi_Fermeture_Diffusion_Nafion_included */
