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
// File      : Source_Term_pemfc_base.h
// Directory : $PEMFC_ROOT/src/Sources
//
/////////////////////////////////////////////////////////////////////////////

#ifndef Source_Term_pemfc_base_included
#define Source_Term_pemfc_base_included

#include <Source_base.h>
#include <Ref_Champ_base.h>
#include <Probleme_base.h>
#include <Matrice_Morse.h>
#include <DoubleTab.h>
#include <Zone_dis.h>
#include <Zone_Cl_dis.h>

/////////////////////////////////////////////////////////////////////////////
//
// .DESCRIPTION : class Source_Term_pemfc_base
//  Modele de source de dissolution des especes (H2, N2, O2) de type:
//  Sa = Da * gamma_CL / e_ionomer * (Ceq - C)
//  avec
//  Ceq = c * R * T * H
// <Description of class Source_Term_pemfc_base>
//
/////////////////////////////////////////////////////////////////////////////

class Source_Term_pemfc_base : public Source_base
{

  Declare_base( Source_Term_pemfc_base ) ;

public :
  virtual DoubleTab& ajouter(DoubleTab& ) const;
  virtual DoubleTab& calculer(DoubleTab& ) const;
  virtual void mettre_a_jour(double temps);
  virtual void contribuer_a_avec(const DoubleTab&, Matrice_Morse&) const;
  //virtual void associer_zones(const Zone_dis& ,const Zone_Cl_dis& )=0;
  //virtual void associer_pb(const Probleme_base& )=0;
protected :
  Nom nom_espece_;						// nom d'espece dans la liste { H2 O2 N2 vap H2O}
  REF(Champ_base)  T_;   				// champ_inc Temperature du probleme couple T
  REF(Champ_base) Da_;					// champ conductivite Da
  REF(Champ_base)  C_;					// Champ_Inc dissolved concentration
  REF(Champ_base)  c_;					// Champ_Inc concentration du problem couple c
  double epsilon_ionomer_;				// ionomer proportion
  double epsilon_	; 					// porosity
  double gamma_CL_ ; 					// specific surface m2/m3
  double thickness_ionomer_;			// = (1-epsilon)*epsilon_ionomer / gamma_CL

  DoubleVect volumes_;
  virtual void remplir_volumes()=0;
  double eval_f(double diffu, double Ci, double ci, double T) const;
  double eval_derivee_f(double diffu) const;
};

#endif /* Source_Term_pemfc_base_included */
