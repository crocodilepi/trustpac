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
//////////////////////////////////////////////////////////////////////////////
//
// File:        Convection_Concentration_Gaz.h
// Directory:   $TRUST_ROOT/src/ThHyd
// Version:     /main/19
//
//////////////////////////////////////////////////////////////////////////////


#ifndef Convection_Concentration_Gaz_included
#define Convection_Concentration_Gaz_included

#include <Equation_base.h>
#include <Operateur_Conv.h>
#include <Schema_Temps.h>
#include <Ref_Champ_base.h>
class Fluide_Incompressible;
#include <Ref_Fluide_Quasi_Compressible.h>

class Champ_Don;

//////////////////////////////////////////////////////////////////////////////
//
// .DESCRIPTION
//     classe Convection_Concentration_Gaz
//     Cette classe est la base des equations modelisant le transport
//     d'un scalaire.
//     Cette classe porte les termes communs a l'equation
//     de transport d'un scalaire en regime laminaire.
//     On se place sous l'hypothese de fluide incompressible.
//         DT/dt = div(terme visqueux) + termes sources/rho_0
//     avec DT/dt : derivee particulaire du scalaire
//     Rq: l'implementation de la classe permet bien sur de negliger
//         certains termes de l'equation (le terme diffusif,
//         le terme convectif,tel ou tel terme source).
// .SECTION voir aussi
//     Equation_base
//     Classe abstraite
//     Methodes abstraites:
//       Entree& lire(const Motcle&, Entree&)
//       const Champ_Inc& inconnue() const
//       Champ_Inc& inconnue()
//////////////////////////////////////////////////////////////////////////////
class Convection_Concentration_Gaz : public Equation_base
{
  Declare_instanciable(Convection_Concentration_Gaz);

public :

  void set_param(Param& titi);
  int lire_motcle_non_standard(const Motcle&, Entree&);
  int nombre_d_operateurs() const;
  const Operateur& operateur(int) const;
  Operateur& operateur(int);
  inline void associer_vitesse(const Champ_base& );
  inline const Champ_base& vitesse_transportante() const;
  inline const Champ_Inc& inconnue() const ;
  inline Champ_Inc& inconnue() ;
  virtual const Champ_base& vitesse_pour_transport();

  void discretiser();


  const Milieu_base& milieu() const;
  Milieu_base& milieu();
  void associer_milieu_base(const Milieu_base& );

protected :

  REF(Champ_base) la_vitesse_transportante;
  Operateur_Conv terme_convectif;
  Champ_Inc la_concentration_;
  REF(Fluide_Quasi_Compressible) le_fluide_;
};


// Description:
//    Renvoie une reference sur le champ representant
//    la vitesse transportante.
// Precondition:
// Parametre:
//    Signification:
//    Valeurs par defaut:
//    Contraintes:
//    Acces:
// Retour: Champ_Inc&
//    Signification: le champ representant la vitesse transportante
//    Contraintes: reference constante
// Exception:
// Effets de bord:
// Postcondition:
inline const Champ_base& Convection_Concentration_Gaz::vitesse_transportante() const
{
  return la_vitesse_transportante.valeur();
}


// Description:
//    Associe la vitesse transportante a l'equation.
// Precondition:
// Parametre: Champ_Inc& vit
//    Signification: le champ a affecter a la vitesse transportante
//    Valeurs par defaut:
//    Contraintes: reference constante
//    Acces: entree
// Retour:
//    Signification:
//    Contraintes:
// Exception:
// Effets de bord:
// Postcondition: une vitesse transportante a ete associe a l'equation
inline void Convection_Concentration_Gaz::associer_vitesse(const Champ_base& vit)
{
  la_vitesse_transportante = ref_cast(Champ_Inc_base,vit);
}


inline const Champ_Inc& Convection_Concentration_Gaz::inconnue() const
{
  return la_concentration_;
}
inline Champ_Inc& Convection_Concentration_Gaz::inconnue()
{
  return la_concentration_;
}


#endif
