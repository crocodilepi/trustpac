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
// File:        Convection_Concentration_Gaz.cpp
// Directory:   $TRUST_ROOT/src/ThHyd
// Version:     /main/17
//
//////////////////////////////////////////////////////////////////////////////

#include <Convection_Concentration_Gaz.h>
#include <Probleme_base.h>
#include <Milieu_base.h>
#include <Discretisation_base.h>
#include <Param.h>
#include <Fluide_Quasi_Compressible.h>
#include <Loi_Etat_Melange_GP_Fraction_Molaire.h>

Implemente_instanciable(Convection_Concentration_Gaz,"Convection_Concentration_Gaz",Equation_base);
// XD convection_concentration_gaz eqn_base convection_concentration_gaz -1 equation of cg
// XD attr diffusion suppress_param diffusion 1 no diffusion

// Description:
//    Simple appel a Equation_base::printOn(Sortie&)
// Precondition:
// Parametre: Sortie& is
//    Signification: un flot de sortie
//    Valeurs par defaut:
//    Contraintes:
//    Acces: entree/sortie
// Retour: Sortie&
//    Signification: le flot de sortie modifie
//    Contraintes:
// Exception:
// Effets de bord:
// Postcondition: la methode ne modifie pas l'objet
Sortie& Convection_Concentration_Gaz::printOn(Sortie& is) const
{
  return Equation_base::printOn(is);
}


// Description:
//    cf Equation_base::readOn(Entree&)
// Precondition:
// Parametre: Entree& is
//    Signification: un flot d'entree
//    Valeurs par defaut:
//    Contraintes:
//    Acces: entree/sortie
// Retour: Entree&
//    Signification: le flot d'entree modifie
//    Contraintes:
// Exception: terme convectif non specifie
// Effets de bord:
// Postcondition:
Entree& Convection_Concentration_Gaz::readOn(Entree& is)
{
  Equation_base::readOn(is);
  Nom nom="Convection_";
  nom+=inconnue().le_nom(); // On ajoute le nom de l'inconnue pour prevoir une equation de scalaires passifs
  terme_convectif.set_fichier(nom);
  terme_convectif.set_description((Nom)"Convective mass transfer rate=Integral(-C*u*ndS)[m"+(Nom)(dimension+bidim_axi)+".Mol.s-1]");
  return is;
}

void Convection_Concentration_Gaz::set_param(Param& param)
{
  Equation_base::set_param(param);
  param.ajouter_non_std("convection",(this));
  param.ajouter_condition("is_read_convection","The convection operator must be read, select negligeable type if you want to neglect it.");
}


int Convection_Concentration_Gaz::lire_motcle_non_standard(const Motcle& mot, Entree& is)
{
  if (mot=="convection")
    {
      Cerr << "Reading and typing of the convection operator : " << finl;
      //   const Champ_base& ch_vitesse_transportante = vitesse_pour_transport();

      const Probleme_base& pb = probleme();
      const Champ_base& ch_vitesse_transportante =pb.get_champ("Um");
      //const Champ_base& vit_transportante =pb.get_champ("vitesse");


      associer_vitesse(ch_vitesse_transportante);


      terme_convectif.associer_vitesse(ch_vitesse_transportante);
      is >> terme_convectif;
      return 1;
    }
  else
    return Equation_base::lire_motcle_non_standard(mot,is);
  return 1;
}
// Description:
//    Renvoie le nombre d'operateurs de l'equation:
//    2 pour une equation de diffusion.
// Precondition:
// Parametre:
//    Signification:
//    Valeurs par defaut:
//    Contraintes:
//    Acces:
// Retour: int
//    Signification: le nombre d'operateurs de l'equation
//    Contraintes: toujours egal a 2
// Exception:
// Effets de bord:
// Postcondition: la methode ne modifie pas l'objet
int Convection_Concentration_Gaz::nombre_d_operateurs() const
{
  return 1;
}

// Description:
//    Renvoie l'operateur specifie par son index:
//     renvoie terme_convectif si i = 1
//     exit si i>1
//    (version const)
// Precondition:
// Parametre: int i
//    Signification: l'index de l'operateur a renvoyer
//    Valeurs par defaut:
//    Contraintes: 0 <= i <= 1
//    Acces: entree
// Retour: Operateur&
//    Signification: l'operateur specifie
//    Contraintes: reference constante
// Exception: l'equation n'a pas plus de 2 operateurs
// Effets de bord:
// Postcondition: la methode ne modifie pas l'objet
const Operateur& Convection_Concentration_Gaz::operateur(int i) const
{
  switch(i)
    {
    case 0:
      return terme_convectif;
    default :
      Cerr << "Error for Convection_Concentration_Gaz::operateur(int i)" << finl;
      Cerr << "Convection_Concentration_Gaz has " << nombre_d_operateurs() <<" operators "<<finl;
      Cerr << "and you are trying to access the " << i <<" th one."<< finl;
      exit();
    }
  throw;
}

// Description:
//    Renvoie l'operateur specifie par son index:
//     renvoie terme_diffusif si i = 0
//     renvoie terme_convectif si i = 1
//     exit si i>1
// Precondition:
// Parametre: int i
//    Signification: l'index de l'operateur a renvoyer
//    Valeurs par defaut:
//    Contraintes: 0 <= i <= 1
//    Acces: entree
// Retour: Operateur&
//    Signification: l'operateur specifie
//    Contraintes:
// Exception: l'equation n'a pas plus de 2 operateurs
// Effets de bord:
// Postcondition:
Operateur& Convection_Concentration_Gaz::operateur(int i)
{
  switch(i)
    {
    case 0:
      return terme_convectif;
    default :
      Cerr << "Error for Convection_Concentration_Gaz::operateur(int i)" << finl;
      Cerr << "Convection_Concentration_Gaz has " << nombre_d_operateurs() <<" operators "<<finl;
      Cerr << "and you are trying to access the " << i <<" th one."<< finl;
      exit();
    }
  // Pour les compilos!!
  throw ;
}



void  Convection_Concentration_Gaz::discretiser()
{
  discretisation().discretiser_champ("Temperature",zone_dis().valeur(),"cg","1", 1, schema_temps().nb_valeurs_temporelles() ,0.,la_concentration_);

  champs_compris_.ajoute_champ(la_concentration_);
  Equation_base::discretiser();
  const Fluide_Quasi_Compressible& un_fluideQC = ref_cast(Fluide_Quasi_Compressible,le_fluide_.valeur());
  Loi_Etat_Melange_GP_Fraction_Molaire& loi_etat = ref_cast_non_const(Loi_Etat_Melange_GP_Fraction_Molaire,un_fluideQC.loi_etat().valeur());
  loi_etat.associer_cg(inconnue());
}


const Champ_base& Convection_Concentration_Gaz::vitesse_pour_transport()
{
  return probleme().equation(0).inconnue();
}

void Convection_Concentration_Gaz::associer_milieu_base(const Milieu_base& un_milieu)
{
  const Fluide_Quasi_Compressible& un_fluideQC = ref_cast(Fluide_Quasi_Compressible,un_milieu);
  le_fluide_=(un_fluideQC);
}
const Milieu_base& Convection_Concentration_Gaz::milieu() const
{
  return le_fluide_.valeur();
}
Milieu_base& Convection_Concentration_Gaz::milieu()
{
  return le_fluide_.valeur();
}
