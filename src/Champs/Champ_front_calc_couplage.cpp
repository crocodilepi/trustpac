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
// File      : Champ_front_calc_couplage.cpp
// Directory : $PEMFC_ROOT/src/Champs
//
/////////////////////////////////////////////////////////////////////////////

#include <Champ_front_calc_couplage.h>
#include <Loi_Fermeture_base.h>
#include <Interprete.h>
#include <Champ_Inc_base.h>
#include <Frontiere_dis_base.h>
#include <Frontiere.h>

Implemente_instanciable( Champ_front_calc_couplage, "Champ_front_calc_couplage", Champ_front_calc ) ;

Sortie& Champ_front_calc_couplage::printOn( Sortie& os ) const
{
  Champ_front_calc::printOn( os );
  return os;
}

Entree& Champ_front_calc_couplage::readOn( Entree& is )
{
  Champ_front_calc::readOn( is );
  //is >> nom_loi_;
  //loi_ = ref_cast( Loi_Fermeture_base, interprete( ).objet( nom_loi_ ) );
  return is;
}

void Champ_front_calc_couplage::mettre_a_jour(double temps)
{
  Cerr << "Champ_front_calc_couplage::mettre_a_jour temps=" << temps << finl;
  Cerr << "Debug Champ_Inc temps_present= " << l_inconnue.valeur().temps() << ", champ_front nom=" << l_inconnue.valeur().le_nom() << ", frontiere nom=" << front_dis().frontiere().le_nom() << ", temps_defaut=" << get_temps_defaut() << ", les_valeurs=" << les_valeurs <<", valeurs="<< l_inconnue.valeur().valeurs(temps) << finl;
  //loi_.valeur().mettre_a_jour(temps);
  Champ_front_calc::mettre_a_jour(temps);
  Cerr << "Champ_front apres MAJ, au temps=" << temps << ", valeurs=" << valeurs_au_temps(temps) << finl;
}

int Champ_front_calc_couplage::initialiser(double temps,const Champ_Inc_base& inco)
{
  int initialized = 0;
  if(l_inconnue.non_nul())
    {
      initialized = Champ_front_calc::initialiser(temps, l_inconnue.valeur());
      Cerr << "Champ_front_calc_couplage::initialiser with inco associe=" << l_inconnue.valeur() << finl;
    }
  else
    {
      initialized = Champ_front_calc::initialiser(temps, inco);
      Cerr << "Champ_front_calc_couplage::initialiser with inco de l'equation=" << inco << finl;
    }

  // debug
  //Cerr<< "Valeur apres initialise: " << valeurs() << finl;

  if(!initialized)
    return 0;
  return 1;
}

DoubleTab& Champ_front_calc_couplage::valeurs_au_temps(double temps)
{
  // DEBUG
//  int nbcases = les_valeurs->nb_cases();
//  Cerr << "Debug Champ_front_calc_couplage::valeurs_au_temps nom=" << l_inconnue.valeur().le_nom() << ", frontiere nom=" << front_dis().frontiere().le_nom() << "requested time=" << temps << ", availables times="<< nbcases << ", times available=";
//  for(int i=0; i<nbcases; i++)
//    Cerr << "  " << les_valeurs[i].temps();
//  Cerr <<", time default="<< get_temps_defaut() << finl;
//  for(int i=0; i<nbcases; i++)
//    {
//      Cerr << "Debug case=" << i << ", temps=" << les_valeurs[i].temps() << ", valeurs=" <<  les_valeurs[i].valeurs() << finl;
//    }
  // END DEBUG
  for(int i=0; i<les_valeurs->nb_cases(); i++)
    {
      if(est_egal(les_valeurs[i].temps(),temps))
        return les_valeurs[i].valeurs();
    }

  return valeurs(); // For compilers
}

// Description:
//    Renvoie les valeurs au temps desire.
//    Sinon, sort en erreur.
const DoubleTab& Champ_front_calc_couplage::valeurs_au_temps(double temps) const
{
  // DEBUG
//  int nbcases = les_valeurs->nb_cases();
//  Cerr << "Debug Champ_front_calc_couplage::valeurs_au_temps nom=" << l_inconnue.valeur().le_nom() << ", frontiere nom=" << front_dis().frontiere().le_nom() << "requested time=" << temps << ", availables times="<< nbcases << ", times available=";
//  for(int i=0; i<nbcases; i++)
//    Cerr << "  " << les_valeurs[i].temps();
//  Cerr <<", time default="<< get_temps_defaut() << finl;
//  for(int i=0; i<nbcases; i++)
//    {
//      Cerr << "Debug case=" << i << ", temps=" << les_valeurs[i].temps() << ", valeurs=" <<  les_valeurs[i].valeurs() << finl;
//    }
  // END DEBUG
  for(int i=0; i<les_valeurs->nb_cases(); i++)
    {
      if(est_egal(les_valeurs[i].temps(),temps))
        return les_valeurs[i].valeurs();
    }

  return valeurs(); // For compilers
}
