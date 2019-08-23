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
// File      : Loi_Fermeture_PEMFC_Cas2.cpp
// Directory : $PEMFC_ROOT/src
//
/////////////////////////////////////////////////////////////////////////////

#include <Loi_Fermeture_PEMFC_Cas2.h>
#include <Probleme_base.h>
#include <Discretisation_base.h>
#include <Equation_base.h>
#include <Param.h>
#include <Interprete.h>

Implemente_instanciable( Loi_Fermeture_PEMFC_Cas2, "Loi_Fermeture_PEMFC_Cas2", Loi_Fermeture_PEMFC_base ) ;
// XD loi_fermeture_pemfc_cas2 loi_fermeture_base loi_fermeture_pemfc_cas2 -1 Loi for test only

Sortie& Loi_Fermeture_PEMFC_Cas2::printOn( Sortie& os ) const
{
  Loi_Fermeture_PEMFC_base::printOn( os );
  return os;
}

Entree& Loi_Fermeture_PEMFC_Cas2::readOn( Entree& is )
{
  Loi_Fermeture_PEMFC_base::readOn( is );

  return is;
}

void Loi_Fermeture_PEMFC_Cas2::discretiser(const Discretisation_base& dis)
{
  Loi_Fermeture_PEMFC_base::discretiser(dis);
  ref_equation_=mon_probleme().get_equation_by_name("Convection_diffusion_Concentration");

  int nb_comp = equation().inconnue().valeur().nb_comp();

  dis.discretiser_champ("champ_elem",equation().zone_dis().valeur(),"pemfc_cas2","unit", nb_comp*nb_comp, 0., diffu_);
  champs_compris_.ajoute_champ(diffu_);
  diffu_.valeur().fixer_nature_du_champ(multi_scalaire);

  dis.discretiser_champ("champ_elem",equation().zone_dis().valeur(),"Ni","mol/m2/s",nb_comp*dimension,0.,Ni_);
  Ni_->fixer_nature_du_champ(multi_scalaire);
  champs_compris_.ajoute_champ(Ni_);

  dis.discretiser_champ("champ_elem",equation().zone_dis().valeur(),"ud","m/s", dimension,0.,ud_);
  ud_->fixer_nature_du_champ(vectoriel);
  champs_compris_.ajoute_champ(ud_);

  // vitesse (GDL)
  dis.discretiser_champ("vitesse",equation().zone_dis().valeur(),"ug","m/s", dimension,1 /* une case en temps */,0.,ug_gdl_);
  ug_gdl_->fixer_nature_du_champ(vectoriel);
  champs_compris_.ajoute_champ(ug_gdl_);
  ug_gdl_.associer_eqn(equation());

  dis.discretiser_champ("temperature",equation().zone_dis().valeur(),"cg","mol/m3",1,1 /* une case en temps */,0.,cg_gdl_);
  cg_gdl_->fixer_nature_du_champ(scalaire);
  champs_compris_.ajoute_champ(cg_gdl_);
  cg_gdl_.associer_eqn(equation());

  // fraction molaire (GDL)
  dis.discretiser_champ("temperature",equation().zone_dis().valeur(),"Xi","sans_dimension", nb_comp,1 /* une case en temps */,0.,Xi_gdl_);
  Xi_gdl_->fixer_nature_du_champ(multi_scalaire);
  champs_compris_.ajoute_champ(Xi_gdl_);
  Xi_gdl_.associer_eqn(equation());
  //Cerr << "DEBUG valeur Xi_GDL apres discretisation " << Xi_gdl_.valeurs() << finl;

  // Flux surfacique de concentration sortant (GDL)
  dis.discretiser_champ("temperature",equation().zone_dis().valeur(),"Ns","mol/s", nb_comp,1 /* une case en temps */,0.,Ns_gdl_);
  Ns_gdl_->fixer_nature_du_champ(multi_scalaire);
  champs_compris_.ajoute_champ(Ns_gdl_);
  Ns_gdl_.associer_eqn(equation());
}

void Loi_Fermeture_PEMFC_Cas2::set_param(Param& param)
{
  Loi_Fermeture_PEMFC_base::set_param(param);
  param.ajouter("Nom_pb_T",&nom_pb_T_, Param::OPTIONAL); // XD_ADD_P string pb_T nom du probleme de conduction thermique (inconnu Temperature)
  param.ajouter("por",&ch_por_, Param::OPTIONAL);		// XD_ADD_P champ_base champ de porosite
  param.ajouter("tor",&ch_tor_, Param::OPTIONAL);		// XD_ADD_P champ_base champ de tortuosite
  param.ajouter("Rp",&ch_Rp_, Param::OPTIONAL);			// XD_ADD_P champ_base pore radius [m]
  param.ajouter("K",&ch_K_, Param::OPTIONAL);			// XD_ADD_P champ_base permeabilite [m2]
  param.ajouter("masses_molaires",&Mi_,Param::REQUIRED);  // XD_ADD_P list masses molaires [kg/mol]
}

void Loi_Fermeture_PEMFC_Cas2::completer()
{
  Loi_Fermeture_PEMFC_base::completer();
}

void Loi_Fermeture_PEMFC_Cas2::mettre_a_jour(double temps)
{
  Loi_Fermeture_PEMFC_base::mettre_a_jour(temps);
}
