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
// File      : Loi_Fermeture_PEMFC_Cas4.cpp
// Directory : $PEMFC_ROOT/src
//
/////////////////////////////////////////////////////////////////////////////

#include <Loi_Fermeture_PEMFC_Cas4.h>
#include <Probleme_base.h>
#include <Discretisation_base.h>
#include <Equation_base.h>
#include <Schema_Temps_base.h>
#include <Param.h>

Implemente_instanciable_sans_constructeur( Loi_Fermeture_PEMFC_Cas4, "Loi_Fermeture_PEMFC_Cas4", Loi_Fermeture_PEMFC_base ) ;
// XD loi_fermeture_pemfc_cas4 loi_fermeture_base loi_fermeture_pemfc_cas4 -1 Loi for test only

Loi_Fermeture_PEMFC_Cas4::Loi_Fermeture_PEMFC_Cas4()
{
  is_cas4_=(1);
}
Sortie& Loi_Fermeture_PEMFC_Cas4::printOn( Sortie& os ) const
{
  Loi_Fermeture_PEMFC_base::printOn( os );
  return os;
}

Entree& Loi_Fermeture_PEMFC_Cas4::readOn( Entree& is )
{
  Loi_Fermeture_PEMFC_base::readOn( is );

  is_cas4_=1;
#if 0
  if (VGDL_.size_array()!=dimension)
    {
      Cerr<<" VGDL faux "<<finl;
      exit();
    }
#endif
  return is;
}

void Loi_Fermeture_PEMFC_Cas4::preparer_calcul()
{
  Loi_Fermeture_PEMFC_base::preparer_calcul();
  if (mon_probleme().reprise_effectuee())
    {
      Equation_base& NS=ref_equation_.valeur().probleme().equation(0);
      double temps=NS.schema_temps().temps_courant();

      NS.zone_Cl_dis().imposer_cond_lim(NS.inconnue(),temps);
      NS.inconnue().valeurs().echange_espace_virtuel();
    }
}

void Loi_Fermeture_PEMFC_Cas4::set_param(Param& param)
{
  Loi_Fermeture_PEMFC_base::set_param(param);
  //param.ajouter("VGDL",&VGDL_); // XD_ADD_P list  VGDL pour verification
  //param.ajouter("eps_sur_tau2",&eps_sur_tau2_); // XD_ADD_P field_base eps/tau/tau
  param.ajouter("por",&ch_por_, Param::OPTIONAL);		// XD_ADD_P champ_base champ de porosite
  param.ajouter("tor",&ch_tor_, Param::OPTIONAL);		// XD_ADD_P champ_base champ de tortuosite
}


void Loi_Fermeture_PEMFC_Cas4::discretiser(const Discretisation_base& dis)
{
  Loi_Fermeture_PEMFC_base::discretiser(dis);
  ref_equation_=mon_probleme().get_equation_by_name("Convection_Diffusion_fraction_molaire_QC");

  int nb_comp = equation().inconnue().valeur().nb_comp();

  dis.discretiser_champ("champ_elem",equation().zone_dis().valeur(),"pemfc_cas4","unit",nb_comp*nb_comp,0.,diffu_);
  champs_compris_.ajoute_champ(diffu_);

  dis.discretiser_champ("champ_elem",equation().zone_dis().valeur(),"Ni","unit",nb_comp*dimension,0.,Ni_);
  Ni_->fixer_nature_du_champ(multi_scalaire);
  champs_compris_.ajoute_champ(Ni_);

  dis.discretiser_champ("champ_elem",equation().zone_dis().valeur(),"ud","m/s", dimension,0.,ud_);
  ud_->fixer_nature_du_champ(vectoriel);
  champs_compris_.ajoute_champ(ud_);

  dis.discretiser_champ("vitesse",equation().zone_dis().valeur(),"UM","m/s", dimension,1 /* une case en temps */,0.,Um_);
  Um_->fixer_nature_du_champ(vectoriel);
  champs_compris_.ajoute_champ(Um_);

  dis.discretiser_champ("temperature",equation().zone_dis().valeur(),"ci","mol/m3", nb_comp,1 /* une case en temps */,0.,ci_ch_);
  ci_ch_->fixer_nature_du_champ(multi_scalaire);
  champs_compris_.ajoute_champ(ci_ch_);
  ci_ch_.associer_eqn(equation());

  dis.discretiser_champ("temperature",equation().zone_dis().valeur(),"Ns","mol/s", nb_comp,1 /* une case en temps */,0.,Ns_ch_);
  Ns_ch_->fixer_nature_du_champ(multi_scalaire);
  champs_compris_.ajoute_champ(Ns_ch_);
  Ns_ch_.associer_eqn(equation());
}
