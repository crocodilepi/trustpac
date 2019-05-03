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
// File      : Loi_Fermeture_transport_ionique.cpp
// Directory : $PEMFC_ROOT/src/Loi_Fermeture
//
/////////////////////////////////////////////////////////////////////////////

#include <Loi_Fermeture_transport_ionique.h>
#include <Champ_base.h>
#include <Champ_P1NC.h>
#include <Interprete.h>
#include <Probleme_base.h>
#include <Operateur.h>
#include <Operateur_base.h>
#include <Zone_VF.h>
#include <Domaine.h>
#include <Milieu_base.h>

Implemente_instanciable( Loi_Fermeture_transport_ionique, "Loi_Fermeture_transport_ionique", Loi_Fermeture_base ) ;

Sortie& Loi_Fermeture_transport_ionique::printOn( Sortie& os ) const
{
  Loi_Fermeture_base::printOn( os );
  return os;
}

Entree& Loi_Fermeture_transport_ionique::readOn( Entree& is )
{
  Loi_Fermeture_base::readOn( is );

  const Domaine& dom = equation().probleme().domaine();
  dom.creer_tableau_elements(T_);
  dom.creer_tableau_elements(C_);

  return is;
}

void Loi_Fermeture_transport_ionique::set_param(Param& param)
{
  param.ajouter("T_0", &T_0_, Param::REQUIRED);
  param.ajouter("CSO3", &C_SO3_, Param::REQUIRED);
  param.ajouter("por_naf", &por_naf_, Param::REQUIRED);
  param.ajouter("eps_naf", &eps_naf_, Param::REQUIRED);
  param.ajouter("tor_naf", &tor_naf_, Param::REQUIRED);
  param.ajouter("nom_ssz_MB", &nom_ssz_MB_, Param::REQUIRED);
  param.ajouter("nom_ssz_CLc", &nom_ssz_CLc_, Param::REQUIRED);
  param.ajouter("nom_ssz_CLa", &nom_ssz_CLa_, Param::REQUIRED);
  param.ajouter("nom_pb_T",&nom_pb_T_,Param::REQUIRED);
  param.ajouter("nom_champ_T",&nom_champ_T_,Param::REQUIRED);
  param.ajouter("nom_pb_C_H2O",&nom_pb_C_,Param::REQUIRED);
  param.ajouter("nom_champ_C_H2O",&nom_champ_C_,Param::REQUIRED);
}


void Loi_Fermeture_transport_ionique::discretiser( const Discretisation_base& dis)
{
  Loi_Fermeture_base::discretiser(dis);
  ref_equation_=mon_probleme().get_equation_by_name("Conduction");
  dis.discretiser_champ("champ_elem",equation().zone_dis().valeur(),"kappa","unit", 1 ,0. , kappa_);
  champs_compris_.ajoute_champ(kappa_);

  dis.discretiser_champ("champ_elem",equation().zone_dis().valeur(),"I_i","unit", dimension ,0. , I_i_);
  I_i_ -> fixer_nature_du_champ(vectoriel);
  champs_compris_.ajoute_champ(I_i_);
}

inline void Loi_Fermeture_transport_ionique::completer()
{
  Loi_Fermeture_base::completer();
  // get the reference to the coupling fields
  Probleme_base& pb_T = ref_cast(Probleme_base,interprete().objet(nom_pb_T_));
  ch_T_ = pb_T.get_champ(nom_champ_T_);
  assert(ch_T_.valeur().que_suis_je().find("P1NC") != -1);
  Probleme_base& pb_phi = ref_cast(Probleme_base,interprete().objet(nom_pb_C_));
  ch_C_ = pb_phi.get_champ(nom_champ_C_);
  assert(ch_C_.valeur().que_suis_je().find("P1NC") != -1);
  ch_phi_ = equation().inconnue();
  MB_ = mon_probleme().domaine().ss_zone(nom_ssz_MB_);
  CL_a_ = mon_probleme().domaine().ss_zone(nom_ssz_CLa_);
  CL_c_ =  mon_probleme().domaine().ss_zone(nom_ssz_CLc_);
}

inline void Loi_Fermeture_transport_ionique::mettre_a_jour(double temps)
{
  // mettre a jour les champs couples -> interpoler vers P0
  const Zone_VF& la_zone = ref_cast(Zone_VF, equation().zone_dis().valeur());
  const DoubleTab& xp=la_zone.xp(); // Recuperation des centre de gravite des elements pour P0

  ch_T_.valeur().mettre_a_jour(temps);
  ch_T_.valeur().valeur_aux(xp, T_);

  ch_C_.valeur().mettre_a_jour(temps);
  ch_C_.valeur().valeur_aux(xp, C_);

  // mettre a jour le champ kappa (seul membrane) P0
  DoubleTab& kappa = kappa_.valeurs();
  const DoubleTab& diffu = mon_probleme().milieu().diffusivite().valeurs();				// reprende champ diffusivite
  for (int poly = 0; poly < CL_a_.valeur().nb_elem_tot(); poly++)
    {
      int elem = CL_a_.valeur()(poly);
      kappa(elem) = diffu(elem,0);
    }
  for (int poly = 0; poly < CL_c_.valeur().nb_elem_tot(); poly++)
    {
      int elem = CL_c_.valeur()(poly);
      kappa(elem) = diffu(elem,0);
    }
  for (int poly = 0; poly < MB_.valeur().nb_elem_tot(); poly++)
    {
      int elem = MB_.valeur()(poly);
      kappa(elem) = f_kappa(T_(elem), C_(elem));
    }

  // mettre a jour le champ I_i_ P0
  DoubleTab& I_i = I_i_.valeurs();
  // calcul du gradient
  DoubleTab grad(0, 1, dimension);
  equation().zone_dis().zone().creer_tableau_elements(grad);
  const DoubleTab& nu=kappa_.valeurs();
  Champ_P1NC::calcul_gradient(ch_phi_.valeur().valeurs(),grad,ref_cast(Zone_Cl_VEF,equation().zone_Cl_dis().valeur()));
  grad.echange_espace_virtuel();
  for (int elem=0; elem<grad.dimension(0); elem++)
    {
      for (int j=0; j<dimension; j++)
        {
          I_i(elem,j)=-nu(elem)*grad(elem,0,j);
        }
    }

  Cerr << "Loi_Fermeture_transport_ionique::mettre_a_jour" << finl;
}




