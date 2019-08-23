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
#include <Sous_Zone.h>
#include <Zone_VEF.h>
#include <Zone.h>
#include <Scatter.h>
#include <PEMFC_ToolBox.h>

Implemente_instanciable_sans_constructeur( Loi_Fermeture_transport_ionique, "Loi_Fermeture_transport_ionique", Loi_Fermeture_base ) ;

Loi_Fermeture_transport_ionique::Loi_Fermeture_transport_ionique( void )
{
  newton_max_iter_ = 1000;
  newton_threshold_ = 1e-6;
  minimal_perturbation_value_ = 1e-12;
  relative_perturbation_for_derivatives_ = 0.01 ;
}

Sortie& Loi_Fermeture_transport_ionique::printOn( Sortie& os ) const
{
  Loi_Fermeture_base::printOn( os );
  return os;
}

Entree& Loi_Fermeture_transport_ionique::readOn( Entree& is )
{
  Loi_Fermeture_base::readOn( is );

  const Domaine& dom = mon_probleme().domaine();
//  dom = equation().probleme().domaine();
  CL_a_ = dom.ss_zone(nom_ssz_CLa_);
  CL_c_ = dom.ss_zone(nom_ssz_CLc_);
  dom.creer_tableau_elements(T_);
  dom.creer_tableau_elements(C_H2_);
  dom.creer_tableau_elements(C_O2_);
  dom.creer_tableau_elements(C_H2O_);
  dom.creer_tableau_elements(psi_);
  dom.creer_tableau_elements(phi_);

  assert(nom_pb_T_ !="??" || temperature_.non_nul());
  assert(nom_pb_C_O2_ !="??" || Co_.non_nul());
  assert(nom_pb_C_H2_ !="??" || Ch_.non_nul());
  assert(nom_pb_C_H2O_ !="??" || Ce_.non_nul());

  return is;
}

void Loi_Fermeture_transport_ionique::set_param(Param& param)
{
  //param.ajouter("T_0", &T_0_, Param::OPTIONAL);
  //param.ajouter("CSO3", &C_SO3_, Param::OPTIONAL);
  param.ajouter("por_naf", &por_naf_, Param::REQUIRED);
  param.ajouter("eps_naf", &eps_naf_, Param::REQUIRED);
  param.ajouter("tor_naf", &tor_naf_, Param::REQUIRED);
  param.ajouter("nom_ssz_CLc", &nom_ssz_CLc_, Param::REQUIRED);
  param.ajouter("nom_ssz_CLa", &nom_ssz_CLa_, Param::REQUIRED);
  param.ajouter("nom_pb_phi", &nom_pb_phi_, Param::REQUIRED);
  param.ajouter("nom_champ_phi", &nom_champ_phi_, Param::REQUIRED);
  param.ajouter("nom_pb_psi", &nom_pb_psi_, Param::REQUIRED);
  param.ajouter("nom_champ_psi", &nom_champ_psi_, Param::REQUIRED);

  param.ajouter("nom_pb_T", &nom_pb_T_, Param::OPTIONAL);
  param.ajouter("nom_champ_T", &nom_champ_T_, Param::OPTIONAL);
  param.ajouter("nom_pb_C_H2", &nom_pb_C_H2_, Param::OPTIONAL);
  param.ajouter("nom_champ_C_H2", &nom_champ_C_H2_, Param::OPTIONAL);
  param.ajouter("nom_pb_C_O2", &nom_pb_C_O2_, Param::OPTIONAL);
  param.ajouter("nom_champ_C_O2", &nom_champ_C_O2_, Param::OPTIONAL);
  param.ajouter("nom_pb_C_H2O", &nom_pb_C_H2O_, Param::OPTIONAL);
  param.ajouter("nom_champ_C_H2O", &nom_champ_C_H2O_, Param::OPTIONAL);

  param.ajouter("temperature", &temperature_, Param::OPTIONAL);
  param.ajouter("Co", &Co_, Param::OPTIONAL);
  param.ajouter("Ch", &Ch_, Param::OPTIONAL);
  param.ajouter("Ce", &Ce_, Param::OPTIONAL);

  param.ajouter("i0_field_override", &i0_field_override_, Param::OPTIONAL);
  param.ajouter("Erev_field_override", &Erev_field_override_, Param::OPTIONAL);
  param.ajouter("newton_max_iter", &newton_max_iter_, Param::OPTIONAL);
  param.ajouter("newton_threshold", &newton_threshold_, Param::OPTIONAL);
  param.ajouter("minimal_perturbation_value", &minimal_perturbation_value_, Param::OPTIONAL);
  param.ajouter("relative_perturbation_for_derivatives", &relative_perturbation_for_derivatives_, Param::OPTIONAL);
}


void Loi_Fermeture_transport_ionique::discretiser( const Discretisation_base& dis)
{
  Loi_Fermeture_base::discretiser(dis);

  ref_equation_=mon_probleme().get_equation_by_name("Conduction");
  const Zone_dis_base& zd = equation().zone_dis().valeur();

  dis.discretiser_champ("champ_elem", zd, "kappa","S/m", 1 ,0. , ch_kappa_);
  champs_compris_.ajoute_champ(ch_kappa_);

  dis.discretiser_champ("champ_elem",zd,"I_i","A", dimension ,0. , ch_Ii_);
  ch_Ii_ -> fixer_nature_du_champ(vectoriel);
  champs_compris_.ajoute_champ(ch_Ii_);

  dis.discretiser_champ("champ_elem",zd,"Erev","V", 1 , 0. , ch_erev_);
  champs_compris_.ajoute_champ(ch_erev_);

  dis.discretiser_champ("champ_elem",zd,"eta","V", 1 , 0. , ch_eta_);
  champs_compris_.ajoute_champ(ch_eta_);

  dis.discretiser_champ("champ_elem",zd,"io","A/m3", 1 , 0. , ch_io_);
  champs_compris_.ajoute_champ(ch_io_);

  dis.discretiser_champ("champ_elem",zd,"ir","A/m3", 1 , 0. , ch_ir_);
  champs_compris_.ajoute_champ(ch_ir_);

  dis.discretiser_champ("champ_elem",zd,"ip","A/m3", 1 , 0. , ch_ip_);
  champs_compris_.ajoute_champ(ch_ip_);

  dis.discretiser_champ("champ_elem",zd,"Q_reac","W/m3", 1 , 0. , ch_q_reac_);
  champs_compris_.ajoute_champ(ch_q_reac_);

  dis.discretiser_champ("champ_elem",zd,"Q_perm","W/m3", 1 , 0. , ch_q_perm_);
  champs_compris_.ajoute_champ(ch_q_perm_);

  dis.discretiser_champ("champ_elem",equation().zone_dis().valeur(),"dirdphi","unit", 1 , 0. , ch_dir_dphi_);
  champs_compris_.ajoute_champ(ch_dir_dphi_);

  dis.discretiser_champ("temperature",zd,"Erev_","V", 1 , 0. , Erev_field_);
  champs_compris_.ajoute_champ(Erev_field_);

  dis.discretiser_champ("temperature",zd,"ir_","A/m3", 1 , 0. , ir_field_);
  champs_compris_.ajoute_champ(ir_field_);

  dis.discretiser_champ("temperature",zd,"i0_","A/m3", 1 , 0. , i0_field_);
  champs_compris_.ajoute_champ(i0_field_);

  dis.discretiser_champ("champ_elem",zd,"DirDcO2","unit", 1 , 0. , ch_DirDcO2_);
  champs_compris_.ajoute_champ(ch_DirDcO2_);

  dis.discretiser_champ("champ_elem",zd,"DirDcH2O","unit", 1 , 0. , ch_DirDcH2O_);
  champs_compris_.ajoute_champ(ch_DirDcH2O_);

  dis.discretiser_champ("champ_elem",zd,"DirDcH2","unit", 1 , 0. , ch_DirDcH2_);
  champs_compris_.ajoute_champ(ch_DirDcH2_);

  dis.discretiser_champ("champ_elem",zd,"DQreacDT","unit", 1 , 0. , ch_DQreacDT_);
  champs_compris_.ajoute_champ(ch_DQreacDT_);
}


// If the pbName is not "none", sets the refField reference to the requested field
void Loi_Fermeture_transport_ionique::fetch_field(const Nom& pbName, const Nom& fieldName, REF(Champ_base) &refField)
{

  if (Motcle(pbName) != "none")
    //if (pbName != "??")
    {
      Probleme_base& pb = ref_cast(Probleme_base,interprete().objet(pbName));
      refField = pb.get_champ(fieldName);
      assert(refField.valeur().que_suis_je().find("P1NC") !=-1);
    }
}

inline void Loi_Fermeture_transport_ionique::completer()
{
  Loi_Fermeture_base::completer();

  // get the reference to the coupling fields
  fetch_field(nom_pb_T_, nom_champ_T_, ch_T_);
  fetch_field(nom_pb_C_H2_, nom_champ_C_H2_, ch_C_H2_);
  fetch_field(nom_pb_C_O2_, nom_champ_C_O2_, ch_C_O2_);
  fetch_field(nom_pb_C_H2O_, nom_champ_C_H2O_, ch_C_H2O_);
  fetch_field(nom_pb_psi_, nom_champ_psi_, ch_psi_);
  fetch_field(nom_pb_phi_, nom_champ_phi_, ch_phi_);

  //  Initialize the "common_faces list"
  Probleme_base& pb_phi = ref_cast(Probleme_base,interprete().objet(nom_pb_phi_));
  Probleme_base& pb_psi = ref_cast(Probleme_base,interprete().objet(nom_pb_psi_));
  Probleme_base& pb_T = ref_cast(Probleme_base,interprete().objet(nom_pb_T_));
  const Zone_VEF& zone_psi = ref_cast(Zone_VEF, pb_psi.domaine_dis().zone_dis(0).valeur());
  const Zone_VEF& zone_phi = ref_cast(Zone_VEF, pb_phi.domaine_dis().zone_dis(0).valeur());
  const Zone_VEF& zone_T = ref_cast(Zone_VEF, pb_T.domaine_dis().zone_dis(0).valeur());

  IntTab commonFacesPsiPhi;
  PEMFC_ToolBox::find_common_faces_indices(zone_psi, zone_phi, commonFacesPsiPhi);
  // for each face i2 in zone_phi, index1in2[i2] = index of same face in zone_T:
  ArrOfInt index1in2;
  PEMFC_ToolBox::find_faces_indices1in2(zone_T, zone_phi, index1in2);
  const int ncommon = commonFacesPsiPhi.dimension(0);
  face_index_psi_.resize_array(ncommon);
  face_index_phi_.resize_array(ncommon);
  face_index_T_.resize_array(ncommon);
  face_is_anode_.resize_array(ncommon);
  ArrOfInt tag_anode(zone_phi.nb_faces_tot());
  {
    // tag faces that belong to the anode zone:
    const Sous_Zone& sszone = zone_phi.zone().domaine().ss_zone(nom_ssz_CLa_);
    int nelem = sszone.nb_elem_tot();
    const IntTab& elemfaces = zone_phi.elem_faces();
    int nfaces_elem = elemfaces.dimension(1); // n faces per element
    for (int i = 0; i < nelem; i++)
      {
        for (int j = 0; j < nfaces_elem; j++)
          {
            int face = elemfaces(sszone[i],j);
            tag_anode[face] = 1;
          }
      }
  }
  for (int i = 0; i < ncommon; i++)
    {
      int iphi = commonFacesPsiPhi(i,1);
      face_index_psi_[i] = commonFacesPsiPhi(i,0);
      face_index_phi_[i] = iphi;
      face_index_T_[i] = index1in2[iphi];
      if (face_index_T_[i]<0)
        {
          double x = zone_phi.xv(iphi,0);
          double y = zone_phi.xv(iphi,1);
          double z = (dimension==3)?zone_phi.xv(iphi,2):0.;
          Cerr << "Error in Loi_Fermeture_transport_ionique::completer():\n"
               << "face " << iphi << "at xyz=(" << x << " " << y << " " << z << ")  not found in temperature domain"<<endl;
        }
      face_is_anode_[i] = tag_anode[iphi];
    }
}

inline void Loi_Fermeture_transport_ionique::mettre_a_jour(double temps)
{
  // mettre a jour les champs couples -> interpoler vers P0
  const Zone_VF& la_zone = ref_cast(Zone_VF, equation().zone_dis().valeur());
  const DoubleTab& xp=la_zone.xp(); // Recuperation des centre de gravite des elements pour P0

  // mettre a jour 4 tableaux de valeurs du champ couple
  /*
  ch_T_.valeur().mettre_a_jour(temps);
  ch_C_H2_.valeur().mettre_a_jour(temps);
  ch_C_O2_.valeur().mettre_a_jour(temps);
  ch_C_H2O_.valeur().mettre_a_jour(temps);
  ch_psi_.valeur().mettre_a_jour(temps);
  ch_phi_.valeur().mettre_a_jour(temps);
  */
  // interpolation vers P0
  if (ch_T_.non_nul())
    {
      ch_T_.valeur().valeur_aux( xp, T_ );
      T_.echange_espace_virtuel();
    }
  if(temperature_.non_nul())
    {
      temperature_.mettre_a_jour(temps);	// // update field (in case it is time dependent)
      temperature_.valeur().valeur_aux( xp, T_ );
      T_.echange_espace_virtuel();
    }
  if (ch_C_H2_.non_nul())
    {
      ch_C_H2_.valeur().valeur_aux( xp, C_H2_ );
      C_H2_.echange_espace_virtuel();
    }
  if(Ch_.non_nul())
    {
      Ch_.mettre_a_jour(temps);
      Ch_.valeur().valeur_aux( xp, C_H2_ );
      C_H2_.echange_espace_virtuel();
    }
  if (ch_C_O2_.non_nul())
    {
      ch_C_O2_.valeur().valeur_aux( xp, C_O2_ );
      C_O2_.echange_espace_virtuel();
    }
  if (Co_.non_nul())
    {
      Co_.mettre_a_jour(temps);
      Co_.valeur().valeur_aux( xp, C_O2_ );
      C_O2_.echange_espace_virtuel();
    }
  if (ch_C_H2O_.non_nul())
    {
      ch_C_H2O_.valeur().valeur_aux( xp, C_H2O_ );
      C_H2O_.echange_espace_virtuel();
    }
  if (Ce_.non_nul())
    {
      Ce_.mettre_a_jour(temps);
      Ce_.valeur().valeur_aux( xp, C_H2O_ );
      C_H2O_.echange_espace_virtuel();
    }

  ch_psi_.valeur().valeur_aux( xp, psi_ );
  psi_.echange_espace_virtuel();
  ch_phi_.valeur().valeur_aux( xp, phi_ );
  phi_.echange_espace_virtuel();

  // Overriding the value of Erev if requested in datafile (for testing purpose)
  if (Erev_field_override_.non_nul())
    {
      // update field (in case it is time dependent)
      Erev_field_.valeur().mettre_a_jour(temps);
      // evaluate field at the position of unknowns:
      Erev_field_.valeur().affecter(Erev_field_override_);
    }
  if (i0_field_override_.non_nul())
    {
      // update field (in case it is time dependent)
      i0_field_.valeur().mettre_a_jour(temps);
      // evaluate field at the position of unknowns:
      i0_field_.valeur().affecter(i0_field_override_);
    }

  // mettre a jour le champ kappa (seul membrane) P0
  DoubleTab& kappa = ch_kappa_.valeurs();
  DoubleTab& por = por_naf_.valeurs();
  DoubleTab& eps = eps_naf_.valeurs();
  DoubleTab& tor = tor_naf_.valeurs();

  int nb_elem = kappa.dimension(0);
  for (int elem = 0; elem < nb_elem; ++elem)
    {
      kappa(elem) = f_kappa(T_(elem), C_H2O_(elem), por(elem, 0), eps(elem,0), tor(elem,0));
    }
  ch_kappa_.mettre_a_jour(temps);

  // mettre a jour le champ I_i_ P0
  DoubleTab& I_i = ch_Ii_.valeurs();
  // calcul du gradient
  DoubleTab grad(0, 1, dimension);
  equation().zone_dis().zone().creer_tableau_elements(grad);
  const DoubleTab& nu=ch_kappa_.valeurs();
  Champ_P1NC::calcul_gradient(ch_phi_.valeur().valeurs(),grad,ref_cast(Zone_Cl_VEF,equation().zone_Cl_dis().valeur()));
  grad.echange_espace_virtuel();
  for (int elem=0; elem<grad.dimension(0); elem++)
    {
      for (int j=0; j<dimension; j++)
        {
          I_i(elem,j)=-nu(elem)*grad(elem,0,j);
        }
    }
  ch_Ii_.mettre_a_jour(temps);

  // mettre a jour des champs compris
  DoubleTab& val_eta = ch_eta_.valeurs();
  DoubleTab& val_erev = ch_erev_.valeurs();
  DoubleTab& val_i0 = ch_io_.valeurs();
  DoubleTab& val_ir = ch_ir_.valeurs();
  DoubleTab& val_ip = ch_ip_.valeurs();
  DoubleTab& val_q_reac = ch_q_reac_.valeurs();
  DoubleTab& val_q_perm = ch_q_perm_.valeurs();
  DoubleTab& val_dir_dphi = ch_dir_dphi_.valeurs();

  DoubleTab& val_DirDcH2O = ch_DirDcH2O_.valeurs();
  DoubleTab& val_DirDcH2 = ch_DirDcH2_.valeurs();
  DoubleTab& val_DirDcO2 = ch_DirDcO2_.valeurs();

  DoubleTab& val_DQreacDT = ch_DQreacDT_.valeurs();

  for (int poly = 0; poly < CL_a_.valeur().nb_elem_tot(); poly++)
    {
      int elem = CL_a_.valeur()(poly);
      double T = T_(elem);
      double a_H2 = C_H2_(elem)/(f_Henry_H2(T)*P_ref);
      double ld = C_H2O_(elem)/C_SO3;
      if(ld <= 0.)
        {
          //Cerr << "ld <= 0. -> set to a_lim" << finl;
          ld = a_lim;
          Cerr << "Concentration H2O negative " << finl;
          abort();
        }
      double a_H = f_lambda(1.) / (ld);							// ATTENTION divise by zero
      double phi = phi_(elem);
      double psi = psi_(elem);

      double erev = eval_erev_anode(T,a_H2, a_H);
      double eta  = eval_eta(psi, phi , erev);
      double i0   = eval_i0_anode(T, a_H2, a_H);
      double ir	  = eval_ir_anode(i0, eta, T);
      double ip   = 0;											// TO-DO: need to evaluating
      double q_reac = eval_q_reac_anode(psi, phi, ir);
      double q_perm = eval_q_perm_anode(ip);

      double dirdphi = eval_dirdphi_anode(i0, eta, T);

      // derivatives evaluation
      double DirDcH2 = 0.0 ;
      eval_derivatives_anode_H2(elem, T, C_H2_(elem), a_H, DirDcH2 );

      // derivee dQreacdT
      double derevdT = eval_derevdT_anode(T,a_H2, a_H);
      double dirdT = eval_dirdT_anode(i0,eta,T, derevdT);
      double dQdT = eval_q_reac_anode(psi, phi, dirdT);

      val_erev(elem)  = erev;
      val_eta(elem) 	= eta;
      val_i0(elem)	= i0;
      val_ir(elem) 	= ir;
      val_ip(elem)	= ip;
      val_q_reac(elem)= q_reac;
      val_q_perm(elem)= q_perm;

      val_dir_dphi(elem) = dirdphi;

      val_DirDcH2 ( elem ) = DirDcH2;

      val_DQreacDT(elem) = dQdT;
    }
  for (int poly = 0; poly < CL_c_.valeur().nb_elem_tot(); poly++)
    {
      int elem = CL_c_.valeur()(poly);
      double T = T_(elem);
      double a_O2 = C_O2_(elem)/(f_Henry_O2(T)*P_ref);
      double ld = C_H2O_(elem)/C_SO3;
      if(ld <= 0.)
        {
          //Cerr << "ld <= 0. -> set to a_lim" << finl;
          ld = a_lim;
          Cerr << "Concentration H2O negative " << finl;
          abort();
        }
      double a_H2O = f_lambda_inv(ld);				// a_H20 = f_lambda_inv(ld) -> need testing -> dvq: tested OK
      double a_H = f_lambda(1.) / (ld);										// ATTENTION divise by zero
      double phi = phi_(elem);
      double psi = psi_(elem);

      double erev = eval_erev_cathode(T, a_O2 ,a_H2O, a_H);
      double eta  = eval_eta(psi, phi, erev);
      double i0   = eval_i0_cathode(T, a_O2, a_H2O, a_H);
      double ir	  = eval_ir_cathode(i0, eta, T);
      double ip   = 0;														// TO-DO: need to evaluating
      double q_reac = eval_q_reac_cathode(psi, phi, ir);
      double q_perm = eval_q_perm_cathode(ip);

      double dirdphi = eval_dirdphi_cathode(i0,eta,T);

      double DirDcO2 = 0.0 ;
      eval_derivatives_cathode_O2(elem, T, C_O2_(elem), a_H2O, a_H, DirDcO2 );
      double DirDcH2O = 0.0;
      eval_derivatives_cathode_H2O(elem, T, C_O2_(elem), C_H2O_(elem), a_O2, a_H, DirDcH2O );

      // derivee dQreacdT
      double derevdT = eval_derevdT_cathode(T, a_O2 ,a_H2O, a_H);
      double dirdT = eval_dirdT_cathode(i0,eta,T,derevdT);
      double dQdT = eval_q_reac_cathode(psi, phi, dirdT);

      val_erev(elem) 	= erev;
      val_eta(elem) 	= eta;
      val_i0(elem)	= i0;
      val_ir(elem) 	= ir;
      val_ip(elem)	= ip;
      val_q_reac(elem)	= q_reac;
      val_q_perm(elem)	= q_perm;

      val_dir_dphi(elem) = dirdphi;
      val_DirDcO2( elem ) = DirDcO2;
      val_DirDcH2O( elem ) = DirDcH2O;

      val_DQreacDT(elem) = dQdT;
    }
  ch_eta_.mettre_a_jour(temps);
  ch_erev_.mettre_a_jour(temps);
  ch_io_.mettre_a_jour(temps);
  ch_ir_.mettre_a_jour(temps);
  ch_ip_.mettre_a_jour(temps);
  ch_q_reac_.mettre_a_jour(temps);
  ch_q_perm_.mettre_a_jour(temps);
  ch_dir_dphi_.mettre_a_jour(temps);
  ch_DirDcH2O_.mettre_a_jour(temps);
  ch_DirDcH2_.mettre_a_jour(temps);
  ch_DirDcO2_.mettre_a_jour(temps);
  ch_DQreacDT_.mettre_a_jour(temps);


  Cerr << "Loi_Fermeture_transport_ionique::mettre_a_jour" << finl;
  Cerr << "champ E_rev [Volt] min max " << mp_min_vect(val_erev) << " " << mp_max_vect(val_erev) << finl;
  Cerr << "champ eta [Volt] min max " << mp_min_vect(val_eta) << " " << mp_max_vect(val_eta) << finl;
  Cerr << "champ i0 [A/m2] min max " << mp_min_vect(val_i0) << " " << mp_max_vect(val_i0) << finl;
  Cerr << "champ ir [A/m3] min max " << mp_min_vect(val_ir) << " " << mp_max_vect(val_ir) << finl;
  Cerr << "champ ip [A/m3] min max " << mp_min_vect(val_ip) << " " << mp_max_vect(val_ip) << finl;
  Cerr << "champ q_reac [W/m3] min max " << mp_min_vect(val_q_reac) << " " << mp_max_vect(val_q_reac) << finl;
  Cerr << "champ dqreacdT min max " << mp_min_vect(val_DQreacDT) << " " << mp_max_vect(val_DQreacDT) << finl;
  Cerr << "champ q_perm [W/m3] min max " << mp_min_vect(val_q_perm) << " " << mp_max_vect(val_q_perm) << finl;
}

double Loi_Fermeture_transport_ionique::compute_perturbation( const double& c , const Nom& infoString ) const
{
  double perturbation = dabs( relative_perturbation_for_derivatives_ * c );
  return perturbation;
}

void Loi_Fermeture_transport_ionique::eval_derivatives_anode_H2(const int& elem, const double& T, const double& C_H2, const double& a_H, double& Dir ) const
{

  const double perturbation = compute_perturbation( C_H2 , "ir_in_anode_for_H2" );

  DoubleTab c( 2 );
  c[ 0 ] = C_H2 + perturbation ;
  c[ 1 ] = C_H2 - perturbation ;

  DoubleTab evals( 2 ) ;

  for( int i=0; i<2; i++ )
    {
      const double a_H2 = c[ i ]/( f_Henry_H2( T ) * P_ref );

      const double phi = phi_(elem);
      const double psi = psi_(elem);

      const double erev = eval_erev_anode(T,a_H2, a_H);
      const double eta  = eval_eta(psi, phi , erev);
      const double i0  = eval_i0_anode(T, a_H2, a_H);
      evals[ i ] = eval_ir_anode(i0, eta, T);
    }

  Dir = ( evals[ 0 ] - evals[ 1 ] ) / ( 2 * perturbation );
}

void Loi_Fermeture_transport_ionique::eval_derivatives_cathode_O2(const int& elem, const double& T, const double& C_O2, const double& a_H2O, const double& a_H, double& Dir ) const
{

  const double perturbation = compute_perturbation( C_O2 , "ir_in_cathode_for_O2");

  DoubleTab c( 2 );
  c[ 0 ] = C_O2 + perturbation ;
  c[ 1 ] = C_O2 - perturbation ;

  DoubleTab evals( 2 );

  for( int i=0; i<2; i++ )
    {
      const double a_O2 = c[ i ]/(f_Henry_O2(T)*P_ref);

      const double phi = phi_(elem);
      const double psi = psi_(elem);

      const double erev = eval_erev_cathode(T, a_O2 ,a_H2O, a_H);
      const double eta  = eval_eta(psi, phi, erev);
      const double i0   = eval_i0_cathode(T, a_O2, a_H2O, a_H);
      evals[ i ] = eval_ir_cathode(i0, eta, T);
    }

  Dir = ( evals[ 0 ] - evals[ 1 ] ) / ( 2 * perturbation );


}

void Loi_Fermeture_transport_ionique::eval_derivatives_cathode_H2O(const int& elem, const double& T, const double& C_O2, const double& C_H2O, const double& a_O2, const double& a_H, double& Dir ) const
{

  const double perturbation = compute_perturbation( C_H2O ,"ir_in_cathode_for_H2O"  );

  DoubleTab c( 2 );
  c[ 0 ] = C_H2O + perturbation ;
  c[ 1 ] = C_H2O - perturbation ;

  DoubleTab evals( 2 );

  for( int i=0; i<2; i++ )
    {

      double ld = c[ i ]/C_SO3;
      if(ld <= 0.)
        {
          ld = a_lim;
          Cerr << "Concentration H2O negative " << finl;
          abort();
        }
      const double a_H2O = f_lambda_inv( ld );

      const double phi = phi_(elem);
      const double psi = psi_(elem);

      const double erev = eval_erev_cathode(T, a_O2 ,a_H2O, a_H);
      const double eta  = eval_eta(psi, phi, erev);
      const double i0   = eval_i0_cathode(T, a_O2, a_H2O, a_H);
      evals[ i ] = eval_ir_cathode(i0, eta, T);

    }

  Dir = ( evals[ 0 ] - evals[ 1 ] ) / ( 2 * perturbation );

}

// -dG/(n*F)-R*T*log(a_X2_lim^nu_H2*a_H^nu_H_a)/(n*F)
double Loi_Fermeture_transport_ionique::eval_erev_anode(double T, double a_H2, double a_H) const
{
  //double nF = 2. * 96500.;
  double RTsurnF = R * T / (n_a*F);
  //double a_H2_lim = max(a_H2, a_lim);
  //double RTsurnFlog = RTsurnF * log(pow(a_H2_lim, nu_H2)*pow(a_H,nu_H_a));
  double RTsurnFlog = RTsurnF * log(pow(a_H2, nu_H2)*pow(a_H,nu_H_a));
  //double dHox0_a = 25e3; 			// TO-DO: dHox0_a   = 25e3; dHox0_c   = 167.9e3;
  //double dSox0_a = -172; 			// TO-DO: dSox0_a   = -172; dSox0_c   = -205.6;
  double dG_0 = dH0_a - T * dS0_a;
  return -dG_0/(n_a*F) - RTsurnFlog;
}

double Loi_Fermeture_transport_ionique::eval_derevdT_anode(double T, double a_H2, double a_H) const
{
  double RsurnF = R / (n_a*F);
  //double a_H2_lim = max(a_H2, a_lim);
  //double RTsurnFlog = RsurnF * log(pow(a_H2_lim, nu_H2)*pow(a_H,nu_H_a));
  double RTsurnFlog = RsurnF * log(pow(a_H2, nu_H2)*pow(a_H,nu_H_a));
  return -RTsurnFlog;
}

double Loi_Fermeture_transport_ionique::eval_derevdT_cathode(double T, double a_O2, double a_H2O, double a_H) const
{
  double RsurnF = R / (n_c * F);
  //double a_O2_lim = max(a_O2, a_lim);
  //double a_H2O_lim = max(a_H2O,a_lim);
  //double aij = pow(a_O2_lim, nu_O2)*pow(a_H,nu_H_c)*pow(a_H2O_lim,nu_H2O);
  double aij = pow(a_O2, nu_O2)*pow(a_H,nu_H_c)*pow(a_H2O,nu_H2O);
  double loga = log(aij);
  double RTsurnFlog = RsurnF * loga;
  return -RTsurnFlog;
}

// -dG/(n*F)-R*T*log(a_X2_lim^nu_O2*a_H^nu_H_c*a_H2O^nu_H2O)/(n*F)
double Loi_Fermeture_transport_ionique::eval_erev_cathode(double T, double a_O2, double a_H2O, double a_H) const
{
  //double nF = 2. * 96500.;
  double RTsurnF = R * T / (n_c * F);
  //double a_O2_lim = max(a_O2, a_lim);
  //double a_H2O_lim = max(a_H2O,a_lim);
  //double aij = pow(a_O2_lim, nu_O2)*pow(a_H,nu_H_c)*pow(a_H2O_lim,nu_H2O);
  double aij = pow(a_O2, nu_O2)*pow(a_H,nu_H_c)*pow(a_H2O,nu_H2O);
  double loga = log(aij);
  double RTsurnFlog = RTsurnF * loga;
  //double dHox0_c = 167.9e3; 			// TO-DO: dHox0_a   = 25e3; dHox0_c   = 167.9e3;
  //double dSox0_c = -205.6; 			// TO-DO: dSox0_a   = -172; dSox0_c   = -205.6;
  double dG_0 = dH0_c - T * dS0_c;
  return -dG_0/(n_c*F) - RTsurnFlog;
}

double Loi_Fermeture_transport_ionique::eval_eta(double psi, double phi, double erev) const
{
  return psi - phi - erev;
}

// n*F*kox0^(1-alpha)*kred0^alpha
double Loi_Fermeture_transport_ionique::eval_i0_anode(double T, double a_H2, double a_H) const
{
  double k0 = kB/(s0*NA*h);
  double RT = R * T;
  double dGox0_a = dHox0_a - T*dSox0_a;
  double dG0_a = dH0_a - T * dS0_a;
  double dGred0_a = dGox0_a + dG0_a;
  double kox0_a = k0*T*exp(-dGox0_a/RT);
  double kred0_a = k0*T*exp(-dGred0_a/RT);
  double i00 = n_a*F*pow(kox0_a, (1-alpha_a))*pow(kred0_a,alpha_a);
  //double a_H2_lim = max(a_H2, a_lim);
  // i00*max(a_X2/a_X2_lim,0)*a_X2_lim^((1-alpha_a)*nu_H2)*a_H^(-alpha_a*nu_H_a)
  //return i00*max(a_H2/a_H2_lim,0.)*pow(a_H2_lim, (1-alpha_a)*nu_H2)*pow(a_H, -alpha_a*nu_H_a);
  return i00*pow(a_H2, (1-alpha_a)*nu_H2)*pow(a_H, -alpha_a*nu_H_a);
}

double Loi_Fermeture_transport_ionique::eval_i0_cathode(double T, double a_O2, double a_H2O, double a_H) const
{
  double k0 = kB/(s0*NA*h);
  double RT = R * T;
  double dGox0_c = dHox0_c - T*dSox0_c;
  double dG0_c = dH0_c - T * dS0_c;
  double dGred0_c = dGox0_c + dG0_c;
  double kox0_c = k0*T*exp(-dGox0_c/RT);
  double kred0_c = k0*T*exp(-dGred0_c/RT);
  double i00 = n_c*F*pow(kox0_c, (1-alpha_c))*pow(kred0_c,alpha_c);
  //double a_O2_lim = max(a_O2, a_lim);
  //double a_H2O_lim = max(a_H2O,a_lim);
  // i00*max(a_X2/a_X2_lim,0)*a_X2_lim^(-alpha_c*nu_O2)*a_H^(-alpha_c*nu_H_c)*a_H2O^((1-alpha_c)*nu_H2O)
  //return i00*max(a_O2/a_O2_lim,0.)*pow(a_O2_lim, (1-alpha_c)*nu_O2)*pow(a_H, -alpha_c*nu_H_c)*pow(a_H2O_lim, (1-alpha_c)*nu_H2O);
  return i00*pow(a_O2, (1-alpha_c)*nu_O2)*pow(a_H, -alpha_c*nu_H_c)*pow(a_H2O, (1-alpha_c)*nu_H2O);
}

double Loi_Fermeture_transport_ionique::eval_ir_anode(double io, double eta, double T) const
{
  double nFsurRT = n_a * F / (R * T);
  double res = 0;
  double x1 = alpha_a * nFsurRT * eta;
  //if (x1>50)
  //x1=50;													// VERIFIER
  res += exp(x1);
  double x2 = -(1 - alpha_a) * nFsurRT * eta;
  //if (x2>50)
  //x2=50;													// VERIFIER
  res -= exp(x2);
  res *= io * gamma_CL_a;
  return res;
}

double Loi_Fermeture_transport_ionique::eval_ir_cathode(double io, double eta, double T) const
{
  double nFsurRT = n_c * F / (R * T);
  double res = 0;
  double x1 = alpha_c * nFsurRT * eta;
  //if (x1>50)
  //x1=50;												// VERIFIER
  res += exp(x1);
  double x2 = -(1 - alpha_c) * nFsurRT * eta;
  //if (x2>50)
  //x2=50;													// VERIFIER
  res -= exp(x2);
  res *= io * gamma_CL_c;
  return res;
}

double Loi_Fermeture_transport_ionique::eval_q_reac_anode(double psi, double phi, double ir) const
{
  return (psi - phi + dH0_a / (n_a*F))*ir;
}

double Loi_Fermeture_transport_ionique::eval_q_reac_cathode(double psi, double phi, double ir) const
{
  return (psi - phi + dH0_c / (n_c*F))*ir;
}

double Loi_Fermeture_transport_ionique::eval_q_perm_anode(double ip) const
{
  return -(dH0_a+dH0_c)/(n_a*F)*ip;
}

double Loi_Fermeture_transport_ionique::eval_q_perm_cathode(double ip) const
{
  return -(dH0_a+dH0_c)/(n_c*F)*ip;
}

// dirdphi = -dirdpsi
double Loi_Fermeture_transport_ionique::eval_dirdphi_anode(double io, double eta, double T) const
{
  return -eval_dirdpsi_anode(io, eta, T);
}

// i0.gamma[alpha.nF/RT.exp(alpha.nF/RT.eta)+(1-alpha)nF/RT.exp(-(1-alpha)nF/RT.eta)]
double Loi_Fermeture_transport_ionique::eval_dirdpsi_anode(double io, double eta, double T) const
{
  double nFsurRT = n_a * F / (R * T);
  double res = 0;
  double x1 = alpha_a * nFsurRT * eta;
  //if (x1>50)
  //x1=50;													// VERIFIER
  res += alpha_a * nFsurRT * exp(x1);
  double x2 = -(1 - alpha_a) * nFsurRT * eta;
  //if (x2>50)
  //x2=50;													// VERIFIER
  res += (1-alpha_a)*nFsurRT * exp(x2);
  res *= io * gamma_CL_a;
  return res;
}

// dirdphi = -dirdpsi
double Loi_Fermeture_transport_ionique::eval_dirdphi_cathode(double io, double eta, double T) const
{
  return - eval_dirdpsi_cathode(io, eta, T);
}

// i0.gamma[alpha.nF/RT.exp(alpha.nF/RT.eta)+(1-alpha)nF/RT.exp(-(1-alpha)nF/RT.eta)]
double Loi_Fermeture_transport_ionique::eval_dirdpsi_cathode(double io, double eta, double T) const
{
  double nFsurRT = n_c * F / (R * T);
  double res = 0;
  double x1 = alpha_c * nFsurRT * eta;
  //if (x1>50)
  //x1=50;												// VERIFIER
  res += alpha_c * nFsurRT * exp(x1);
  double x2 = -(1 - alpha_c) * nFsurRT * eta;
  //if (x2>50)
  //x2=50;													// VERIFIER
  res += (1-alpha_c) * nFsurRT * exp(x2);
  res *= io * gamma_CL_c;
  return res;
}
// Function called by the Newton solver for the coupled equations phi-psi:
//  assumes that E_rev is uptodate and that the ir exchange term should be
//  recomputed for the given value of phi and psi. We compute also the derivative
//  to build the jacobian matrix for the coupled problem...
// Inputs: psi=values of psi on the faces of the psi domain
//  phi: same for phi, on the phi domain
// Outputs:
//  ir: exchange current*control volume on common faces to phi and psi.
//      array with same indexing as face_index_xxx_ array (indices of common faces)
//
//  DirDpsi: same structure, filled with the (d(ir)/d(psi))*control volume
void Loi_Fermeture_transport_ionique::compute_ir_DirDpsi(const DoubleTab& psi_values, const DoubleTab& phi_values, DoubleTab& ir, DoubleTab& DirDpsi)
{
  int i;
  Probleme_base& pb_T = ref_cast(Probleme_base,interprete().objet(nom_pb_T_));
  const DoubleTab& T_at_faces = pb_T.get_champ(nom_champ_T_).valeurs();
  const int nfaces = face_index_psi_.size_array();
  ir.resize(nfaces);
  DirDpsi.resize(nfaces);
  const DoubleTab& Erev_field_values = Erev_field_.valeurs();
  const DoubleTab& i0_field_values = i0_field_.valeurs();
  DoubleTab& ir_field_values = ir_field_.valeurs();  // output
  for (i = 0; i < nfaces; i++)
    {
      const int indexPhi = face_index_phi_[i];
      const int indexPsi = face_index_psi_[i];
      double T = T_at_faces[face_index_T_[i]];
      double psi = psi_values[indexPsi];
      double phi = phi_values[indexPhi];
      // Erev_field is discretized on the domain of the phi problem
      double Erev = Erev_field_values[indexPhi];
      double i0 = i0_field_values[indexPhi];
      double n, alpha, gamma_CL;
      if (face_is_anode_[i])
        {
          n = n_a;
          alpha = alpha_a;
          gamma_CL = gamma_CL_a;
        }
      else
        {
          n = n_c;
          alpha = alpha_c;
          gamma_CL = gamma_CL_c;
        }
      double eta = psi - phi - Erev; //Erev_common_faces_[i];
      double nFsurRT = n * F / (R * T);
      double x1 = alpha * nFsurRT;
      double ex1 = exp(x1*eta);
      double x2 = -(1 - alpha) * nFsurRT;
      double ex2 = exp(x2*eta);
      double ir_value = (ex1 - ex2) * i0 * gamma_CL;
      ir[i] = ir_value;
      ir_field_values[indexPhi] = ir[i];
      DirDpsi[i] = (x1 * ex1 - x2 * ex2) * i0 * gamma_CL;
    }
}

// assume que dirdT = i0.gamma.[alpha.nF/R.exp(alphanF/RT.eta)+(1-alpha)nF/R.exp(-(1-alpha)nF/RT.eta)] * d(eta/T)/dT
// d(eta/T)/dT = -eta/T^2-dErev/dT/T
double Loi_Fermeture_transport_ionique::eval_dirdT_anode(const double& io, const double& eta, const double& T, const double& dErevdT) const
{
  double nFsurR = n_a * F / R ;
  double nFsurRT = nFsurR/T;
  double res = 0;
  double x1 = alpha_a * nFsurRT * eta;
  //if (x1>50)
  //x1=50;													// VERIFIER
  res += alpha_a * nFsurR * exp(x1);
  double x2 = -(1 - alpha_a) * nFsurRT * eta;
  //if (x2>50)
  //x2=50;													// VERIFIER
  res += (1-alpha_a)*nFsurR * exp(x2);
  res *= io * gamma_CL_c;
  double dhdT = - eta/(T*T) - dErevdT/T;	// fonction h(T) = eta(T)/T
  res *= dhdT;
  return res;
}

double Loi_Fermeture_transport_ionique::eval_dirdT_cathode(const double& io,const double& eta,const double& T, const double& dErevdT) const
{
  double nFsurR = n_c * F / R ;
  double nFsurRT = nFsurR/T;
  double res = 0;
  double x1 = alpha_c * nFsurRT * eta;
  //if (x1>50)
  // x1=50;													// VERIFIER
  res += alpha_c * nFsurR * exp(x1);
  double x2 = -(1 - alpha_c) * nFsurRT * eta;
  //if (x2>50)
  //x2=50;													// VERIFIER
  res += (1-alpha_c)*nFsurR * exp(x2);
  res *= io * gamma_CL_c;
  double dhdT = - eta/(T*T)- dErevdT/T;	// fonction h(T) = eta(T)/T
  res *= dhdT;
  return res;
}
