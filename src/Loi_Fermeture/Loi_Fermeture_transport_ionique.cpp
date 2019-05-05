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
  return is;
}

void Loi_Fermeture_transport_ionique::set_param(Param& param)
{
  param.ajouter("T_0", &T_0_, Param::REQUIRED);
  param.ajouter("CSO3", &C_SO3_, Param::REQUIRED);
  param.ajouter("por_naf", &por_naf_, Param::REQUIRED);
  param.ajouter("eps_naf", &eps_naf_, Param::REQUIRED);
  param.ajouter("tor_naf", &tor_naf_, Param::REQUIRED);
  param.ajouter("nom_ssz_CLc", &nom_ssz_CLc_, Param::REQUIRED);
  param.ajouter("nom_ssz_CLa", &nom_ssz_CLa_, Param::REQUIRED);
  param.ajouter("nom_pb_phi", &nom_pb_phi_, Param::REQUIRED);
  param.ajouter("nom_champ_phi", &nom_champ_phi_, Param::REQUIRED);
  param.ajouter("nom_pb_psi", &nom_pb_psi_, Param::REQUIRED);
  param.ajouter("nom_champ_psi", &nom_champ_psi_, Param::REQUIRED);
  param.ajouter("nom_pb_T", &nom_pb_T_, Param::REQUIRED);
  param.ajouter("nom_champ_T", &nom_champ_T_, Param::REQUIRED);
  param.ajouter("nom_pb_C_H2", &nom_pb_C_H2_, Param::REQUIRED);
  param.ajouter("nom_champ_C_H2", &nom_champ_C_H2_, Param::REQUIRED);
  param.ajouter("nom_pb_C_O2", &nom_pb_C_O2_, Param::REQUIRED);
  param.ajouter("nom_champ_C_O2", &nom_champ_C_O2_, Param::REQUIRED);
  param.ajouter("nom_pb_C_H2O", &nom_pb_C_H20_, Param::REQUIRED);
  param.ajouter("nom_champ_C_H2O", &nom_champ_C_H20_, Param::REQUIRED);
}


void Loi_Fermeture_transport_ionique::discretiser( const Discretisation_base& dis)
{
  Loi_Fermeture_base::discretiser(dis);
  ref_equation_=mon_probleme().get_equation_by_name("Conduction");
  dis.discretiser_champ("champ_elem",equation().zone_dis().valeur(),"kappa","unit", 1 ,0. , ch_kappa_);
  champs_compris_.ajoute_champ(ch_kappa_);

  dis.discretiser_champ("champ_elem",equation().zone_dis().valeur(),"I_i","unit", dimension ,0. , ch_Ii_);
  ch_Ii_ -> fixer_nature_du_champ(vectoriel);
  champs_compris_.ajoute_champ(ch_Ii_);

  dis.discretiser_champ("champ_elem",equation().zone_dis().valeur(),"Erev","unit", 1 , 0. , ch_erev_);
  champs_compris_.ajoute_champ(ch_erev_);
  dis.discretiser_champ("champ_elem",equation().zone_dis().valeur(),"eta","unit", 1 , 0. , ch_eta_);
  champs_compris_.ajoute_champ(ch_eta_);
  dis.discretiser_champ("champ_elem",equation().zone_dis().valeur(),"io","unit", 1 , 0. , ch_io_);
  champs_compris_.ajoute_champ(ch_io_);
  dis.discretiser_champ("champ_elem",equation().zone_dis().valeur(),"ir","unit", 1 , 0. , ch_ir_);
  champs_compris_.ajoute_champ(ch_ir_);
  dis.discretiser_champ("champ_elem",equation().zone_dis().valeur(),"ip","unit", 1 , 0. , ch_ip_);
  champs_compris_.ajoute_champ(ch_ip_);
  dis.discretiser_champ("champ_elem",equation().zone_dis().valeur(),"Q_reac","unit", 1 , 0. , ch_q_reac_);
  champs_compris_.ajoute_champ(ch_q_reac_);
  dis.discretiser_champ("champ_elem",equation().zone_dis().valeur(),"Q_perm","unit", 1 , 0. , ch_q_perm_);
  champs_compris_.ajoute_champ(ch_q_perm_);
}

inline void Loi_Fermeture_transport_ionique::completer()
{
  Loi_Fermeture_base::completer();
  // get the reference to the coupling fields
  Probleme_base& pb_T = ref_cast(Probleme_base,interprete().objet(nom_pb_T_));
  ch_T_ = pb_T.get_champ(nom_champ_T_);
  assert(ch_T_.valeur().que_suis_je().find("P1NC") !=-1);

  Probleme_base& pb_H2 = ref_cast(Probleme_base,interprete().objet(nom_pb_C_H2_));
  ch_C_H2_ = pb_H2.get_champ(nom_champ_C_H2_);
  assert(ch_C_H2_.valeur().que_suis_je().find("P1NC") !=-1);

  Probleme_base& pb_O2 = ref_cast(Probleme_base,interprete().objet(nom_pb_C_O2_));
  ch_C_O2_ = pb_O2.get_champ(nom_champ_C_O2_);
  assert(ch_C_O2_.valeur().que_suis_je().find("P1NC") !=-1);

  Probleme_base& pb_H2O = ref_cast(Probleme_base,interprete().objet(nom_pb_C_H20_));
  ch_C_H20_ = pb_H2O.get_champ(nom_champ_C_H20_);
  assert(ch_C_H20_.valeur().que_suis_je().find("P1NC") !=-1);

  Probleme_base& pb_psi = ref_cast(Probleme_base,interprete().objet(nom_pb_psi_));
  ch_psi_ = pb_psi.get_champ(nom_champ_psi_);
  assert(ch_psi_.valeur().que_suis_je().find("P1NC") !=-1);

  Probleme_base& pb_phi = ref_cast(Probleme_base,interprete().objet(nom_pb_phi_));
  ch_phi_ = pb_phi.get_champ(nom_champ_phi_);
  assert(ch_phi_.valeur().que_suis_je().find("P1NC") !=-1);

}

inline void Loi_Fermeture_transport_ionique::mettre_a_jour(double temps)
{
  // mettre a jour les champs couples -> interpoler vers P0
  const Zone_VF& la_zone = ref_cast(Zone_VF, equation().zone_dis().valeur());
  const DoubleTab& xp=la_zone.xp(); // Recuperation des centre de gravite des elements pour P0

  // mettre a jour 4 tableaux de valeurs du champ couple
  ch_T_.valeur().mettre_a_jour(temps);
  ch_C_H2_.valeur().mettre_a_jour(temps);
  ch_C_O2_.valeur().mettre_a_jour(temps);
  ch_C_H20_.valeur().mettre_a_jour(temps);
  ch_psi_.valeur().mettre_a_jour(temps);
  ch_phi_.valeur().mettre_a_jour(temps);

  // interpolation vers P0
  ch_T_.valeur().valeur_aux( xp, T_ );			// T
  ch_C_H2_.valeur().valeur_aux( xp, C_H2_ );		// C_H2
  ch_C_O2_.valeur().valeur_aux( xp, C_O2_ );		// C_O2
  ch_C_H20_.valeur().valeur_aux( xp, C_H2O_ );    // C_H2O
  ch_psi_.valeur().valeur_aux( xp, psi_ );		// psi
  ch_phi_.valeur().valeur_aux( xp, phi_ );		// phi

  // mettre a jour le champ kappa (seul membrane) P0
  DoubleTab& kappa = ch_kappa_.valeurs();
  int nb_elem = kappa.dimension(0);
  for (int elem = 0; elem < nb_elem; ++elem)
    {
      kappa(elem) = f_kappa(T_(elem), C_H2O_(elem));
    }

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

  // mettre a jour des champs compris
  DoubleTab& val_eta = ch_eta_.valeurs();
  DoubleTab& val_erev = ch_erev_.valeurs();
  DoubleTab& val_i0 = ch_io_.valeurs();
  DoubleTab& val_ir = ch_ir_.valeurs();
  DoubleTab& val_ip = ch_ip_.valeurs();
  DoubleTab& val_q_reac = ch_q_reac_.valeurs();
  DoubleTab& val_q_perm = ch_q_perm_.valeurs();

  for (int poly = 0; poly < CL_a_.valeur().nb_elem_tot(); poly++)
    {
      int elem = CL_a_.valeur()(poly);
      double T = T_(elem);
      double a_H2 = C_H2_(elem)/(f_Henry_H2(T)*P_ref);
      double ld = C_H2O_(elem)/C_SO3;
      if(ld <= 0.)
        {
          Cerr << "ld <= 0. -> set to a_lim" << finl;
          ld = a_lim;
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

      val_erev(elem)  = erev;
      val_eta(elem) 	= eta;
      val_i0(elem)	= i0;
      val_ir(elem) 	= ir;
      val_ip(elem)	= ip;
      val_q_reac(elem)= q_reac;
      val_q_perm(elem)= q_perm;
    }
  for (int poly = 0; poly < CL_c_.valeur().nb_elem_tot(); poly++)
    {
      int elem = CL_c_.valeur()(poly);
      double T = T_(elem);
      double a_O2 = C_O2_(elem)/(f_Henry_O2(T)*P_ref);
      double ld = C_H2O_(elem)/C_SO3;
      if(ld <= 0.)
        {
          Cerr << "ld <= 0. -> set to a_lim" << finl;
          ld = a_lim;
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
      val_erev(elem) 	= erev;
      val_eta(elem) 	= eta;
      val_i0(elem)		= i0;
      val_ir(elem) 		= ir;
      val_ip(elem)		= ip;
      val_q_reac(elem)	= q_reac;
      val_q_perm(elem)	= q_perm;

    }

  // DEBUG
  //Cerr << "champ_erev" << val_erev_ << finl;
  //Cerr << "champ_psi" << psi_ << finl;
  //Cerr << "champ_phi" << phi_ << finl;
  //Cerr << "champ_eta" << val_eta << finl;

  Cerr << "Loi_Fermeture_transport_ionique::mettre_a_jour" << finl;
}

// -dG/(n*F)-R*T*log(a_X2_lim^nu_H2*a_H^nu_H_a)/(n*F)
double Loi_Fermeture_transport_ionique::eval_erev_anode(double T, double a_H2, double a_H)
{
  double nF = 2. * 96500.;
  double RTsurnF = 8.314 * T / nF;
  double a_H2_lim = max(a_H2, a_lim);
  double RTsurnFlog = RTsurnF * log(pow(a_H2_lim, nu_H2)*pow(a_H,nu_H_a));
  //double dHox0_a = 25e3; 			// TO-DO: dHox0_a   = 25e3; dHox0_c   = 167.9e3;
  //double dSox0_a = -172; 			// TO-DO: dSox0_a   = -172; dSox0_c   = -205.6;
  double dG_0 = dHox0_a - T * dSox0_a;
  return -dG_0/nF - RTsurnFlog;
}

// -dG/(n*F)-R*T*log(a_X2_lim^nu_O2*a_H^nu_H_c*a_H2O^nu_H2O)/(n*F)
double Loi_Fermeture_transport_ionique::eval_erev_cathode(double T, double a_O2, double a_H2O, double a_H)
{
  double nF = 2. * 96500.;
  double RTsurnF = 8.314 * T / nF;
  double a_O2_lim = max(a_O2, a_lim);
  double a_H2O_lim = max(a_H2O,a_lim);
  double aij = pow(a_O2_lim, nu_O2)*pow(a_H,nu_H_c)*pow(a_H2O_lim,nu_H2O);
  double loga = log(aij);
  double RTsurnFlog = RTsurnF * loga;
  //double dHox0_c = 167.9e3; 			// TO-DO: dHox0_a   = 25e3; dHox0_c   = 167.9e3;
  //double dSox0_c = -205.6; 			// TO-DO: dSox0_a   = -172; dSox0_c   = -205.6;
  double dG_0 = dHox0_c - T * dSox0_c;
  return -dG_0/nF - RTsurnFlog;
}

double Loi_Fermeture_transport_ionique::eval_eta(double psi, double phi, double erev)
{
  return psi - phi - erev;
}

// n*F*kox0^(1-alpha)*kred0^alpha
double Loi_Fermeture_transport_ionique::eval_i0_anode(double T, double a_H2, double a_H)
{
  double k0 = kB*T/(s0*NA*h);
  double RT = R * T;
  double dGox0_a = dHox0_a - T*dSox0_a;
  double dG0_a = dH0_a - T * dS0_a;
  double dGred0 = dGox0_a + dG0_a;
  double kox0 = k0*T*exp(-dGox0_a/RT);
  double kred0 = k0*T*exp(-dGred0/RT);
  double i00 = n_a*F*pow(kox0, (1-alpha_a))*pow(kred0,alpha_a);
  double a_H2_lim = max(a_H2, a_lim);
  // i00*max(a_X2/a_X2_lim,0)*a_X2_lim^((1-alpha_a)*nu_H2)*a_H^(-alpha_a*nu_H_a)
  return i00*max(a_H2/a_H2_lim,0.)*pow(a_H2_lim, (1-alpha_a)*nu_H2)*pow(a_H, -alpha_a*nu_H_a);
}

double Loi_Fermeture_transport_ionique::eval_i0_cathode(double T, double a_O2, double a_H2O, double a_H)
{
  double k0 = kB*T/(s0*NA*h);
  double RT = R * T;
  double dGox0_c = dHox0_c - T*dSox0_c;
  double dG0_c = dH0_c - T * dS0_c;
  double dGred0 = dGox0_c + dG0_c;
  double kox0 = k0*T*exp(-dGox0_c/RT);
  double kred0 = k0*T*exp(-dGred0/RT);
  double i00 = n_c*F*pow(kox0, (1-alpha_a))*pow(kred0,alpha_a);
  double a_O2_lim = max(a_O2, a_lim);
  double a_H2O_lim = max(a_H2O,a_lim);
  // i00*max(a_X2/a_X2_lim,0)*a_X2_lim^(-alpha_c*nu_O2)*a_H^(-alpha_c*nu_H_c)*a_H2O^((1-alpha_c)*nu_H2O)
  return i00*max(a_O2/a_O2_lim,0.)*pow(a_O2_lim, (1-alpha_c)*nu_O2)*pow(a_H, -alpha_c*nu_H_c)*pow(a_H2O_lim, (1-alpha_c)*nu_H2O);
}

double Loi_Fermeture_transport_ionique::eval_ir_anode(double io, double eta, double T)
{
  double nFsurRT = n_a * F / (R * T);
  double res = 0;
  double x1 = alpha_a * nFsurRT * eta;
  if (x1>50)
    x1=50;													// VERIFIER
  res += exp(x1);
  double x2 = -(1 - alpha_a) * nFsurRT * eta;
  if (x2>50)
    x2=50;													// VERIFIER
  res -= exp(x2);
  res *= io * gamma_CL;
  return res;
}

double Loi_Fermeture_transport_ionique::eval_ir_cathode(double io, double eta, double T)
{
  double nFsurRT = n_c * F / (R * T);
  double res = 0;
  double x1 = alpha_c * nFsurRT * eta;
  if (x1>50)
    x1=50;												// VERIFIER
  res += exp(x1);
  double x2 = -(1 - alpha_c) * nFsurRT * eta;
  if (x2>50)
    x2=50;													// VERIFIER
  res -= exp(x2);
  res *= io * gamma_CL;
  return res;
}

double Loi_Fermeture_transport_ionique::eval_q_reac_anode(double psi, double phi, double ir)
{
  return (psi - phi + dH0_a / (n_a*F))*ir;
}

double Loi_Fermeture_transport_ionique::eval_q_reac_cathode(double psi, double phi, double ir)
{
  return (psi - phi + dH0_c / (n_c*F))*ir;
}

double Loi_Fermeture_transport_ionique::eval_q_perm_anode(double ip)
{
  return -(dH0_a+dH0_c)/(n_a*F)*ip;
}

double Loi_Fermeture_transport_ionique::eval_q_perm_cathode(double ip)
{
  return -(dH0_a+dH0_c)/(n_c*F)*ip;
}


