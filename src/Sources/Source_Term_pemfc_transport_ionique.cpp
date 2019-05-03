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
// File      : Source_Term_pemfc_transport_ionique.cpp
// Directory : $PEMFC_ROOT/src/Sources
//
/////////////////////////////////////////////////////////////////////////////

#include <Source_Term_pemfc_transport_ionique.h>
#include <Zone_VEF.h>
#include <Zone_Cl_VEF.h>
#include <Probleme_base.h>
#include <Equation_base.h>
#include <Interprete.h>
#include <Domaine.h>

Implemente_instanciable( Source_Term_pemfc_transport_ionique, "Source_Term_pemfc_transport_ionique_VEF_P1NC", Source_base ) ;

Sortie& Source_Term_pemfc_transport_ionique::printOn( Sortie& os ) const
{
  Source_base::printOn( os );
  return os;
}

Entree& Source_Term_pemfc_transport_ionique::readOn( Entree& is )
{
  Source_base::readOn( is );
  Cerr << " Source_Term_pemfc_transport_ionique::readOn " << finl  ;
  Param param(que_suis_je());
  set_param(param);
  param.lire_avec_accolades(is);

  dom_ = equation().probleme().domaine();
  CL_a_ = dom_.valeur().ss_zone(nom_ssz_CLa_);
  CL_c_ = dom_.valeur().ss_zone(nom_ssz_CLc_);
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

  assert((nom_pb_psi_ != "??") || (nom_pb_phi_ != "??"));

  if(nom_pb_psi_ != "??" && (nom_pb_phi_ == "??"))
    {
      // transport ionique

      Probleme_base& pb_psi = ref_cast(Probleme_base,interprete().objet(nom_pb_psi_));
      ch_psi_ = pb_psi.get_champ(nom_champ_psi_);
      assert(ch_psi_.valeur().que_suis_je().find("P1NC") !=-1);
      ch_phi_ = equation().inconnue();
      assert(ch_phi_.valeur().que_suis_je().find("P1NC") !=-1);
      phi_.ref(ch_phi_.valeur().valeurs());
    }
  else if(nom_pb_phi_ != "??" && (nom_pb_psi_ == "??"))
    {
      // transport electrique

      Probleme_base& pb_phi = ref_cast(Probleme_base,interprete().objet(nom_pb_phi_));
      ch_phi_ = pb_phi.get_champ(nom_champ_phi_);
      assert(ch_phi_.valeur().que_suis_je().find("P1NC") !=-1);
      ch_psi_ = equation().inconnue();
      assert(ch_psi_.valeur().que_suis_je().find("P1NC") !=-1);
      psi_.ref(ch_psi_.valeur().valeurs());
    }
  else
    {
      Cerr << "On ne connait pas le probleme actuel & couple" << finl;
      exit();
    }

  // discretiser les champs
  discretiser(equation().discretisation());

  return is;
}

void Source_Term_pemfc_transport_ionique::set_param(Param& param)
{
  //param.ajouter("nom_domaine", &nom_domaine_, Param::REQUIRED);		// DVQ: pas necessaire car associe a l'equation -> get domain
  param.ajouter("nom_ssz_CLc", &nom_ssz_CLc_, Param::REQUIRED);
  param.ajouter("nom_ssz_CLa", &nom_ssz_CLa_, Param::REQUIRED);
  param.ajouter("nom_pb_psi", &nom_pb_psi_, Param::OPTIONAL);
  param.ajouter("nom_champ_psi", &nom_champ_psi_, Param::OPTIONAL);
  param.ajouter("nom_pb_phi", &nom_pb_phi_, Param::OPTIONAL);
  param.ajouter("nom_champ_phi", &nom_champ_phi_, Param::OPTIONAL);
  param.ajouter("nom_pb_T", &nom_pb_T_, Param::REQUIRED);
  param.ajouter("nom_champ_T", &nom_champ_T_, Param::REQUIRED);
  param.ajouter("nom_pb_C_H2", &nom_pb_C_H2_, Param::REQUIRED);
  param.ajouter("nom_champ_C_H2", &nom_champ_C_H2_, Param::REQUIRED);
  param.ajouter("nom_pb_C_O2", &nom_pb_C_O2_, Param::REQUIRED);
  param.ajouter("nom_champ_C_O2", &nom_champ_C_O2_, Param::REQUIRED);
  param.ajouter("nom_pb_C_H2O", &nom_pb_C_H20_, Param::REQUIRED);
  param.ajouter("nom_champ_C_H2O", &nom_champ_C_H20_, Param::REQUIRED);
}

void Source_Term_pemfc_transport_ionique::associer_pb(const Probleme_base& pb)
{
  // recuperer le milieu_base -> nothing to do
  Cerr << " Source_Term_pemfc_transport_ionique::associer_pb " << finl ;
  assert(pb.que_suis_je() == "Pb_Conduction");
//  int ok = 0;
//  const Equation_base& eqn = pb.equation(0);
//  assert(eqn.que_suis_je() == "Conduction");
//  if  (eqn.que_suis_je() == "Conduction")
//    {
//      associer_zones(eqn.zone_dis(),eqn.zone_Cl_dis());
//      ok = 1;
//    }
//  if (!ok)
//    {
//      Cerr << "Erreur TRUST dans Source_Term_pemfc_transport_ionique::associer_pb()" << finl;
//      Cerr << "On ne trouve pas d'equation de conduction dans le probleme" << finl;
//      exit();
//    }

}

void Source_Term_pemfc_transport_ionique::completer()
{
  // get all references to the coupling fields
  Source_base::completer();

//  int dim = 1;		// all fields are scalar
//  T_.resize(0, dim);
//  a_H2_.resize(0, dim);
//  a_O2_.resize(0, dim);
//  a_H2O_.resize(0, dim);
//  a_H_.resize(0, dim);
  if(nom_pb_psi_ != "??") 		// transport ionique couple avec psi
    {
//      psi_.resize(0, dim);
      la_zone_VEF.valeur().creer_tableau_faces(psi_);
      phi_.ref(ch_phi_.valeur().valeurs());			// phi
    }
  if(nom_pb_phi_ != "??") 		// transport electrique couplage avec phi
    {
//      phi_.resize(0, dim);
      la_zone_VEF.valeur().creer_tableau_faces(phi_);
      psi_.ref(ch_psi_.valeur().valeurs());			// psi
    }

  la_zone_VEF.valeur().creer_tableau_faces(T_);
  la_zone_VEF.valeur().creer_tableau_faces(C_H2_);
  la_zone_VEF.valeur().creer_tableau_faces(C_O2_);
  la_zone_VEF.valeur().creer_tableau_faces(C_H2O_);
  //la_zone_VEF.valeur().creer_tableau_faces(a_H_);

  // ATTENTION: needs to verifying carefully the initial value -> function with exp() et pow() easy to cause a floating point error

  //T_ = T_ref;
  //C_H2_ = a_lim;
  //C_O2_ = a_lim;
  //C_H2O_ = a_lim;
  //a_H_ = a_lim;
}

void Source_Term_pemfc_transport_ionique::associer_zones(
  const Zone_dis& zone_dis, const Zone_Cl_dis& zcl_dis)
{
  la_zone_VEF = ref_cast(Zone_VEF,zone_dis.valeur());
  la_zcl_VEF = ref_cast(Zone_Cl_VEF,zcl_dis.valeur());
  remplir_volumes();
}

void Source_Term_pemfc_transport_ionique::discretiser(const Discretisation_base& dis)
{
  dis.discretiser_champ("temperature",equation().zone_dis().valeur(),"Erev","unit", 1 , 0. , ch_erev_);
  champs_compris_.ajoute_champ(ch_erev_);
  dis.discretiser_champ("temperature",equation().zone_dis().valeur(),"eta","unit", 1 , 0. , ch_eta_);
  champs_compris_.ajoute_champ(ch_eta_);
  dis.discretiser_champ("temperature",equation().zone_dis().valeur(),"io","unit", 1 , 0. , ch_io_);
  champs_compris_.ajoute_champ(ch_io_);
  dis.discretiser_champ("temperature",equation().zone_dis().valeur(),"ir","unit", 1 , 0. , ch_ir_);
  champs_compris_.ajoute_champ(ch_ir_);
  dis.discretiser_champ("temperature",equation().zone_dis().valeur(),"ip","unit", 1 , 0. , ch_ip_);
  champs_compris_.ajoute_champ(ch_ip_);
  dis.discretiser_champ("temperature",equation().zone_dis().valeur(),"Q_reac","unit", 1 , 0. , ch_q_reac_);
  champs_compris_.ajoute_champ(ch_q_reac_);
  dis.discretiser_champ("temperature",equation().zone_dis().valeur(),"Q_perm","unit", 1 , 0. , ch_q_perm_);
  champs_compris_.ajoute_champ(ch_q_perm_);
}

void Source_Term_pemfc_transport_ionique::mettre_a_jour(double temps)
{
  // mettre a jour 4 tableaux de valeurs du champ couple
  ch_T_.valeur().mettre_a_jour(temps);
  ch_C_H2_.valeur().mettre_a_jour(temps);
  ch_C_O2_.valeur().mettre_a_jour(temps);
  ch_C_H20_.valeur().mettre_a_jour(temps);
  ch_psi_.valeur().mettre_a_jour(temps);
  ch_phi_.valeur().mettre_a_jour(temps);

  const DoubleTab& xv=la_zone_VEF.valeur().xv(); // centre de gravite des faces pour P1NC

  // interpolation vers P1NC
  ch_T_.valeur().valeur_aux( xv, T_ );
  ch_C_H2_.valeur().valeur_aux( xv, C_H2_ );		// C_H2
  ch_C_O2_.valeur().valeur_aux( xv, C_O2_ );		// C_O2
  ch_C_H20_.valeur().valeur_aux( xv, C_H2O_ );		// C_H2O

  if(nom_pb_psi_ != "??") 		// transport ionique couple avec psi
    {
      ch_psi_.valeur().valeur_aux( xv, psi_ );		// psi interpole
      phi_.ref(ch_phi_.valeur().valeurs());			// phi
    }
  if(nom_pb_phi_ != "??") 		// transport electrique couplage avec phi
    {
      ch_phi_.valeur().valeur_aux( xv, phi_ );		// phi interpole
      psi_.ref(ch_psi_.valeur().valeurs());			// phi
    }

  /*
  int nb_faces = la_zone_VEF.valeur().nb_faces();
  // convertir le champ C -> le champ activite
  for (int face = 0; face < nb_faces; ++face)
    {
      double Tf = T_(face);
      a_H2_(face) /= f_Henry_H2(Tf)*P_ref;			// a_H2 = C_H2 / (Henry_H2 * Pref)
      a_O2_(face) /= f_Henry_O2(Tf)*P_ref;			// a_02 = C_02 / (Henry_O2 * Pref)
      double ld = a_H2O_(face)/C_SO3;				// ld = C_H20 / C_SO3
      a_H2O_(face) = f_lambda_inv(ld);				// a_H20 = f_lambda_inv(ld) -> need testing -> dvq: tested OK
      a_H_(face) = f_lambda(1.) / (ld);				// ATTENTION divise by zero
      //a_H_(face) = 1.;								// approximation 1.
    }
  */
  // mettre a jour des champs compris
  DoubleTab& val_eta = ch_eta_.valeurs();
  DoubleTab& val_erev_ = ch_erev_.valeurs();
  DoubleTab& val_i0 = ch_io_.valeurs();
  DoubleTab& val_ir = ch_ir_.valeurs();
  DoubleTab& val_ip = ch_ip_.valeurs();
  DoubleTab& val_q_reac = ch_q_reac_.valeurs();
  DoubleTab& val_q_perm = ch_q_perm_.valeurs();

  IntTab faces_ssz;	// faces belong to the sous_zone -> flag = 1, if not, flag = 0
  la_zone_VEF.valeur().creer_tableau_faces(faces_ssz);
  faces_ssz = 0;		// init with no flag (all faces are unchecked)

  for (int poly = 0; poly < CL_a_.valeur().nb_elem_tot(); poly++)
    {
      int elem = CL_a_.valeur()(poly);
      for (int f = 0; f < la_zone_VEF.valeur().zone().nb_faces_elem(0); ++f)
        {
          int face = la_zone_VEF.valeur().elem_faces(elem, f);
          if(!faces_ssz(face))
            {
              double T = T_(face);
              double a_H2 = C_H2_(face)/(f_Henry_H2(T)*P_ref);
              double ld = C_H2O_(face)/C_SO3;
              if(ld <= 0.)
                {
                  Cerr << "ld <= 0. -> set to a_lim" << finl;
                  ld = a_lim;
                }
              double a_H = f_lambda(1.) / (ld);										// ATTENTION divise by zero
              double phi = phi_(face);
              double psi = psi_(face);

              double erev = eval_erev_anode(T,a_H2, a_H);
              double eta  = eval_eta(psi, phi , erev);
              double i0   = eval_i0_anode(T, a_H2, a_H);
              double ir	  = eval_ir_anode(i0, eta, T);
              double ip   = 0;														// TO-DO: need to evaluating
              double q_reac = eval_q_reac_anode(psi, phi, ir);
              double q_perm = eval_q_perm_anode(ip);

              val_erev_(face) 	= erev;
              val_eta(face) 	= eta;
              val_i0(face)		= i0;
              val_ir(face) 		= ir;
              val_ip(face)		= ip;
              val_q_reac		= q_reac;
              val_q_perm		= q_perm;
              faces_ssz(face) = 1;		// marquer comme deja traite
            }
        }
    }
  faces_ssz = 0;		// init with no flag (all faces are unchecked)
  for (int poly = 0; poly < CL_c_.valeur().nb_elem_tot(); poly++)
    {
      int elem = CL_c_.valeur()(poly);
      for (int f = 0; f < la_zone_VEF.valeur().zone().nb_faces_elem(0); ++f)
        {
          int face = la_zone_VEF.valeur().elem_faces(elem, f);
          if(!faces_ssz(face))
            {
              double T = T_(face);
              double a_O2 = C_O2_(face)/(f_Henry_O2(T)*P_ref);
              double ld = C_H2O_(face)/C_SO3;
              if(ld <= 0.)
                {
                  Cerr << "ld <= 0. -> set to a_lim" << finl;
                  ld = a_lim;
                }
              double a_H2O = f_lambda_inv(ld);				// a_H20 = f_lambda_inv(ld) -> need testing -> dvq: tested OK
              double a_H = f_lambda(1.) / (ld);										// ATTENTION divise by zero
              double phi = phi_(face);
              double psi = psi_(face);

              double erev = eval_erev_cathode(T, a_O2 ,a_H2O, a_H);
              double eta  = eval_eta(psi, phi, erev);
              double i0   = eval_i0_cathode(T, a_O2, a_H2O, a_H);
              double ir	  = eval_ir_cathode(i0, eta, T);
              double ip   = 0;														// TO-DO: need to evaluating
              double q_reac = eval_q_reac_cathode(psi, phi, ir);
              double q_perm = eval_q_perm_cathode(ip);
              val_erev_(face) 	= erev;
              val_eta(face) 	= eta;
              val_i0(face)		= i0;
              val_ir(face) 		= ir;
              val_ip(face)		= ip;
              val_q_reac		= q_reac;
              val_q_perm		= q_perm;
              faces_ssz(face) = 1;		// marquer comme deja traite
            }
        }
    }





  // DEBUG
  //Cerr << "champ_erev" << val_erev_ << finl;
  //Cerr << "champ_psi" << psi_ << finl;
  //Cerr << "champ_phi" << phi_ << finl;
  //Cerr << "champ_eta" << val_eta << finl;
}

// -dG/(n*F)-R*T*log(a_X2_lim^nu_H2*a_H^nu_H_a)/(n*F)
double Source_Term_pemfc_transport_ionique::eval_erev_anode(double T, double a_H2, double a_H)
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
double Source_Term_pemfc_transport_ionique::eval_erev_cathode(double T, double a_O2, double a_H2O, double a_H)
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

double Source_Term_pemfc_transport_ionique::eval_eta(double psi, double phi, double erev)
{
  return psi - phi - erev;
}

// n*F*kox0^(1-alpha)*kred0^alpha
double Source_Term_pemfc_transport_ionique::eval_i0_anode(double T, double a_H2, double a_H)
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

double Source_Term_pemfc_transport_ionique::eval_i0_cathode(double T, double a_O2, double a_H2O, double a_H)
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

double Source_Term_pemfc_transport_ionique::eval_ir_anode(double io, double eta, double T)
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

double Source_Term_pemfc_transport_ionique::eval_ir_cathode(double io, double eta, double T)
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

// source = ie = ir+ip si ionique, source = -ie = -(ir+ip) si electrique
DoubleTab& Source_Term_pemfc_transport_ionique::ajouter(DoubleTab& resu) const
{
  assert(resu.dimension(0)==volumes_.size());
  assert(resu.dimension(0)==T_.size());
  assert(resu.dimension(0)==C_H2_.size());
  assert(resu.dimension(0)==C_O2_.size());
  assert(resu.dimension(0)==C_H2O_.size());
  assert(resu.dimension(0)==psi_.size());
  assert(resu.dimension(0)==phi_.size());

  int signe = (nom_champ_psi_ != "??")?1:-1;		// ionique (phi) = 1, electrique (psi) = -1

  int nb_faces = la_zone_VEF.valeur().zone().nb_faces_elem(0);
  for (int face = 0; face < nb_faces; ++face)
    {
      // champ ir doit etre mettre_a_jour() avant ajouter()
      resu(face) += signe * (ch_ir_.valeurs()(face) + ch_ip_.valeurs()(face)) * volumes_(face);				// -> NEEDS VERIFYING
    }
  return resu;
}

DoubleTab& Source_Term_pemfc_transport_ionique::calculer(DoubleTab& resu) const
{
  resu = 0;
  ajouter(resu);
  return resu;
}

double Source_Term_pemfc_transport_ionique::eval_q_reac_anode(double psi, double phi, double ir)
{
  return (psi - phi + dH0_a / (n_a*F))*ir;
}

double Source_Term_pemfc_transport_ionique::eval_q_reac_cathode(double psi, double phi, double ir)
{
  return (psi - phi + dH0_c / (n_c*F))*ir;
}

double Source_Term_pemfc_transport_ionique::eval_q_perm_anode(double ip)
{
  return -(dH0_a+dH0_c)/(n_a*F)*ip;
}

double Source_Term_pemfc_transport_ionique::eval_q_perm_cathode(double ip)
{
  return -(dH0_a+dH0_c)/(n_c*F)*ip;
}

void Source_Term_pemfc_transport_ionique::remplir_volumes()
{
  volumes_.ref(la_zone_VEF.valeur().volumes_entrelaces());
}
