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
// File      : Loi_Fermeture_transport_ionique.h
// Directory : $PEMFC_ROOT/src/Loi_Fermeture
//
/////////////////////////////////////////////////////////////////////////////

#ifndef Loi_Fermeture_transport_ionique_included
#define Loi_Fermeture_transport_ionique_included

#include <Loi_Fermeture_base.h>
#include <Param.h>
#include <Discretisation_base.h>
#include <Champ_Fonc.h>
#include <Equation_base.h>
#include <Ref_Champ_base.h>
#include <Ref_Champ_Inc.h>
#include <Ref_Operateur_base.h>
#include <DoubleTab.h>
#include <Ref_Sous_Zone.h>
#include <Champ_Don.h>


/////////////////////////////////////////////////////////////////////////////
//
// .DESCRIPTION : class Loi_Fermeture_transport_ionique
// cette clase permet de creer / mettre a jour les champs suivants partout:
// - champ scalaire de coefficient de conductivite inonique kappa
//   qui est dependant a la temperature T et l'absoption de l'eau C_H20:
//   kappa = f_kappa(T,ld) avec ld = C_H20 / C_SO3
// - champ de courant ionique Ii = - kappa*grad(phi)
//   avec phi est le champ inconnu de l'equation de transport ionique
// Les champs electro-chimique suivants sont existants dans seul CLa et CLc:
// + champ de potentiel standard E_rev
// + eta 	psi - phi - E_rev
// + exchange current density i0 different selon CA anode ou cathode
// + le courant de reaction total ir a la surface du catalyseur, <0 a la cathode
// + le courant de permeation ip
// + Chaleur liee a la reaction Q_reac
// + Chaleur liee au courant de permeation Q_perm
// <Description of class Loi_Fermeture_transport_ionique>
//
/////////////////////////////////////////////////////////////////////////////

class Loi_Fermeture_transport_ionique : public Loi_Fermeture_base
{

  Declare_instanciable( Loi_Fermeture_transport_ionique ) ;

public :
  void discretiser(const Discretisation_base& );
  void completer();
  void mettre_a_jour(double temps);
  void set_param(Param& param);
  void compute_ir_DirDpsi(const DoubleTab& psi, const DoubleTab& phi, DoubleTab& ir, DoubleTab& DirDpsi);
  const ArrOfInt& face_index_psi() const
  {
    return face_index_psi_;
  }
  const ArrOfInt& face_index_phi() const
  {
    return face_index_phi_;
  }
  const ArrOfInt& face_index_T() const
  {
    return face_index_T_;
  }
  const ArrOfInt& face_is_anode() const
  {
    return face_is_anode_;
  }
protected :
  void fetch_field(const Nom& pbName, const Nom& fieldName, REF(Champ_base) & refField);

  Champ_Fonc ch_kappa_;	// P0 champ de conductivite (scalaire)
  Champ_Fonc ch_Ii_;	// P0 champ de courant (vectoriel)
  Champ_Fonc ch_erev_;	// P0 reference potential
  Champ_Fonc ch_eta_;	// P0 eta
  Champ_Fonc ch_io_;	// P0 exchange current density
  Champ_Fonc ch_ir_;	// P0 champ electrochemique courant ir = io.gamma_CL.[exp()-exp()]
  Champ_Fonc ch_ip_;  	// P0 champ electrochemique courant ip = N_H2.2.F/e_CL
  Champ_Fonc ch_q_reac_;// P0 chaleur
  Champ_Fonc ch_q_perm_;// P0 chaleur
  // This one is optionaly read in the data file:
  Champ_Don  Erev_field_override_; // if provided, overrides the value computed from activities
  Champ_Don  i0_field_override_;
  // Fields discretized like the unknowns:
  Champ_Fonc Erev_field_;
  Champ_Fonc i0_field_;
  Champ_Fonc ir_field_;

  Champ_Fonc ch_dir_dphi_;	// P0 derivee ir par rapport de phi
  //Champ_Fonc ch_dir_dpsi_;

  Champ_Fonc ch_DirDcH2_;	// P0 derivatives for implicitation
  Champ_Fonc ch_DirDcH2O_;	// P0 derivatives for implicitation
  Champ_Fonc ch_DirDcO2_;	// P0 derivatives for implicitation

  Champ_Fonc ch_DQreacDT_;

  Champ_Don temperature_; // temperature if provided
  Champ_Don Co_;		// dissolved concentration of O2, if provided
  Champ_Don Ch_;		// dissolved concentration of H2, if provided
  Champ_Don Ce_;		// dissolved concentration of H2O, if provided

  // pour readOn
  Nom nom_ssz_CLa_;
  Nom nom_ssz_CLc_;
  Nom nom_pb_phi_;			// transport ionique
  Nom nom_champ_phi_;		// champ potentiel ionique
  Nom nom_pb_psi_;			// transport electrique
  Nom nom_champ_psi_;		// champ potentiel electrique
  Nom nom_pb_T_;			// pb conduction de chaleur
  Nom nom_champ_T_;			// defaut temperature
  Nom nom_pb_C_H2_;			// dissolve H2
  Nom nom_champ_C_H2_;		// champ C_H2
  Nom nom_pb_C_O2_;			// dissolve O2
  Nom nom_champ_C_O2_;		// champ C_O2
  Nom nom_pb_C_H2O_;		// dissolve H2O
  Nom nom_champ_C_H2O_;		// champ C_H2O

  //double T_0_;				// dans le cas T constant
  //double C_SO3_;				//
  Champ_Don por_naf_;			// porosite de Nafion
  Champ_Don eps_naf_;			// ionomer proportionnel de Nafion
  Champ_Don tor_naf_;			// tortuosite de Nafion
  REF(Equation_base) ref_equation_;
  REF(Sous_Zone) CL_a_;						// sous_zone anode
  REF(Sous_Zone) CL_c_;						// sous_zone cathode
  const Equation_base& equation() const
  {
    return ref_equation_.valeur();
  } ;

  REF(Champ_base) ch_T_;	// champ temperature (importe depuis pb_conduction_T)
  REF(Champ_base) ch_C_H2_;	// champ concentration absorbe H2 dans Nafion (importe depuis pb_dissolve_H2)
  REF(Champ_base) ch_C_O2_;	// champ concentration absorbe O2 dans Nafion (importe depuis pb_dissolve_O2)
  REF(Champ_base) ch_C_H2O_;// champ concentration absorbe H20 dans Nafion (importe depuis pb_dissolve_H2O)
  REF(Champ_base) ch_phi_;	// champ de potentiel ionique
  REF(Champ_base) ch_psi_;	// champ de potentiel electrique

  DoubleTab T_;			// P0 pour VDF, P1NC pour VEF
  DoubleTab C_H2_;		// P0 pour VDF, P1NC pour VEF
  DoubleTab C_O2_;		// P0 pour VDF, P1NC pour VEF
  DoubleTab C_H2O_;		// P0 pour VDF, P1NC pour VEF
  DoubleTab psi_;		// P0 pour VDF, P1NC pour VEF
  DoubleTab phi_;		// P0 pour VDF, P1NC pour VEF

  double newton_threshold_ ;
  int newton_max_iter_ ;
  double minimal_perturbation_value_;
  double relative_perturbation_for_derivatives_ ;

  inline double f_kappa(const double& T, const double& C, const double& por, const double& eps, const double& tor) const;
  inline  double f_Henry_H2(const double& T) const;		// Constante de Henry H2 dans Nafion
  inline  double f_Henry_O2(const double& T) const ;		// Constante de Henry O2 dans Nafion
  inline  double f_Henry_N2(const double& T) const ;		// Constante de Henry N2 dans Nafion
  inline  double f_Psat(const double& T) const ;			// Pression vapeur saturante
  inline  double f_lambda(const double& a) const ;		// lambda en fonction de l'activite
  inline  double f_lambda_inv(const double& ld) const ;	// activite en fonction de lambda

  double eval_erev_anode(double T, double a_H2, double a_H) const;
  double eval_erev_cathode(double T, double a_O2, double a_H2O, double a_H) const;
  double eval_eta(double psi, double phi, double erev) const;
  double eval_i0_anode(double T, double a_H2, double a_H) const;
  double eval_i0_cathode(double T, double a_O2, double a_H2O, double a_H) const;
  double eval_ir_anode(double io, double eta, double T) const;
  double eval_ir_cathode(double io, double eta, double T) const;
  double eval_q_reac_anode(double psi, double phi, double ir) const;
  double eval_q_reac_cathode(double psi, double phi, double ir) const;
  double eval_q_perm_anode(double ip) const;
  double eval_q_perm_cathode(double ip) const;

  // dvq 28/05/19: ajouter le calcul de derive de ir par rapport phi et psi
  double eval_dirdphi_anode(double io, double eta, double T) const;
  double eval_dirdpsi_anode(double io, double eta, double T) const;
  double eval_dirdphi_cathode(double io, double eta, double T) const;
  double eval_dirdpsi_cathode(double io, double eta, double T) const;

  double compute_perturbation( const double& f, const Nom& info ) const;
  void eval_derivatives_anode_H2(const int& elem, const double& T,
                                 const double& C_H2, const double& a_H, double& Dir ) const;
  void eval_derivatives_cathode_O2(const int& elem, const double& T,
                                   const double& C_O2, const double& a_H2O, const double& a_H, double& Dir ) const;
  void eval_derivatives_cathode_H2O(const int& elem, const double& T,
                                    const double& C_O2, const double& C_H2O, const double& a_O2, const double& a_H, double& Dir ) const;

  // dvq 18/06/19 ajouter le calcul de derivee de ir par rapport T
  double eval_derevdT_anode(double T, double a_H2, double a_H) const;
  double eval_derevdT_cathode(double T, double a_O2, double a_H2O, double a_H) const;
  double eval_dirdT_anode(const double& io, const double& eta, const double& T, const double& dErevdT) const;
  double eval_dirdT_cathode(const double& io, const double& eta, const double& T, const double& dErevdT) const;

  // For each common unknown between the two domains phi and psi (1 common unknown = 1 line):
  ArrOfInt face_index_psi_; // index of face in psi domain
  ArrOfInt face_index_phi_; // index of face in phi domain
  ArrOfInt face_index_T_; // index of face in temperature domain
  ArrOfInt face_is_anode_; // 1 if yes, 0 if no.
  // Control volume associated with the unknowns common to phi/psi.
  // The control volume only counts the part in the intersection
  // (eg: for faces on the boundary of phi domain, only the volume inside)
  // Number of lines in the array = number of lines in face_index_psi_
  //DoubleTab psi_phi_control_volumes_;

};

const double a_lim = 1e-3;				// Activit minimale en dessous de laquelle son effet est linearise
const double v0_max = 100;				// Vitesse maximale de reaction (par unite de surface de Pt) [mol/m^2/s]
const int n_a = 2;						// n_c Nombre d''electrons echanges a l''anode'
const int n_c = 2;						// n_c Nombre d''electrons echanges a la cathode
const double kB = 1.38064852e-23;	 	// Boltzman constant [J/K]
const double s0 = 6.41e-20;				// average Pt surface per reaction site [m2]
const double gamma_CL_a = 1.67e7;			// specific surface [m2/m3]
const double gamma_CL_c = 2.5e7;
const double NA = 6.022140857e23;		// Avogadro constant [/mol]
const double h  = 6.62607015e-34;		// Planck constant [Js]
const double R = 8.314;					// gas constant J/mol/K
const double F = 96500;					// Faraday constant C/mol
const double nu_H2 = 1.;		// 0.5*n_a Coefficient stoechiometrique de la reaction pour l''hydrogene
const double nu_H_a = -2.;		// -n_a Coefficient stoechiometrique de la reaction pour les protons a l''anode
const double nu_O2 = -0.5;		// -0.25*n_c Coefficient stoechiometrique de la reaction pour l''oxygene
const double nu_H2O = 1;		// 0.5*n_c Coefficient stoechiometrique de la reaction pour l''eau
const double nu_H_c = -2.;		// -n_c Coefficient stoechiometrique de la reaction pour les protons a la cathode

//const double Cdl_a = 10.;			// [F/m^2], Capacite de double couche
//const double Cdl_c = 10.;			// [F/m^2], Capacite de double couche

const double alpha_a = 0.5;			// Coefficient de transfert de charge a l''anode
const double alpha_c = 0.5; // 0.216;		// Nombre d''electrons echanges a la cathode
const double dHox0_a = 25.e3;		// [J/mol] Enthalpie de formation du complexe active (dans le sens de l''oxidation)
const double dSox0_a = -172;		// [J/mol/K] Entropie de formation du complexe active (dans le sens de l''oxidation)
const double dHox0_c = 1.679e5;		// [J/mol] Enthalpie de formation du complexe active (dans le sens de l''oxidation)
const double dSox0_c = -205.6;		// [J/mol/K] Entropie de formation du complexe active (dans le sens de l''oxidation)

const double dH0_a = 0;
const double dS0_a = -0.104;
const double dH0_c = -2.858e5;
const double dS0_c = -163.8;

const double C_SO3 = 2036.; 	// [mol/m^3], Concentration en sites sulfones dans le Nafion

const double M_H2O = 18e-3; 	//	[kg/mol]
const double M_O2 = 32e-3; 		// [kg/mol]
const double M_H2 = 2e-3; 		// [kg/mol]
const double M_N2 = 28e-3; 		// [kg/mol]
const double rhol = 1e3; 		// [kg/m^3]
const double P_ref = 1.013e5; 	// [Pa]
const double T_ref = 298; 		// [K]
const double Cp_O2 = 29.75; 	// [J/mol/K]
const double Cp_H2 = 28.86; 	// [J/mol/K]
const double Cp_N2 = 30.2; 		// [J/mol/K]
const double Cp_vap = 34.474; 	// [J/mol/K]
const double Cp_liq = 75.38;	 // [J/mol/K]

inline double Loi_Fermeture_transport_ionique::f_kappa(const double& T, const double& C, const double& por, const double& eps, const double& tor) const
{
  double ld = C / C_SO3;
  double ad = 0.5139*ld-0.326;
  // double ad_lim = max(ad,1.e-3);
  // double k0 = exp(1268*(1./303.-1./T))*ad_lim;
  double k0 = exp(1268*(1./303.-1./T))*ad;
  return k0*(1-por)*eps/(tor*tor);
}

inline double Loi_Fermeture_transport_ionique::f_Henry_H2(const double& T) const
{
  return (1./1.01325e5)*exp(9.05e3/(R*T));
}

inline double Loi_Fermeture_transport_ionique::f_Henry_O2(const double& T) const
{
  return (1/1.01325e5)*exp(5.88e3/(R*T));
}

inline double Loi_Fermeture_transport_ionique::f_Henry_N2(const double& T) const
{
  return 6.4e-6*exp(1300*(1/T-1/298.15));
}

inline double Loi_Fermeture_transport_ionique::f_lambda(const double& a) const
{
  return 0.043 + 17.81*a - 39.85*a*a + 36*a*a*a;
}

inline double Loi_Fermeture_transport_ionique::f_lambda_inv(const double& ld) const
{
  // using a newton raphson iteration
  int it ;
  double a = -0.043/17.81;

  double num = f_lambda(a) - ld;
  double den = 17.81 - 2*39.85*a + 3*36*a*a;
  double da = num/den;
  a -= da;
  it=0;
  while (da/a>newton_threshold_ && it<newton_max_iter_)
    {
      num=f_lambda(a) - ld;
      den=17.81 - 2*39.85*a + 3*36*a*a;
      da = num/den;
      a -= da;
      it++;
    }

  if (da/a>newton_threshold_)
    {
      Cerr<<"Loi_Fermeture_transport_ionique::f_lambda_inv resolution Newton fail "<<finl;
      Cerr<<"ld="<< ld <<" da="<<da<<" a="<<a<<finl;
      exit();
    }
  return a;
}

inline double Loi_Fermeture_transport_ionique::f_Psat(const double& T) const
{
  return exp(23.1961-3816.44/(T-46.13));
}

#endif /* Loi_Fermeture_transport_ionique_included */
