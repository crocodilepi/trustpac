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
// File      : Source_Term_pemfc_transport_ionique.h
// Directory : $PEMFC_ROOT/src/Sources
//
/////////////////////////////////////////////////////////////////////////////

#ifndef Source_Term_pemfc_transport_ionique_included
#define Source_Term_pemfc_transport_ionique_included

#include <Source_base.h>
#include <Ref_Zone_VEF.h>
#include <Ref_Zone_Cl_VEF.h>
#include <Ref_Domaine.h>
#include <Ref_Sous_Zone.h>
#include <Param.h>
#include <DoubleTab.h>
#include <Champ_Fonc.h>
#include <Discretisation_base.h>

/////////////////////////////////////////////////////////////////////////////
//
// .DESCRIPTION : class Source_Term_pemfc_transport_ionique
//
// <Description of class Source_Term_pemfc_transport_ionique>
//
/////////////////////////////////////////////////////////////////////////////

class Source_Term_pemfc_transport_ionique : public Source_base
{

  Declare_instanciable( Source_Term_pemfc_transport_ionique ) ;

public :
  DoubleTab& ajouter(DoubleTab& ) const;
  DoubleTab& calculer(DoubleTab& ) const;
  void mettre_a_jour(double temps);
  void associer_zones(const Zone_dis& ,const Zone_Cl_dis& );
  void associer_pb(const Probleme_base& );
  void completer();
  void set_param(Param& param);
  void discretiser(const Discretisation_base& dis);
protected :

  double eval_erev_anode(double T, double a_H2, double a_H);
  double eval_erev_cathode(double T, double a_O2, double a_H2O, double a_H);
  double eval_eta(double psi, double phi, double erev);
  double eval_i0_anode(double T, double a_H2, double a_H);
  double eval_i0_cathode(double T, double a_O2, double a_H2O, double a_H);
  double eval_ir_anode(double io, double eta, double T);
  double eval_ir_cathode(double io, double eta, double T);
  double f_lambda(double a);

  void remplir_volumes();

  REF(Zone_VEF) la_zone_VEF;
  REF(Zone_Cl_VEF) la_zcl_VEF;
  DoubleVect volumes_;

  Nom nom_domaine_;			// pas necessaire car le terme source est associe a l'equation
  Nom nom_ssz_CLc_;			// le sous zone situant le terme source
  Nom nom_ssz_CLa_;			// le sous zone situant le terme source
  Nom nom_pb_psi_;			// transport electrique
  Nom nom_champ_psi_;			// champ potentiel electrique
  Nom nom_pb_T_;				// conduction de la chaleur
  Nom nom_champ_T_;			// champ temperature
  Nom nom_pb_dissolve_H2_;	// dissolve H2
  Nom nom_champ_C_H2_;		// champ C_H2 -> activite H2 a_H2 = C_H2 / (Henry_H2 * Pref)
  Nom nom_pb_dissolve_O2_;	// dissolve 02
  Nom nom_champ_C_O2_;		// champ C_O2 -> activite 02 a_02 = C_02 / (Henry_O2 * Pref)
  Nom nom_pb_dissolve_H20_;	// dissolve H2O
  Nom nom_champ_C_H20_;		// champ C_H20 -> activite H2O a_H2O = lambda / lambda_eq = C_H2O*R*T / P_sat

  REF(Domaine) dom_;		// domaine de transport ionique (CL_c + MB + CL_a)
  REF(Sous_Zone) CL_a_;	// sous_zone anode
  REF(Sous_Zone) CL_c_;	// sous_zone cathode
  REF(Champ_base) ch_T_;		// champ temperature (importe depuis pb_conduction_T)
  REF(Champ_base) ch_C_H2_;	// champ concentration absorbe H2 dans Nafion (importe depuis pb_dissolve_H2)
  REF(Champ_base) ch_C_O2_;	// champ concentration absorbe 02 dans Nafion (importe depuis pb_dissolve_02)
  REF(Champ_base) ch_C_H20_;	// champ concentration absorbe H20 dans Nafion (importe depuis pb_dissolve_H20)
  REF(Champ_base) ch_phi_;	// champ de potentiel ionique (champ inconnu)
  REF(Champ_base) ch_psi_;	// champ de potentiel electrique (importe depuis pb_transport_electrique)

  DoubleTab T_;			// P0 pour VDF, P1NC pour VEF
  DoubleTab a_H2_;		// P0 pour VDF, P1NC pour VEF
  DoubleTab a_O2_;		// P0 pour VDF, P1NC pour VEF
  DoubleTab a_H2O_;		// P0 pour VDF, P1NC pour VEF
  DoubleTab a_H_;		// P0 pour VDF, P1NC pour VEF
  DoubleTab psi_;		// P0 pour VDF, P1NC pour VEF
  DoubleTab phi_;		// P0 pour VDF, P1NC pour VEF

  Champ_Fonc ch_erev_;	// P0 pour VDF, P1NC pour VEF reference potential
  Champ_Fonc ch_eta_;		// P0 pour VDF, P1NC pour VEF eta
  Champ_Fonc ch_io_;		// P0 pour VDF, P1NC pour VEF exchange current density
  Champ_Fonc ch_ir_;		// P0 pour VDF, P1NC pour VEF champ electrochemique courant ir = io.gamma_CL.[exp()-exp()]
  //Champ_Fonc ch_jp_;  	// P0 pour VDF, P1NC pour VEF champ electrochemique courant ip = N_H2.2.F/e_CL

  // Constante de Henry H2 pour le Nafion
  double f_Henry_H2(double T);
  // Constante de Henry O2 pour le Nafion
  double f_Henry_O2(double T);
  // Constante de Henry N2 pour le Nafion
  double f_Henry_N2(double T);
  // Pression vapeur saturante
  double f_Psat(double T);
};

const double a_lim = 1e-3;				// Activit minimale en dessous de laquelle son effet est linearise
const double v0_max = 100;				// Vitesse maximale de reaction (par unite de surface de Pt) [mol/m^2/s]
const int n_a = 2;						// n_c Nombre d''electrons echanges a l''anode'
const int n_c = 2;						// n_c Nombre d''electrons echanges a la cathode
const double kB = 1.38064852e-23;	 	// Boltzman constant [J/K]
const double s0 = 6.41e-20;				// average Pt surface per reaction site [m2]
const double gamma_CL = 1.67e7;			// specific surface [m2/m3]
const double NA = 6.022140857e23;		// Avogadro constant [/mol]
const double h  = 6.62607015e-34;		// Planck constant [Js]
const double R = 8.314;					// gas constant J/mol/K
const double F = 96500;					// Faraday constant C/mol
const double nu_H2 = 1.;		// 0.5*n_a Coefficient stoechiometrique de la reaction pour l''hydrogene
const double nu_H_a = -2.;		// -n_a Coefficient stoechiometrique de la reaction pour les protons a l''anode
const double nu_O2 = -0.5;		// -0.25*n_c Coefficient stoechiometrique de la reaction pour l''oxygene
const double nu_H2O = 1;		// 0.5*n_c Coefficient stoechiometrique de la reaction pour l''eau
const double nu_H_c = -2.;		// -n_c Coefficient stoechiometrique de la reaction pour les protons a la cathode

const double Cdl_a = 10.;			// [F/m^2], Capacite de double couche
const double Cdl_c = 10.;			// [F/m^2], Capacite de double couche

const double alpha_a = 0.5;			// Coefficient de transfert de charge a l''anode
const double alpha_c = 0.216;		// Nombre d''electrons echanges a la cathode
const double dHox0_a = 24.36e3;		// [J/mol] Enthalpie de formation du complexe active (dans le sens de l''oxidation)
const double dSox0_a = -172.3;		// [J/mol/K] Entropie de formation du complexe active (dans le sens de l''oxidation)
const double dHox0_c = 81.52e3;		// [J/mol] Enthalpie de formation du complexe active (dans le sens de l''oxidation)
const double dSox0_c = -285.;		// [J/mol/K] Entropie de formation du complexe active (dans le sens de l''oxidation)
const double dH0_a = 25e3;
const double dS0_a = -172;
const double dH0_c = 167.9e3;
const double dS0_c = -205.6;

const double C_SO3 = 2036.; 	// [mol/m^3], Concentration en sites sulfones dans le Nafion

const double M_H2O = 18e-3; //	[kg/mol]
const double M_O2 = 32e-3; // [kg/mol]
const double M_H2 = 2e-3; // [kg/mol]
const double M_N2 = 28e-3; // [kg/mol]
const double rhol = 1e3; // [kg/m^3]
const double P_ref = 1.013e5; // [Pa]
const double T_ref = 298; // [K]
const double Cp_O2 = 29.75; // [J/mol/K]
const double Cp_H2 = 28.86; // [J/mol/K]
const double Cp_N2 = 30.2; // [J/mol/K]
const double Cp_vap = 34.474; // [J/mol/K]
const double Cp_liq = 75.38; // [J/mol/K]

inline double Source_Term_pemfc_transport_ionique::f_Henry_H2(double T)
{
  return (1./1.01325e5)*exp(9.05e3/(R*T));
}

inline double Source_Term_pemfc_transport_ionique::f_Henry_O2(double T)
{
  return (1/1.01325e5)*exp(5.88e3/(R*T));
}

inline double Source_Term_pemfc_transport_ionique::f_Henry_N2(double T)
{
  return 6.4e-6*exp(1300*(1/T-1/298.15));
}

inline double Source_Term_pemfc_transport_ionique::f_lambda(double a)
{
  return 0.043 + 17.81*a - 39.85*a*a + 36*a*a*a;
}

inline double Source_Term_pemfc_transport_ionique::f_Psat(double T)
{
  return exp(23.1961-3816.44/(T-46.13));
}

#endif /* Source_Term_pemfc_transport_ionique_included */
