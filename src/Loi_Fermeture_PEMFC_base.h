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
// File      : Loi_Fermeture_PEMFC_base.h
// Directory : $PEMFC_ROOT/src
//
/////////////////////////////////////////////////////////////////////////////

#ifndef Loi_Fermeture_PEMFC_base_included
#define Loi_Fermeture_PEMFC_base_included

#include <Loi_Fermeture_base.h>
#include <Champ_Fonc.h>
#include <Champ_Don.h>
#include <Champ_Inc.h>
#include <Ref_Equation_base.h>
#include <Ref_Champ_base.h>
#include <DoubleTab.h>
#include <Zone_VF.h>

/////////////////////////////////////////////////////////////////////////////
//
// .DESCRIPTION : class Loi_Fermeture_PEMFC_base
//
// <Description of class Loi_Fermeture_PEMFC_base>
//
/////////////////////////////////////////////////////////////////////////////

class Loi_Fermeture_PEMFC_base : public Loi_Fermeture_base
{

  Declare_base( Loi_Fermeture_PEMFC_base ) ;

public :

  void discretiser(const Discretisation_base& );
  void mettre_a_jour(double);
  // void set_param(Param& param);
  void completer();
  void calculer_Ni_ud(const double& temps) ;		// pour maj Ni, ud, Um cas4
  void calculer_Ni_ug(const double& temps) ;		// pour maj Ni, ug cas2
  const Equation_base& equation() const
  {
    return ref_equation_.valeur();
  } ;
protected :
  // champ de diffusivite
  Champ_Fonc diffu_;

  // champ de flux diffusive = -nu.grad(ci)
  Champ_Fonc Ni_;

  // champ vitesse diffusive
  Champ_Fonc ud_;

  // champ vitesse massique
  Champ_Inc Um_;

  // Stefan Maxwell dans GDL (cas 2)
  // concentration total cg = sum(ci)
  Champ_Inc cg_gdl_;
  // vitesse diffusive ug = sum(Ni)/cg, Ni molar flux [mol/m2/s], cg total concentration [mol/m3], ug velocity m/s
  Champ_Inc ug_gdl_;
  // fraction molaire Xi = ci/cg [sans_unit]
  Champ_Inc Xi_gdl_;
  // flux_surfacique = integral(Ni, norma)
  Champ_Inc Ns_gdl_;

  // Stefan MAxwell dans CH (cas4)
  // concentration multi-composant ci = cg.Xi
  Champ_Inc ci_ch_;
  // flux_surfacique de concentration  Nx = integral((-nu.grad(ci)+ciUm),norma)
  Champ_Inc Ns_ch_;

  int is_cas4_;

  REF(Equation_base) ref_equation_;

  ArrOfDouble VGDL_;		// ?? pas compris

  // quang: 12/06/19 coupler avec le champ de temperature (optionnel), sinon, T=T_0=const
  Nom nom_pb_T_;
  REF(Champ_base) ch_T_;		// field temperature
  DoubleTab T_;
  // quang: 13/06 ajouter les champ donne porosite eps, tortuosite tau, rayon de pore Rp et permeabilite K
  Champ_Don ch_tor_;	// tortuosite
  Champ_Don ch_por_;	// porosite
  Champ_Don ch_K_;		// permeabilite
  Champ_Don ch_Rp_;		// pore radius
  ArrOfDouble Mi_;		// masse molaires des gaz melanges

  double eval_matrice_diffusion_cas2(const double& cX2,const double& cvap,const double& cN2, DoubleTab& Diff,int elem);
  double eval_matrice_diffusion_cas4(const double& cX2,const double& cvap,const double& cN2, DoubleTab& Diff,int elem);
  void bidouille_nu(DoubleTab& nu,const DoubleTab& inconnue_org,const Zone_VF& zone_VF);
  void discretisation_tools_cells_to_faces(const Champ_base& He,  Champ_base& Hf);

  double f_D_H2_vap(const double& T,const double& P) const;	// Diffusion H2/vapeur
  double f_D_O2_vap(const double& T,const double& P) const;	// Diffusion O2/vapeur
  double f_D_N2_vap(const double& T,const double& P) const;	// Diffusion N2/vapeur

  double f_D_H2_N2(const double& T,const double& P) const;	// Diffusion H2/N2
  double f_D_O2_N2(const double& T,const double& P) const;	// Diffusion O2/N2
  double f_D_H2_O2(const double& T,const double& P) const;	// Diffusion H2/O2
};

const double T_0 = 353.15;			// default temperature if not coupling
const double  R = 8.314;
const double  pi=3.1415926;					// math.pi

inline double Loi_Fermeture_PEMFC_base::f_D_H2_vap(const double& T,const double& Pg) const
{
  return 2.16e-5*pow( T ,2.334)/Pg;
}

inline double Loi_Fermeture_PEMFC_base::f_D_O2_vap(const double& T,const double& Pg) const
{
  return 4.26e-6*pow( T ,2.334)/Pg;
}

inline double Loi_Fermeture_PEMFC_base::f_D_N2_vap(const double& T,const double& Pg) const
{
  return 4.45e-6*pow( T ,2.334)/Pg;
}

inline double Loi_Fermeture_PEMFC_base::f_D_H2_N2(const double& T,const double& Pg) const
{
  return 2.46e-4*pow( T ,1.823)/Pg;
}

inline double Loi_Fermeture_PEMFC_base::f_D_O2_N2(const double& T,const double& Pg) const
{
  return 6.43e-5*pow( T ,1.823)/Pg;
}

inline double Loi_Fermeture_PEMFC_base::f_D_H2_O2(const double& T,const double& Pg) const
{
  return 1.60e-5*pow( T ,1.823)/Pg;
}



#endif /* Loi_Fermeture_PEMFC_base_included */
