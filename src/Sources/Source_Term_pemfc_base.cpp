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
// File      : Source_Term_pemfc_base.cpp
// Directory : $PEMFC_ROOT/src/Sources
//
/////////////////////////////////////////////////////////////////////////////

#include <Source_Term_pemfc_base.h>
#include <Probleme_base.h>
#include <Conduction.h>
#include <Param.h>
#include <Interprete.h>
#include <Matrice_Morse.h>
#include <Champ_base.h>
#include <Champ_Inc.h>
#include <Domaine.h>


Implemente_base( Source_Term_pemfc_base, "Source_Term_pemfc_base", Source_base ) ;

Sortie& Source_Term_pemfc_base::printOn( Sortie& os ) const
{
  Source_base::printOn( os );
  return os;
}

Entree& Source_Term_pemfc_base::readOn( Entree& is )
{
  //Source_base::readOn( is );
  Cerr << " Source_Term_pemfc_base::readOn " << finl  ;
  Param param(que_suis_je());
  param.ajouter("nom_espece",&nom_espece_,Param::REQUIRED);
  param.ajouter("nom_domaine",&nom_domaine_,Param::REQUIRED);		// pas necessaire ? dom_ = equation().problem().domaine() ?
  param.ajouter("nom_ssz",&nom_ssz_,Param::REQUIRED);
  param.ajouter("epsilon", &epsilon_, Param::REQUIRED);
  param.ajouter("epsilon_ionomer", &epsilon_ionomer_, Param::REQUIRED);
  param.ajouter("gamma_CL", &gamma_CL_, Param::REQUIRED);
  param.ajouter("champ_D", &nom_champ_D_, Param::REQUIRED);
  param.ajouter("champ_T", &nom_champ_T_, Param::REQUIRED);
  param.ajouter("champ_C", &nom_champ_C_, Param::REQUIRED);
  param.ajouter("champ_ci", &nom_champ_c_, Param::REQUIRED);
  param.lire_avec_accolades(is);

  Motcles nom_especes_compris_(5);
  nom_especes_compris_[0] = "H2";
  nom_especes_compris_[1] = "O2";
  nom_especes_compris_[2] = "H2O";
  nom_especes_compris_[3] = "N2";
  nom_especes_compris_[4] = "vap";
  if(nom_especes_compris_.search(nom_espece_) == -1)
    {
      Cerr <<" unknown species in the list "<<finl;
      Process::exit();
    }
  thickness_ionomer_ = (1-epsilon_)*epsilon_ionomer_ / gamma_CL_;
  dom_ = equation().probleme().domaine();
  ssz_ = dom_.valeur().ss_zone(nom_ssz_);
  return is;
}

double Source_Term_pemfc_base::eval_f(double diffu, double Ci, double ci, double T) const
{
  double R = 8.314;
  double Ceq;
  if (nom_espece_ == "H2O" || nom_espece_ == "vap")
    {
      double Psat = exp(23.1961-3816.44/(T-46.13));
      double activ = ci * R * T / Psat;
      double lambda_eq = 0.043+17.81*activ-39.85*activ*activ+36*activ*activ*activ;
      double CSO3 = 2036;
      Ceq = lambda_eq * CSO3;
    }
  else
    {
      double H = 0.; 			// H Henri constant (mol.m3/Pa)
      if (nom_espece_ == "H2")
        {
          H = 1. / 4.5e3;
        }
      else if (nom_espece_ == "02")
        {
          H = 1. / 1.33e5 * exp(666/T);
        }
      else if (nom_espece_ == "N2")
        {
          H = 6.4e-6 * exp(1300*(1./T - 1/298.15));
        }
      Ceq = ci * R * T * H;
    }
  return diffu * gamma_CL_ / thickness_ionomer_ * (Ceq - Ci);
}

double Source_Term_pemfc_base::eval_derivee_f(double diffu) const
{
  // expression_derivee_par_rapport_inconnue
  return (- diffu * gamma_CL_ / thickness_ionomer_);
}

DoubleTab& Source_Term_pemfc_base::calculer(DoubleTab& resu) const
{
  resu = 0;
  ajouter(resu);
  return resu;
}

void Source_Term_pemfc_base::mettre_a_jour(double temps)
{
  // update the involving fields -> to-do: nothing
}

void Source_Term_pemfc_base::completer()
{
  Source_base::completer();
}
