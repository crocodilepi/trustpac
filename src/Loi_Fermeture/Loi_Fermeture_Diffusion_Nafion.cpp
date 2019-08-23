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
// File      : Loi_Fermeture_Diffusion_Nafion.cpp
// Directory : $PEMFC_ROOT/src/Loi_Fermeture
//
/////////////////////////////////////////////////////////////////////////////

#include <Loi_Fermeture_Diffusion_Nafion.h>
#include <Pb_Conduction.h>
#include <Param.h>
#include <Discretisation_base.h>
#include <Champ_base.h>
#include <Champ_P1NC.h>
#include <Zone_VF.h>
#include <Zone_Cl_VEF.h>
#include <Interprete.h>

Implemente_instanciable( Loi_Fermeture_Diffusion_Nafion, "Loi_Fermeture_Diffusion_Nafion", Loi_Fermeture_base ) ;

// XD Loi_Fermeture_Diffusion_Nafion Loi_Fermeture_base

const double C_SO3 = 2036.; 	// [mol/m^3], Concentration en sites sulfones dans le Nafion

Sortie& Loi_Fermeture_Diffusion_Nafion::printOn( Sortie& os ) const
{
  Loi_Fermeture_base::printOn( os );
  return os;
}

Entree& Loi_Fermeture_Diffusion_Nafion::readOn( Entree& is )
{
  Loi_Fermeture_base::readOn( is );
  Cerr << "Loi_Fermeture_Diffusion_Nafion::readOn " << nom_espece_ << finl;
  // check param input
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
  equation().zone_dis().zone().creer_tableau_elements(T_);
  equation().zone_dis().zone().creer_tableau_elements(C_);
  return is;
}

void Loi_Fermeture_Diffusion_Nafion::associer_pb_base(const Probleme_base& pb)
{
  assert(sub_type(Pb_Conduction, pb));
  Loi_Fermeture_base::associer_pb_base(pb);
}

void Loi_Fermeture_Diffusion_Nafion::discretiser(const Discretisation_base& dis)
{
  Loi_Fermeture_base::discretiser(dis);
  ref_equation_=mon_probleme().get_equation_by_name("Conduction");
  dis.discretiser_champ("champ_elem",equation().zone_dis().valeur(),"diffusion_nafion","unit", 1 ,0. , diffu_);
  champs_compris_.ajoute_champ(diffu_);

  dis.discretiser_champ("champ_elem",equation().zone_dis().valeur(),"D_i_naf","unit", 1 ,0. , D_i_naf_);
  champs_compris_.ajoute_champ(D_i_naf_);
  /*
    dis.discretiser_champ("champ_elem",equation().zone_dis().valeur(),"N_i_naf","unit", dimension ,0. , N_i_naf_);
    N_i_naf_->fixer_nature_du_champ(vectoriel);
    champs_compris_.ajoute_champ(N_i_naf_);*/
}

void Loi_Fermeture_Diffusion_Nafion::completer()
{
  Loi_Fermeture_base::completer();
  // get the reference to the coupling fields
  if(nom_champ_T_ != "??" && nom_pb_T_ != "??")
    {
      assert(!temperature_.non_nul());
      Probleme_base& pb_T = ref_cast(Probleme_base,interprete().objet(nom_pb_T_));
      ch_T_ = pb_T.get_champ(nom_champ_T_);
    }
  /*  if(nom_pb_phi_ != "??" && nom_champ_I_ != "??")
      {
        assert(nom_espece_ == "vap" || nom_espece_ == "H2O" );
        Probleme_base& pb_phi = ref_cast(Probleme_base,interprete().objet(nom_pb_phi_));
        ch_I_ = pb_phi.get_champ(nom_champ_I_);
      }*/
  ch_C_ = equation().inconnue();

  if(nom_champ_Ceq_ != "??")
    ch_Ceq_ = equation().probleme().get_champ(nom_champ_Ceq_);

  mettre_a_jour(0);
}

void Loi_Fermeture_Diffusion_Nafion::preparer_calcul()
{
  if(ch_Ceq_.non_nul())
    {
      Cerr << "Loi_Fermeture_Diffusion_Nafion::preparer_calcul " << equation().probleme().le_nom() << finl;
      DoubleTab& C_0 = ch_C_.valeur().valeurs();
      DoubleTab& Ceq_0 = ch_Ceq_.valeur().valeurs();
      assert(C_0.dimension(0)==Ceq_0.dimension(0));
      for (int i = 0; i < C_0.dimension(0); ++i)
        {
          C_0(i) = Ceq_0(i);
        }
      C_0.echange_espace_virtuel();
      Cerr << "initialize C_0 with Ceq_0 " << mp_min_vect(Ceq_0) << " " << mp_max_vect(Ceq_0) << finl;
      Loi_Fermeture_base::preparer_calcul();
      // quang: mettre a jour les sources avec la valeur initiale de champ inc
      equation().sources().mettre_a_jour(0);
    }
  else
    {
      Loi_Fermeture_base::preparer_calcul();
    }
}

void Loi_Fermeture_Diffusion_Nafion::mettre_a_jour(double temps)
{
  Cerr << "Loi_Fermeture_Diffusion_Nafion::mettre_a_jour " << equation().probleme().le_nom() << finl;
  // mettre a jour les champs couples -> interpoler vers P0
  const Zone_VF& la_zone = ref_cast(Zone_VF, equation().zone_dis().valeur());
  const DoubleTab& xp=la_zone.xp(); // Recuperation des centre de gravite des elements pour P0
  if(ch_T_.non_nul())
    {
      ch_T_.valeur().mettre_a_jour(temps);
      ch_T_.valeur().valeur_aux(xp, T_);
      T_.echange_espace_virtuel();
    }
  else if(temperature_.non_nul())
    {
      temperature_.valeur().mettre_a_jour(temps);
      temperature_.valeur().valeur_aux(xp, T_);
      T_.echange_espace_virtuel();
    }
  /*  if(ch_I_.non_nul())
      {
        ch_I_.valeur().mettre_a_jour(temps);
        ch_I_.valeur().valeur_aux(xp, I_);
      }
  // interpolation champ_inc P1NC -> P0 if necessaire
  if(ch_C_.valeur().que_suis_je().find("P1NC") !=-1)
  assert(sub_type(Champ_P1NC, ch_C_.valeur().valeur()));
  if(sub_type(Champ_P1NC, ch_C_.valeur().valeur()))
    {
      Champ_P1NC& C_faces = ref_cast(Champ_P1NC, ch_C_.valeur());
      C_faces.valeur_aux(xp, C_);
      C_.echange_espace_virtuel();
    }
  */

  ch_C_.valeur().valeur().valeur_aux(xp, C_);
  C_.echange_espace_virtuel();
  // mettre a jour le champ diffusivite "effective"
  DoubleTab& diffu = diffu_.valeurs();
  DoubleTab& por = por_naf_.valeurs();
  DoubleTab& eps = eps_naf_.valeurs();
  DoubleTab& tor = tor_naf_.valeurs();
  int nb_elem = diffu.dimension(0);
  for (int elem = 0; elem < nb_elem; ++elem)
    {
      diffu(elem) = eval_D_i_eff(T_(elem), C_(elem), por(elem,0), eps(elem,0), tor(elem,0));
    }
  diffu_.mettre_a_jour(temps);

  // mettre a jour le champ diffusivite "normal" D_i_naf
  DoubleTab& D_i_naf = D_i_naf_.valeurs();
  for (int elem = 0; elem < nb_elem; ++elem)
    {
      D_i_naf(elem) = eval_D_i_naf(T_(elem), C_(elem));
    }
  D_i_naf_.mettre_a_jour(temps);

  /*
    // mettre a jour le champ N_i_naf
    DoubleTab& N_i_naf = N_i_naf_.valeurs();
    // calcul du gradient
    DoubleTab grad(0, 1, Objet_U::dimension);
    equation().zone_dis().zone().creer_tableau_elements(grad);
    const DoubleTab& nu=diffu_.valeurs();
    Champ_P1NC::calcul_gradient(ch_C_.valeur().valeurs(),grad,ref_cast(Zone_Cl_VEF,equation().zone_Cl_dis().valeur()));
    grad.echange_espace_virtuel();
    if(nom_espece_ != "vap" && nom_espece_ != "H2O")
      {
        for (int elem=0; elem<grad.dimension(0); elem++)
          {
            for (int j=0; j<dimension; j++)
              {
                N_i_naf(elem,j)= - nu(elem)*grad(elem,0,j);			// flux N_i = -D_i_eff.grad(C)
              }
          }
      }
    else
      {

        for (int elem=0; elem<grad.dimension(0); elem++)
          {
            double ld = C_(elem)/CSO3_;
            double F = 96500;	// nombre de Faraday
            double nd = 1.0 + 0.0028*ld + 0.0026*ld*ld;
            for (int j=0; j<dimension; j++)
              {
                N_i_naf(elem,j) = nd/F*I_(elem, j) - nu(elem)*grad(elem,0,j);	// flux N_H2O = nd/F*I_i - D_w.grad(C)
              }
          }
      }
  */

  //Cerr << "Loi_Fermeture_Diffusion_Nafion::mettre_a_jour" << finl;
  Cerr << "Coefficient diffusion effective D_i_eff min max " << mp_min_vect(diffu) << " " << mp_max_vect(diffu) << finl;
  Cerr << "Coefficient diffusion defaut D_i_naf min max " << mp_min_vect(D_i_naf) << " " << mp_max_vect(D_i_naf) << finl;
}

void Loi_Fermeture_Diffusion_Nafion::set_param(Param& param)
{
  param.ajouter("nom_espece",&nom_espece_,Param::REQUIRED);	// XD_ADD_P chaine in list of 'O2' 'H2' 'H2O' 'vap' 'N2'
  param.ajouter("por_naf", &por_naf_, Param::REQUIRED);		// XD_ADD_P Champ_base porosity
  param.ajouter("eps_naf", &eps_naf_, Param::REQUIRED);	// XD_ADD_P Champ_base ionomer proportion %
  param.ajouter("tor_naf", &tor_naf_, Param::REQUIRED); // XD_ADD_P Champ_base tortuosity
  param.ajouter("nom_pb_T",&nom_pb_T_,Param::OPTIONAL); // XD_ADD_P chaine problem of temperature
  param.ajouter("nom_champ_T",&nom_champ_T_,Param::OPTIONAL);	// XD_ADD_P chaine default 'temperature'
  param.ajouter("temperature", &temperature_, Param::OPTIONAL); // XD_ADD_P Champ_Don if temperature is given
  //param.ajouter("CSO3", &CSO3_, Param::OPTIONAL);	// XD_ADD_P double default value of concentration mol/m3
  param.ajouter("nom_champ_Ceq", &nom_champ_Ceq_, Param::OPTIONAL);	// XD_ADD_P chaine default 'Ceq' in case of initialisation Co=Ceq
}

double Loi_Fermeture_Diffusion_Nafion::eval_D_i_naf(double T, double C)
{
  if (nom_espece_ == "H2")
    {
      return 4.1e-7*exp((-2602.)/T);
    }
  else if (nom_espece_ == "O2")
    {
      return 4.24e-6*exp(-2246. / T);
    }
  else if (nom_espece_ == "N2")
    {
      return 3.1e-7 * exp(-2768. /T);
    }
  else if (nom_espece_ == "H2O" || nom_espece_ == "vap")
    {
      double ld = C / C_SO3;
      return (6.707e-8*ld + 6.387e-7)*exp(-2416. / T);
    }
  else
    {
      Cerr <<" unknown species in the list "<<finl;
      Process::exit();
    }
  // pour compilos.
  return 0.;
}

double Loi_Fermeture_Diffusion_Nafion::eval_D_i_eff(double T, double C, double por, double eps, double tor)
{
  if (nom_espece_ == "H2O" || nom_espece_ == "vap")
    {
      return eval_D_i_naf(T,C);
    }

  double coef = (1. - por)*eps / (tor*tor);
  return eval_D_i_naf(T,C)*coef;
}
