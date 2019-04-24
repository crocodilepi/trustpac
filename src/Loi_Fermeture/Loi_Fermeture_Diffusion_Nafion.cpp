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
//  T_0_ = 353.15;				// dans le cas T constant
//  CSO3_ = 2036.;				// au cas ou "H2O" ou "vap"
//  por_naf_ = 0.47;			// porosite de Nafion
//  eps_naf_ = 0.2;			// ionomer proportionnel de Nafion
//  tor_naf_ = 1.;

  T_.resize(0, 1);			// scalaire
  C_.resize(0, 1);			// scalaire
  I_.resize(0, dimension);	// vectoriel
  equation().zone_dis().zone().creer_tableau_elements(T_);
  T_ = T_0_;
  equation().zone_dis().zone().creer_tableau_elements(C_);
  C_ = CSO3_;
  equation().zone_dis().zone().creer_tableau_elements(I_);		// TO-DO: initialize this array with a default value?
  I_ = 0.;
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

  dis.discretiser_champ("champ_elem",equation().zone_dis().valeur(),"N_i_naf","unit", dimension ,0. , N_i_naf_);
  N_i_naf_->fixer_nature_du_champ(vectoriel);
  champs_compris_.ajoute_champ(N_i_naf_);
}

void Loi_Fermeture_Diffusion_Nafion::completer()
{
  Loi_Fermeture_base::completer();
  // get the reference to the coupling fields
  if(nom_champ_T_ != "??" && nom_pb_T_ != "??")
    {
      Probleme_base& pb_T = ref_cast(Probleme_base,interprete().objet(nom_pb_T_));
      ch_T_ = pb_T.get_champ(nom_champ_T_);
    }
  if(nom_pb_phi_ != "??" && nom_champ_I_ != "??")
    {
      assert(nom_espece_ == "vap" || nom_espece_ == "H2O" );
      Probleme_base& pb_phi = ref_cast(Probleme_base,interprete().objet(nom_pb_phi_));
      ch_I_ = pb_phi.get_champ(nom_champ_I_);
    }
  ch_C_ = equation().inconnue();
}

void Loi_Fermeture_Diffusion_Nafion::preparer_calcul()
{
  Loi_Fermeture_base::preparer_calcul();
}

void Loi_Fermeture_Diffusion_Nafion::mettre_a_jour(double temps)
{
  // mettre a jour les champs couples -> interpoler vers P0
  const Zone_VF& la_zone = ref_cast(Zone_VF, equation().zone_dis());
  const DoubleTab& xp=la_zone.xp(); // Recuperation des centre de gravite des elements pour P0
  if(ch_T_.non_nul())
    {
      ch_T_.valeur().mettre_a_jour(temps);
      ch_T_.valeur().valeur_aux(xp, T_);
    }
  if(ch_I_.non_nul())
    {
      ch_I_.valeur().mettre_a_jour(temps);
      ch_I_.valeur().valeur_aux(xp, I_);
    }
  // interpolation champ_inc P1NC -> P0 if necessaire
  if(ch_C_.que_suis_je().find("P1NC") !=-1)
    {
      Champ_P1NC& C_faces = ref_cast(Champ_P1NC, ch_C_.valeur());
      C_faces.valeur_aux(xp, C_);
    }

  // mettre a jour le champ diffusivite "effective"
  DoubleTab& diffu = diffu_.valeurs();
  int nb_elem = diffu.dimension(0);
  for (int elem = 0; elem < nb_elem; ++elem)
    {
      diffu(elem) = eval_diffu_(T_(elem), C_(elem));
    }

  // mettre a jour le champ D_i_naf
  DoubleTab& D_i_naf = D_i_naf_.valeurs();
  for (int elem = 0; elem < nb_elem; ++elem)
    {
      D_i_naf(elem) = eval_D_i_naf(T_(elem), C_(elem));
    }

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
              N_i_naf(elem,j)=nu(elem)*grad(elem,0,j);
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
              N_i_naf(elem,j)=nd/F*I_(elem, j)-nu(elem)*grad(elem,0,j);
            }
        }
    }


  Cerr << "Loi_Fermeture_Diffusion_Nafion::mettre_a_jour" << finl;
  Cerr << "diffu min max " << mp_min_vect(diffu)<< " " <<mp_max_vect(diffu) << finl;
}

void Loi_Fermeture_Diffusion_Nafion::set_param(Param& param)
{
  param.ajouter("nom_espece",&nom_espece_,Param::REQUIRED);
  param.ajouter("T_0", &T_0_, Param::REQUIRED);
  param.ajouter("por_naf", &por_naf_, Param::REQUIRED);
  param.ajouter("eps_naf", &eps_naf_, Param::REQUIRED);
  param.ajouter("tor_naf", &tor_naf_, Param::REQUIRED);
  param.ajouter("nom_pb_T",&nom_pb_T_,Param::OPTIONAL);
  param.ajouter("nom_champ_T",&nom_champ_T_,Param::OPTIONAL);
  param.ajouter("CSO3", &CSO3_, Param::OPTIONAL);
  param.ajouter("nom_pb_phi",&nom_pb_phi_,Param::OPTIONAL);
  param.ajouter("nom_champ_I",&nom_champ_I_,Param::OPTIONAL);
}

double Loi_Fermeture_Diffusion_Nafion::eval_D_i_naf(double T, double C)
{
  if (nom_espece_ == "H2")
    {
      return 4.1e-7*exp((-2602.)/T);
    }
  else if (nom_espece_ == "02")
    {
      return 4.24e-6*exp(-2246. / T);
    }
  else if (nom_espece_ == "N2")
    {
      return 3.1e-7 * exp(-2768. /T);
    }
  else if (nom_espece_ == "H2O" || nom_espece_ == "vap")
    {
      double ld = C / CSO3_;
      return (6.707e-8*ld + 6.387e-7)*exp(-2416. / T);
    }
  // pour compilos.
  return 0.;
}

double Loi_Fermeture_Diffusion_Nafion::eval_diffu_(double T, double C)
{
  if (nom_espece_ == "H2O" || nom_espece_ == "vap")
    {
      return -eval_D_i_naf(T,C);
    }

  double coef = (1. - por_naf_)*eps_naf_ / (tor_naf_*tor_naf_);
  return eval_D_i_naf(T,C)*coef;
}
