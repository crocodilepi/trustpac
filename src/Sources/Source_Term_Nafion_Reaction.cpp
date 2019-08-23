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
// File      : Source_Term_Nafion_Reaction.cpp
// Directory : $PEMFC_ROOT/src/Sources
//
/////////////////////////////////////////////////////////////////////////////

#include <Source_Term_Nafion_Reaction.h>
#include <Probleme_base.h>
#include <Equation_base.h>
#include <Interprete.h>
#include <Domaine.h>
#include <Zone_VF.h>
#include <Matrice_Morse.h>

Implemente_instanciable( Source_Term_Nafion_Reaction, "Source_Term_Nafion_Reaction_VEF_P1NC", Source_base ) ;

Sortie& Source_Term_Nafion_Reaction::printOn( Sortie& os ) const
{
  Source_base::printOn( os );
  return os;
}

Entree& Source_Term_Nafion_Reaction::readOn( Entree& is )
{
  Source_base::readOn( is );
  Cerr << " Source_Term_Nafion_Reaction::readOn " << finl  ;
  Param param(que_suis_je());
  set_param(param);
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

  dom_ = equation().probleme().domaine();
  assert(!(nom_pb_phi_!="??" && reacCurrent_.non_nul() && permCurrent_.non_nul()));	// either coupling either independence

  if(nom_espece_ == "H2")
    {
      // check
      assert(nom_ssz_CLa_ != "??");
    }
  else if (nom_espece_ == "O2")
    {
      assert(nom_ssz_CLc_ != "??");
    }
  else
    {
      assert(nom_ssz_CLa_ != "??");
      assert(nom_ssz_CLc_ != "??");
    }

  if(nom_ssz_CLa_ != "??")
    CL_a_ = dom_.valeur().ss_zone(nom_ssz_CLa_);
  if(nom_ssz_CLc_ != "??")
    CL_c_ = dom_.valeur().ss_zone(nom_ssz_CLc_);

  dom_.valeur().creer_tableau_elements(ir_);
  dom_.valeur().creer_tableau_elements(ip_);
  dom_.valeur().creer_tableau_elements(DirDcH2_);
  dom_.valeur().creer_tableau_elements(DirDcO2_);
  dom_.valeur().creer_tableau_elements(DirDcH2O_);
  if(reacCurrent_.non_nul())
    {
      assert(permCurrent_.non_nul());
      DirDcH2_ = 0;
      DirDcO2_ = 0;
      DirDcH2O_ = 0;
    }
  return is;
}

void Source_Term_Nafion_Reaction::associer_pb(const Probleme_base& pb)
{
  Cerr << " Source_Term_Nafion_Reaction::associer_pb " << finl ;
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
//      Cerr << "Erreur TRUST dans Source_Term_Nafion_Reaction::associer_pb()" << finl;
//      Cerr << "On ne trouve pas d'equation de conduction dans le probleme" << finl;
//      exit();
//    }
}

void Source_Term_Nafion_Reaction::associer_zones(const Zone_dis& zone_dis, const Zone_Cl_dis& zcl_dis)
{
  la_zone_ = ref_cast(Zone_VF,zone_dis.valeur());
}

void Source_Term_Nafion_Reaction::set_param(Param& param)
{
  param.ajouter("nom_espece",&nom_espece_,Param::REQUIRED);    // XD_ADD_P chaine in list 'O2' 'H2' 'H2O' 'vap' 'N2'
  param.ajouter("nom_CLa",&nom_ssz_CLa_,Param::OPTIONAL);// XD_ADD_P chaine sub-area where the source exists, required for H2, N2, H2O
  param.ajouter("nom_CLc",&nom_ssz_CLc_,Param::OPTIONAL);// XD_ADD_P chaine sub-area where the source exists, required for O2, N2, H2O
  param.ajouter("nom_pb_phi", &nom_pb_phi_, Param::OPTIONAL);
  param.ajouter("nom_champ_ir", &nom_champ_ir_, Param::OPTIONAL);
  param.ajouter("nom_champ_ip", &nom_champ_ip_, Param::OPTIONAL);
  param.ajouter("reacCurrent", &reacCurrent_, Param::OPTIONAL);
  param.ajouter("permCurrent", &permCurrent_, Param::OPTIONAL);
  param.ajouter("por_naf", &por_naf_, Param::REQUIRED);
  param.ajouter("eps_naf", &eps_naf_, Param::REQUIRED);
}

void Source_Term_Nafion_Reaction::completer()
{
  Cerr << "Source_Term_Nafion_Reaction::completer " << equation().probleme().le_nom() << finl;
  Source_base::completer();
  // get the reference to the coupling fields
  Source_base::completer();
  if(nom_pb_phi_ != "??")
    {
      Probleme_base& pb_phi = ref_cast(Probleme_base,interprete().objet(nom_pb_phi_));
      ch_ir_ = pb_phi.get_champ(nom_champ_ir_);
      assert(ch_ir_.valeur().que_suis_je().find("P0") !=-1);
      ch_ip_ = pb_phi.get_champ(nom_champ_ip_);
      assert(ch_ip_.valeur().que_suis_je().find("P0") !=-1);

      ch_DirDcO2_ = pb_phi.get_champ("DirDcO2");
      ch_DirDcH2_ = pb_phi.get_champ("DirDcH2");
      ch_DirDcH2O_ = pb_phi.get_champ("DirDcH2O");
    }
  F_ = 96500;
}

void Source_Term_Nafion_Reaction::mettre_a_jour(double temps)
{
  const DoubleTab& xp=la_zone_.valeur().xp(); // Recuperation des centre de gravite des elements pour P0
  if(ch_ir_.non_nul())
    {
      ch_ir_.valeur().mettre_a_jour(temps);
      ch_ir_.valeur().valeur_aux( xp, ir_ );			// ir
      ir_.echange_espace_virtuel();
    }
  if(ch_ip_.non_nul())
    {
      ch_ip_.valeur().mettre_a_jour(temps);
      ch_ip_.valeur().valeur_aux( xp, ip_ );		    // ip
      ip_.echange_espace_virtuel();
    }
  if(ch_DirDcO2_.non_nul())
    {
      ch_DirDcO2_.valeur( ).mettre_a_jour( temps );
      ch_DirDcO2_.valeur().valeur_aux( xp, DirDcO2_ );
      DirDcO2_.echange_espace_virtuel();
    }
  if (ch_DirDcH2_.non_nul())
    {
      ch_DirDcH2_.valeur( ).mettre_a_jour( temps );
      ch_DirDcH2_.valeur().valeur_aux( xp, DirDcH2_ );
      DirDcH2_.echange_espace_virtuel();
    }
  if (ch_DirDcH2O_.non_nul())
    {
      ch_DirDcH2O_.valeur( ).mettre_a_jour( temps );
      ch_DirDcH2O_.valeur().valeur_aux( xp, DirDcH2O_ );
      DirDcH2O_.echange_espace_virtuel();
    }
  if(reacCurrent_.non_nul() && permCurrent_.non_nul())
    {
      reacCurrent_.valeur().mettre_a_jour(temps);
      reacCurrent_.valeur().valeur_aux( xp, ir_ );			// ir
      ir_.echange_espace_virtuel();
      permCurrent_.valeur().mettre_a_jour(temps);
      permCurrent_.valeur().valeur_aux( xp, ip_ );		    // ip
      ip_.echange_espace_virtuel();
    }
  Cerr << "ch_ir min max " << mp_min_vect(ir_) << " " << mp_max_vect(ir_) << finl;
  Cerr << "ch_ip min max " << mp_min_vect(ip_) << " " << mp_max_vect(ip_) << finl;
}

DoubleTab& Source_Term_Nafion_Reaction::calculer(DoubleTab& resu) const
{
  resu = 0;
  ajouter(resu);
  return resu;
}

DoubleTab& Source_Term_Nafion_Reaction::ajouter(DoubleTab& resu) const
{
  assert(resu.dimension(0)==la_zone_.valeur().nb_faces());

  DoubleVect vol = la_zone_.valeur().volumes();
  const DoubleTab& por = por_naf_.valeurs();
  const DoubleTab& eps = eps_naf_.valeurs();

  if(CL_a_.non_nul())
    {
      for (int poly = 0; poly < CL_a_.valeur().nb_elem_tot(); ++poly)
        {
          int elem = CL_a_.valeur()(poly);
          double coeff = (1-por(elem,0))*eps(elem,0);
          double S = eval_f_anode(ir_(elem), ip_(elem));

          int nb_face_elem = la_zone_.valeur().zone().nb_faces_elem(0);
          for (int f = 0; f < nb_face_elem; ++f)
            {
              int face = la_zone_.valeur().elem_faces(elem, f);
              resu(face) += S / coeff * vol(elem)/ nb_face_elem;
            }
        }
    }

  if(CL_c_.non_nul())
    {
      for (int poly = 0; poly < CL_c_.valeur().nb_elem_tot(); poly++)
        {
          int elem = CL_c_.valeur()(poly);
          double coeff = (1-por(elem,0))*eps(elem,0);
          double S = eval_f_cathode(ir_(elem), ip_(elem));

          int nb_face_elem = la_zone_.valeur().zone().nb_faces_elem(0);
          for (int f = 0; f < nb_face_elem; ++f)
            {
              int face = la_zone_.valeur().elem_faces(elem, f);
              resu(face) += S / coeff * vol(elem)/ nb_face_elem;
            }
        }
    }

  return resu;
}

double Source_Term_Nafion_Reaction::eval_f_anode(double ir, double ip) const
{
  if (nom_espece_ == "H2")
    {
      return - ir / (2. * F_);
    }
  return 0.;
}

double Source_Term_Nafion_Reaction::eval_f_cathode(double ir, double ip) const
{
  if (nom_espece_ == "O2")
    {
      return (ir + ip) / (4. * F_);
    }
  else if (nom_espece_ == "H2O" || nom_espece_ == "vap")
    {
      return -(ir + ip)/(2. * F_);
    }

  return 0.;
}

void Source_Term_Nafion_Reaction::contribuer_a_avec(const DoubleTab& inco, Matrice_Morse& mat) const
{

  const DoubleTab& por = por_naf_.valeurs();
  const DoubleTab& eps = eps_naf_.valeurs();

  DoubleVect vol = la_zone_.valeur().volumes();

  if(CL_a_.non_nul())
    {
      for (int poly = 0; poly < CL_a_.valeur().nb_elem_tot(); ++poly)
        {
          int elem = CL_a_.valeur()(poly);
          double coeff = (1-por(elem,0))*eps(elem,0);
          double S = eval_df_anode( elem );

          int nb_face_elem = la_zone_.valeur().zone().nb_faces_elem(0);
          for (int f = 0; f < nb_face_elem; ++f)
            {
              int face = la_zone_.valeur().elem_faces(elem, f);
              mat.coef(face,face) += S / coeff * vol(elem)/ nb_face_elem;
            }
        }
    }

  if(CL_c_.non_nul())
    {
      for (int poly = 0; poly < CL_c_.valeur().nb_elem_tot(); poly++)
        {
          int elem = CL_c_.valeur()(poly);
          double coeff = (1-por(elem,0))*eps(elem,0);
          double S = eval_df_cathode( elem );

          int nb_face_elem = la_zone_.valeur().zone().nb_faces_elem(0);
          for (int f = 0; f < nb_face_elem; ++f)
            {
              int face = la_zone_.valeur().elem_faces(elem, f);
              mat.coef(face,face) += S / coeff * vol(elem)/ nb_face_elem;
            }
        }
    }

}

double Source_Term_Nafion_Reaction::eval_df_anode( const int& elem ) const
{
  if (nom_espece_ == "H2")
    {
      return  DirDcH2_( elem )/( 2.*F_ ) ;
    }
  return 0.;
}

double Source_Term_Nafion_Reaction::eval_df_cathode( const int& elem ) const
{
  if (nom_espece_ == "O2")
    {
      return  -DirDcO2_( elem )/(4. * F_);
    }
  else if (nom_espece_ == "H2O" || nom_espece_ == "vap")
    {
      return DirDcH2O_( elem )/(2. * F_);
    }

  return 0.;
}
