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
  return is;
}

void Source_Term_Nafion_Reaction::associer_pb(const Probleme_base& pb)
{
  Cerr << " Source_Term_Nafion_Reaction::associer_pb " << finl ;
  assert(pb.que_suis_je() == "Pb_Conduction");
}

void Source_Term_Nafion_Reaction::associer_zones(const Zone_dis& zone_dis, const Zone_Cl_dis& zcl_dis)
{
  la_zone_ = ref_cast(Zone_VF,zone_dis.valeur());
}

void Source_Term_Nafion_Reaction::set_param(Param& param)
{
  param.ajouter("nom_espece",&nom_espece_,Param::REQUIRED);
  //param.ajouter("nom_domaine",&nom_domaine_,Param::REQUIRED);		// pas necessaire ? dom_ = equation().problem().domaine() ?
  param.ajouter("nom_CLa",&nom_ssz_CLa_,Param::OPTIONAL);					// requis pour H2, N2, H20, sauf O2
  param.ajouter("nom_CLc",&nom_ssz_CLc_,Param::OPTIONAL);					// requis pour O2, N2, H20, sauf H2
  param.ajouter("nom_pb_phi", &nom_pb_phi_, Param::REQUIRED);
  param.ajouter("nom_champ_ir", &nom_champ_ir_, Param::REQUIRED);
  param.ajouter("nom_champ_ip", &nom_champ_ip_, Param::REQUIRED);
  param.ajouter("por_naf", &por_naf_, Param::REQUIRED);
  param.ajouter("eps_naf", &eps_naf_, Param::REQUIRED);
}

void Source_Term_Nafion_Reaction::completer()
{
  Source_base::completer();
  // get the reference to the coupling fields
  Source_base::completer();
  Probleme_base& pb_phi = ref_cast(Probleme_base,interprete().objet(nom_pb_phi_));
  ch_ir_ = pb_phi.get_champ(nom_champ_ir_);
  assert(ch_ir_.valeur().que_suis_je().find("P0") !=-1);
  ch_ip_ = pb_phi.get_champ(nom_champ_ip_);
  assert(ch_ip_.valeur().que_suis_je().find("P0") !=-1);
}

void Source_Term_Nafion_Reaction::mettre_a_jour(double temps)
{
  ch_ir_.valeur().mettre_a_jour(temps);
  ch_ip_.valeur().mettre_a_jour(temps);
  const DoubleTab& xp=la_zone_.valeur().xp(); // Recuperation des centre de gravite des elements pour P0
  ch_ir_.valeur().valeur_aux( xp, ir_ );			// ir
  ch_ip_.valeur().valeur_aux( xp, ip_ );		    // ip

  Cerr << "Source_Term_Nafion_Reaction::mettre_a_jour" << finl;
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
          int nb_face_elem = la_zone_.valeur().zone().nb_faces_elem(0);
          for (int f = 0; f < nb_face_elem; ++f)
            {
              int face = la_zone_.valeur().elem_faces(elem, f);
              resu(face) += eval_f(ir_(elem), ip_(elem)) / coeff * vol(elem)/ nb_face_elem;
            }
        }
    }

  if(CL_c_.non_nul())
    {
      for (int poly = 0; poly < CL_c_.valeur().nb_elem_tot(); poly++)
        {
          int elem = CL_c_.valeur()(poly);
          double coeff = (1-por(elem,0))*eps(elem,0);
          int nb_face_elem = la_zone_.valeur().zone().nb_faces_elem(0);
          for (int f = 0; f < nb_face_elem; ++f)
            {
              int face = la_zone_.valeur().elem_faces(elem, f);
              resu(face) += eval_f(ir_(elem), ip_(elem)) / coeff * vol(elem)/ nb_face_elem;
            }
        }
    }

  return resu;
}

double Source_Term_Nafion_Reaction::eval_f(double ir, double ip) const
{
  double F = 96500;
  if (nom_espece_ == "H2")
    {
      return - ir / (2. * F);
    }
  else if (nom_espece_ == "02")
    {
      return (ir + ip) / (4. * F);
    }
  else if (nom_espece_ == "N2")
    {
      return 0.;
    }
  else if (nom_espece_ == "H2O" || nom_espece_ == "vap")
    {
      return -(ir + ip)/(2. * F);
    }
  return 0.;
}
