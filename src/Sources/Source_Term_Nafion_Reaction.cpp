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
#include <Conduction.h>
#include <Domaine.h>
#include <Zone_VEF.h>
#include <Zone_Cl_VEF.h>
#include <Interprete.h>

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

  return is;
}

void Source_Term_Nafion_Reaction::associer_pb(const Probleme_base& pb)
{
  Cerr << " Source_Term_Nafion_Reaction::associer_pb " << finl ;
  assert(pb.que_suis_je() == "Pb_Conduction");
}

void Source_Term_Nafion_Reaction::associer_zones(const Zone_dis& zone_dis, const Zone_Cl_dis& zcl_dis)
{
  // cast VEF
  la_zone_VEF = ref_cast(Zone_VEF,zone_dis.valeur());
  la_zcl_VEF = ref_cast(Zone_Cl_VEF,zcl_dis.valeur());
  remplir_volumes();
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
  Probleme_base& pb_phi = ref_cast(Probleme_base,interprete().objet(nom_pb_phi_));
  ch_ir_ = pb_phi.get_champ(nom_champ_ir_);
  ch_ip_ = pb_phi.get_champ(nom_champ_ip_);
  la_zone_VEF.valeur().creer_tableau_faces(ir_);
  la_zone_VEF.valeur().creer_tableau_faces(ip_);
}

void Source_Term_Nafion_Reaction::mettre_a_jour(double temps)
{
  ch_ir_.valeur().mettre_a_jour(temps);
  ch_ip_.valeur().mettre_a_jour(temps);
  const DoubleTab& xv=la_zone_VEF.valeur().xv(); // centre de gravite des faces pour P1NC
  ch_ir_.valeur().valeur_aux( xv, ir_ );
  ch_ip_.valeur().valeur_aux( xv, ip_ );
}

DoubleTab& Source_Term_Nafion_Reaction::calculer(DoubleTab& resu) const
{
  resu = 0;
  ajouter(resu);
  return resu;
}

void Source_Term_Nafion_Reaction::remplir_volumes()
{
  volumes_.ref(la_zone_VEF.valeur().volumes_entrelaces());
}

DoubleTab& Source_Term_Nafion_Reaction::ajouter(DoubleTab& resu) const
{
  double inv_rhoCp = 1./((1-por_naf_)*eps_naf_);

  IntTab faces_ssz;	// faces belong to the sous_zone -> flag = 1, if not, flag = 0
  la_zone_VEF.valeur().creer_tableau_faces(faces_ssz);
  faces_ssz = 0;		// init with no flag (all faces are unchecked)
  // to-do
  if(CL_a_.non_nul())
    {
      for (int poly = 0; poly < CL_a_.valeur().nb_elem_tot(); poly++)
        {
          int elem = CL_a_.valeur()(poly);
          for (int f = 0; f < la_zone_VEF.valeur().zone().nb_faces_elem(0); ++f)
            {
              int face = la_zone_VEF.valeur().elem_faces(elem, f);
              if(!faces_ssz(face))
                {
                  resu(face) += eval_f(ir_(face), ip_(face)) * inv_rhoCp * volumes_(face);
                  faces_ssz(face) = 1;		// marquer comme deja traite
                }
            }
        }
      faces_ssz = 0;		// init with no flag (all faces are unchecked)
    }
  if(CL_c_.non_nul())
    {
      for (int poly = 0; poly < CL_c_.valeur().nb_elem_tot(); poly++)
        {
          int elem = CL_c_.valeur()(poly);
          for (int f = 0; f < la_zone_VEF.valeur().zone().nb_faces_elem(0); ++f)
            {
              int face = la_zone_VEF.valeur().elem_faces(elem, f);
              if(!faces_ssz(face))
                {
                  resu(face) += eval_f(ir_(face), ip_(face))* inv_rhoCp * volumes_(face);
                  faces_ssz(face) = 1;		// marquer comme deja traite
                }
            }
        }
      faces_ssz = 0;		// init with no flag (all faces are unchecked)
    }

  // Source = 0 with N2
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
