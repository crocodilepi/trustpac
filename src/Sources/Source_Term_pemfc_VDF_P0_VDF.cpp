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
// File      : Source_Term_pemfc_VDF_P0_VDF.cpp
// Directory : $PEMFC_ROOT/src/Sources
//
/////////////////////////////////////////////////////////////////////////////
#include <Source_Term_pemfc_VDF_P0_VDF.h>
#include <Zone_VDF.h>
#include <Zone_Cl_VDF.h>
#include <Equation_base.h>
#include <Champ_P1NC.h>

Implemente_instanciable( Source_Term_pemfc_VDF_P0_VDF, "Source_Term_pemfc_VDF_P0_VDF", Source_Term_pemfc_base );

Sortie& Source_Term_pemfc_VDF_P0_VDF::printOn(Sortie& os) const
{
  Source_Term_pemfc_base::printOn(os);
  return os;
}

Entree& Source_Term_pemfc_VDF_P0_VDF::readOn(Entree& is)
{
  Source_Term_pemfc_base::readOn(is);

  return is;
}

void Source_Term_pemfc_VDF_P0_VDF::associer_zones(const Zone_dis& zone_dis, const Zone_Cl_dis& zcl_dis)
{
  la_zone_VDF = ref_cast(Zone_VDF,zone_dis.valeur());
  la_zcl_VDF = ref_cast(Zone_Cl_VDF,zcl_dis.valeur());
  remplir_volumes();
}

void Source_Term_pemfc_VDF_P0_VDF::associer_pb(const Probleme_base& pb)
{
  Cerr << " Source_Term_pemfc_VDF_P0_VDF::associer_pb " << finl ;
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
//      Cerr << "Erreur TRUST dans Source_Term_pemfc_VDF_P0_VDF::associer_pb()" << finl;
//      Cerr << "On ne trouve pas d'equation de conduction dans le probleme" << finl;
//      exit();
//    }

}

DoubleTab& Source_Term_pemfc_VDF_P0_VDF::ajouter(DoubleTab& resu) const
{
  assert(resu.dimension(0)==volumes_.size());
  assert(resu.dimension(0)==T_.size());
  assert(resu.dimension(0)==C_.size());
  assert(resu.dimension(0)==diffu_.size());
  assert(resu.dimension(0)==ci_.size());

  double inv_rhoCp = 1./((1-por_naf_)*eps_naf_);

  for (int poly = 0; poly < CL_a_.valeur().nb_elem_tot(); poly++)
    {
      int elem = CL_a_.valeur()(poly);
      resu(elem) = eval_f(diffu_(elem), C_(elem), ci_(elem), T_(elem)) * inv_rhoCp;
    }
  for (int poly = 0; poly < CL_c_.valeur().nb_elem_tot(); poly++)
    {
      int elem = CL_c_.valeur()(poly);
      resu(elem) = eval_f(diffu_(elem), C_(elem), ci_(elem), T_(elem)) * inv_rhoCp;
    }
  return resu;
}

void Source_Term_pemfc_VDF_P0_VDF::contribuer_a_avec(const DoubleTab& inco, Matrice_Morse& mat) const
{
  assert(inco.dimension(0)==volumes_.size());
  assert(inco.dimension(0)==diffu_.size());

  double inv_rhoCp = 1./((1-por_naf_)*eps_naf_);

  for (int poly = 0; poly < CL_a_.valeur().nb_elem_tot(); poly++)
    {
      int elem = CL_a_.valeur()(poly);
      mat.coef(elem,elem) += volumes_(elem) * eval_derivee_f(diffu_(elem))*inv_rhoCp;
    }
  for (int poly = 0; poly < CL_c_.valeur().nb_elem_tot(); poly++)
    {
      int elem = CL_c_.valeur()(poly);
      mat.coef(elem,elem) += volumes_(elem) * eval_derivee_f(diffu_(elem))*inv_rhoCp;
    }
}

void Source_Term_pemfc_VDF_P0_VDF::completer()
{
  Source_Term_pemfc_base::completer();
  T_.resize(0, 1);			// scalaire
  C_.resize(0, 1);			// scalaire
  ci_.resize(0, 1);			// scalaire
  diffu_.resize(0,1);		// scalaire
  ir_.resize(0,1);			// scalaire
  ip_.resize(0,1);			// scalaire

  dom_.valeur().creer_tableau_elements(C_);
  dom_.valeur().creer_tableau_elements(diffu_);
  dom_.valeur().creer_tableau_elements(ci_);
  dom_.valeur().creer_tableau_elements(T_);
  dom_.valeur().creer_tableau_elements(ir_);
  dom_.valeur().creer_tableau_elements(ip_);
}

void Source_Term_pemfc_VDF_P0_VDF::mettre_a_jour(double temps)
{

  Source_Term_pemfc_base::mettre_a_jour(temps);

  const DoubleTab& xp=la_zone_VDF.valeur().xp(); // Recuperation des centre de gravite des elements pour P0

  if(ch_C_.valeur().valeur().que_suis_je().find("P0") !=-1)
    {
      C_.ref(ch_C_.valeur().valeurs());
    }
  else
    {
      Champ_P1NC& ch_C = ref_cast(Champ_P1NC, ch_C_.valeur().valeur());
      ch_C.valeur_aux(xp, C_);
    }

  if(ch_D_i_naf_.valeur().que_suis_je().find("P0") !=-1)
    {
      diffu_.ref(ch_D_i_naf_.valeur().valeurs());
    }
  else
    {
      Champ_P1NC& ch_D = ref_cast(Champ_P1NC, ch_D_i_naf_.valeur());
      ch_D.valeur_aux(xp, diffu_);
    }

  if(ch_T_.non_nul())
    {
      ch_T_.valeur().valeur_aux(xp, T_);
    }
  else
    {
      T_ = T_0_;
    }

  if(ch_ci_cathode_.non_nul() && !ch_ci_anode_.non_nul())
    {
      // case O2
      ch_ci_cathode_.valeur().valeur_aux(xp, ci_);
    }
  else if(ch_ci_anode_.non_nul() && !ch_ci_cathode_.non_nul())
    {
      // case H2
      ch_ci_anode_.valeur().valeur_aux(xp, ci_);
    }
  else
    {
      // case N2, H20
      DoubleTab val_ci_cathode, val_ci_anode;
      ch_ci_cathode_.valeur().valeur_aux(xp, val_ci_cathode);
      ch_ci_anode_.valeur().valeur_aux(xp, val_ci_anode);

      for (int poly = 0; poly < CL_a_.valeur().nb_elem_tot(); poly++)
        {
          int elem = CL_a_.valeur()(poly);
          ci_(elem) = val_ci_anode(elem);
        }

      for (int poly = 0; poly < CL_c_.valeur().nb_elem_tot(); poly++)
        {
          int elem = CL_c_.valeur()(poly);
          ci_(elem)= val_ci_cathode(elem);
        }
    }

  if(ch_ir_.non_nul())
    {
      ch_ir_.valeur().valeur_aux(xp, ir_);
    }
//	  else
//	    {
//	      ir_ = 0.;
//	    }

  if(ch_ip_.non_nul())
    {
      ch_ip_.valeur().valeur_aux(xp, ip_);
    }
//	  else
//	    {
//	      ip_ = 0.;
//	    }

  if(ch_op_.non_nul())
    {
      ch_op_.valeur().valeur_aux(xp, op_);
    }
  //  else
  //    {
  //      op_ = 0.;
  //    }
}

void Source_Term_pemfc_VDF_P0_VDF::remplir_volumes()
{
  volumes_.ref(ref_cast(Zone_VF,equation().zone_dis().valeur()).volumes_entrelaces());
}
