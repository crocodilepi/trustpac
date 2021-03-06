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
// File      : Source_Term_pemfc_VEF_P1NC.cpp
// Directory : $PEMFC_ROOT/src/Sources
//
/////////////////////////////////////////////////////////////////////////////
#include <Source_Term_pemfc_VEF_P1NC.h>
#include <Zone_VEF.h>
#include <Zone_Cl_VEF.h>
#include <Probleme_base.h>
#include <Equation_base.h>
#include <Conduction.h>
#include <Champ_P1NC.h>
#include <Champ_Generique_base.h>
#include <Champ_P0_VEF.h>

Implemente_instanciable( Source_Term_pemfc_VEF_P1NC, "Source_Term_pemfc_VEF_P1NC", Source_Term_pemfc_base );

Sortie& Source_Term_pemfc_VEF_P1NC::printOn(Sortie& os) const
{
  Source_Term_pemfc_base::printOn(os);
  return os;
}

Entree& Source_Term_pemfc_VEF_P1NC::readOn(Entree& is)
{
  Source_Term_pemfc_base::readOn(is);
  return is;
}

void Source_Term_pemfc_VEF_P1NC::associer_zones(const Zone_dis& zone_dis, const Zone_Cl_dis& zcl_dis)
{
  la_zone_VEF = ref_cast(Zone_VEF,zone_dis.valeur());
  la_zcl_VEF = ref_cast(Zone_Cl_VEF,zcl_dis.valeur());
  remplir_volumes();
}

void Source_Term_pemfc_VEF_P1NC::associer_pb(const Probleme_base& pb)
{
  Cerr << " Source_Term_pemfc_VEF_P1NC::associer_pb " << finl ;
  assert(pb.que_suis_je() == "Pb_Conduction");
}

// A VERIFIER
DoubleTab& Source_Term_pemfc_VEF_P1NC::ajouter(DoubleTab& resu) const
{
  assert(resu.dimension(0)==volumes_.size());
  assert(resu.dimension(0)==T_.size());
  assert(resu.dimension(0)==C_.size());
  assert(resu.dimension(0)==diffu_.size());
  assert(resu.dimension(0)==ci_.size());

  double inv_rhoCp = 1./((1-por_naf_)*eps_naf_);

  IntTab faces_ssz;	// faces belong to the sous_zone -> flag = 1, if not, flag = 0
  la_zone_VEF.valeur().creer_tableau_faces(faces_ssz);
  faces_ssz = 0;	// init with no flag (all faces are unchecked)

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
                  resu(face) += eval_f(diffu_(face), C_(face), ci_(face), T_(face)) * volumes_(face) * inv_rhoCp;
                  // necessaire (source*porosite_surf(num_face));
                  faces_ssz(face) = 1;		// marquer comme deja traite
                }
            }
        }
    }

  faces_ssz = 0;		// init with no flag (all faces are unchecked)
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
                  resu(face) += eval_f(diffu_(face), C_(face), ci_(face), T_(face)) * volumes_(face) * inv_rhoCp;
                  // necessaire (source*porosite_surf(num_face));
                  faces_ssz(face) = 1;		// marquer comme deja traite
                }
            }
        }
    }
//
//  if(nom_espece_ == "H2O" || nom_espece_ == "vap")
//    {
//      int nb_faces = la_zone_VEF.valeur().nb_faces();
//      // ajouter un terme source de type: -nd/F*op avec op = div(-kappa.grad(phi))
//      for (int face = 0; face < nb_faces; ++face)
//        {
//          resu(face) += -f_nd(C_(face))/F * op_(face) * volumes_(face) * inv_rhoCp;			// A VERIFIER
//        }
//    }

  return resu;
}

// A VERIFIER, NECESAIRE OU PAS ?
void Source_Term_pemfc_VEF_P1NC::contribuer_a_avec(const DoubleTab& inco, Matrice_Morse& mat) const
{
  assert(inco.dimension(0)==volumes_.size());
  assert(inco.dimension(0)==T_.size());
  assert(inco.dimension(0)==C_.size());
  assert(inco.dimension(0)==diffu_.size());
  assert(inco.dimension(0)==ci_.size());

  IntTab faces_ssz;	// faces belong to the sous_zone -> flag = 1, if not, flag = 0
  la_zone_VEF.valeur().creer_tableau_faces(faces_ssz);
  faces_ssz = 0;		// init with no flag (all faces are unchecked)

  double inv_rhoCp = 1./((1-por_naf_)*eps_naf_);
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
                  mat.coef(face,face) += volumes_(face) * eval_derivee_f(diffu_(face))*inv_rhoCp;
                  // necessaire (source*porosite_surf(num_face));
                  faces_ssz(face) = 1;		// marquer comme deja traite
                }
            }
        }
    }
  faces_ssz = 0;		// init with no flag (all faces are unchecked)
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
                  mat.coef(face,face) += volumes_(face) * eval_derivee_f(diffu_(face))*inv_rhoCp;
                  // necessaire (source*porosite_surf(num_face));
                  faces_ssz(face) = 1;		// marquer comme deja traite
                }
            }
        }
    }
}

void Source_Term_pemfc_VEF_P1NC::completer()
{
  Source_Term_pemfc_base::completer();
  //T_.resize(0, 1);			// scalaire
  //C_.resize(0, 1);			// scalaire
  //ci_.resize(0, 1);			// scalaire
  //diffu_.resize(0,1);		    // scalaire
  //ir_.resize(0,1);			// scalaire
  //ip_.resize(0,1);			// scalaire
  //op_.resize(0,1);			// scalaire

  la_zone_VEF.valeur().creer_tableau_faces(C_);
  la_zone_VEF.valeur().creer_tableau_faces(diffu_);
  la_zone_VEF.valeur().creer_tableau_faces(ci_);
  la_zone_VEF.valeur().creer_tableau_faces(T_);
//  la_zone_VEF.valeur().creer_tableau_faces(ir_);
//  la_zone_VEF.valeur().creer_tableau_faces(ip_);
//  la_zone_VEF.valeur().creer_tableau_faces(op_);
}

// mettre a jour les valeurs suivants: diffu_, C_, T_, ci_, ir_, ip_, op_
void Source_Term_pemfc_VEF_P1NC::mettre_a_jour(double temps)
{

  Source_Term_pemfc_base::mettre_a_jour(temps);

  const DoubleTab& xv=la_zone_VEF.valeur().xv(); // Recuperation des centre de gravite des faces pour P1NC

  if(ch_C_.valeur().valeur().que_suis_je().find("P1NC") !=-1)
//  if(sub_type(Champ_P1NC, ch_C_.valeur().valeur()))
    {
      C_.ref(ch_C_.valeur().valeurs());
    }
  else
    {
      ch_C_.valeur().valeur().valeur_aux(xv, C_);
//      assert(sub_type(Champ_P0_VEF, ch_C_.valeur().valeur()));
//      Champ_P0_VEF& ch_C = ref_cast(Champ_P0_VEF, ch_C_.valeur().valeur());
//      ch_C.valeur_aux(xv, C_);
    }

//  if(sub_type(Champ_P1NC, ch_D_i_naf_.valeur()))
  if(ch_D_i_naf_.valeur().que_suis_je().find("P1NC") !=-1)
    {
      diffu_.ref(ch_D_i_naf_.valeur().valeurs());
    }
  else
    {
      ch_D_i_naf_.valeur().valeur_aux(xv, diffu_);
//      assert(sub_type(Champ_P0_VEF, ch_D_i_naf_.valeur()));
//      Champ_P0_VEF& ch_D = ref_cast(Champ_P0_VEF, ch_D_i_naf_.valeur());
//      ch_D.valeur_aux(xv, diffu_);
    }

  if(ch_T_.non_nul())
    {
      ch_T_.valeur().valeur_aux(xv, T_);
    }
  else
    {
      T_ = T_0_;
    }

  if(ch_ci_cathode_.non_nul() && !ch_ci_anode_.non_nul())
    {
      // case O2
      ch_ci_cathode_.valeur().valeur_aux_compo(xv, ci_, 0);										// ncomp = 0 pour O2
    }
  else if(ch_ci_anode_.non_nul() && !ch_ci_cathode_.non_nul())
    {
      // case H2
      ch_ci_anode_.valeur().valeur_aux_compo(xv, ci_, 0);										// ncomp = 0 pour O2
    }
  else
    {
      // case N2, H20
      DoubleTab val_ci_cathode, val_ci_anode;
      la_zone_VEF.valeur().creer_tableau_faces(val_ci_cathode);
      la_zone_VEF.valeur().creer_tableau_faces(val_ci_anode);
      if(nom_espece_ == "N2")
        {
          ch_ci_cathode_.valeur().valeur_aux_compo(xv, val_ci_cathode, 2);							// ncomp = 1 pour H20, 2 pour N2
          ch_ci_anode_.valeur().valeur_aux_compo(xv, val_ci_anode, 2);
        }
      else
        {
          ch_ci_cathode_.valeur().valeur_aux_compo(xv, val_ci_cathode, 1);							// ncomp = 1 pour H20, 2 pour N2
          ch_ci_anode_.valeur().valeur_aux_compo(xv, val_ci_anode, 1);
        }

      IntTab faces_ssz;	// faces belong to the sous_zone -> flag = 1, if not, flag = 0
      la_zone_VEF.valeur().creer_tableau_faces(faces_ssz);

      faces_ssz = 0;		// init with no flag (all faces are unchecked)
      for (int poly = 0; poly < CL_a_.valeur().nb_elem_tot(); poly++)
        {
          int elem = CL_a_.valeur()(poly);
          for (int f = 0; f < la_zone_VEF.valeur().zone().nb_faces_elem(0); ++f)
            {
              int face = la_zone_VEF.valeur().elem_faces(elem, f);
              if(!faces_ssz(face))
                {
                  ci_(face) = val_ci_anode(face);
                  // necessaire (source*porosite_surf(num_face));
                  faces_ssz(face) = 1;		// marquer comme deja traite
                }
            }
        }
      faces_ssz = 0;		// init with no flag (all faces are unchecked)
      for (int poly = 0; poly < CL_c_.valeur().nb_elem_tot(); poly++)
        {
          int elem = CL_c_.valeur()(poly);
          for (int f = 0; f < la_zone_VEF.valeur().zone().nb_faces_elem(0); ++f)
            {
              int face = la_zone_VEF.valeur().elem_faces(elem, f);
              if(!faces_ssz(face))
                {
                  ci_(face) = val_ci_cathode(face);
                  // necessaire (source*porosite_surf(num_face));
                  faces_ssz(face) = 1;		// marquer comme deja traite
                }
            }
        }
    }

  // mettre a jour ch_S
  DoubleTab& val_S = ch_S_.valeurs();
  if(CL_a_.non_nul())
    {
      IntTab faces_ssz;	// faces belong to the sous_zone -> flag = 1, if not, flag = 0
      la_zone_VEF.valeur().creer_tableau_faces(faces_ssz);
      faces_ssz = 0;	// init with no flag (all faces are unchecked)

      for (int poly = 0; poly < CL_a_.valeur().nb_elem_tot(); poly++)
        {
          int elem = CL_a_.valeur()(poly);
          for (int f = 0; f < la_zone_VEF.valeur().zone().nb_faces_elem(0); ++f)
            {
              int face = la_zone_VEF.valeur().elem_faces(elem, f);
              if(!faces_ssz(face))
                {
                  val_S = eval_f(diffu_(face), C_(face), ci_(face), T_(face));
                  faces_ssz(face) = 1;
                }
            }
        }
    }

  if(CL_c_.non_nul())
    {
      IntTab faces_ssz;	// faces belong to the sous_zone -> flag = 1, if not, flag = 0
      la_zone_VEF.valeur().creer_tableau_faces(faces_ssz);
      faces_ssz = 0;	// init with no flag (all faces are unchecked)
      for (int poly = 0; poly < CL_c_.valeur().nb_elem_tot(); poly++)
        {
          int elem = CL_c_.valeur()(poly);
          for (int f = 0; f < la_zone_VEF.valeur().zone().nb_faces_elem(0); ++f)
            {
              int face = la_zone_VEF.valeur().elem_faces(elem, f);
              if(!faces_ssz(face))
                {
                  val_S = eval_f(diffu_(face), C_(face), ci_(face), T_(face));
                  faces_ssz(face) = 1;		// marquer comme deja traite
                }
            }
        }
    }

//  if(ch_ir_.non_nul())
//    {
//      ch_ir_.valeur().valeur_aux(xv, ir_);
//    }
  //  else
  //    {
  //      ir_ = 0.;
  //    }

//  if(ch_ip_.non_nul())
//    {
//      ch_ip_.valeur().valeur_aux(xv, ip_);
//    }
  //  else
  //    {
  //      ip_ = 0.;
  //    }

//  if(ch_op_.non_nul())
//    {
//      ch_op_.valeur().valeur_aux(xv, op_);
//    }
  //  else
  //    {
  //      op_ = 0.;
  //    }

}

void Source_Term_pemfc_VEF_P1NC::remplir_volumes()
{
  volumes_.ref(la_zone_VEF.valeur().volumes_entrelaces());
}
