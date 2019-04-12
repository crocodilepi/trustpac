/****************************************************************************
* Copyright (c) 2015, CEA
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
//////////////////////////////////////////////////////////////////////////////
//
// File:        Op_Diff_VEF_Face_Matricial.h
// Directory:   $TRUST_ROOT/src/VEF/Operateurs
// Version:     /main/24
//
//////////////////////////////////////////////////////////////////////////////


#ifndef Op_Diff_VEF_Face_Matricial_included
#define Op_Diff_VEF_Face_Matricial_included

#include <Op_Diff_VEF_Face.h>
#include <Champ_base.h>
#include <Champ_Don.h>
#include <Ref_Champ_Don.h>

class Param;

class Op_Diff_VEF_Face_Matricial : public Op_Diff_VEF_Face
{
  Declare_instanciable(Op_Diff_VEF_Face_Matricial);

public:

  virtual void set_param(Param& param);

  void associer_diffusivite(const Champ_base& );

  const Champ_base& diffusivite() const;

  DoubleTab& ajouter(const DoubleTab& ,  DoubleTab& ) const;
  DoubleTab& calculer(const DoubleTab& , DoubleTab& ) const;

  void ajouter_cas_multi_scalaire(const DoubleTab& inconnue,
                                  DoubleTab& resu, DoubleTab& flux_bords,
                                  const DoubleTab& nu,
                                  const Zone_Cl_VEF& zone_Cl_VEF,
                                  const Zone_VEF& zone_VEF,
                                  int nb_comp) const;

  void contribuer_a_avec(const DoubleTab&, Matrice_Morse&) const;
  void check_diffusivity( const int& nb_comp ) const;
  void mettre_a_jour(double);

protected :


  REF(Champ_base) diffusivite_;
  Champ_Don diffusivity_read_from_datafile_;

};


#endif


