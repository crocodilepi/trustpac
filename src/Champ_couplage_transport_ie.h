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
// File      : Champ_couplage_transport_ie.h
// Directory : $PEMFC_ROOT/src
//
/////////////////////////////////////////////////////////////////////////////

#ifndef Champ_couplage_transport_ie_included
#define Champ_couplage_transport_ie_included

#include <Champ_Fonc_P0_base.h>
#include <Ref_Probleme_base.h>
#include <Ref_Domaine.h>
#include <Ref_Sous_Zone.h>


/////////////////////////////////////////////////////////////////////////////
//
// .DESCRIPTION : class Champ_couplage_transport_ie
//
// <Description of class Champ_couplage_transport_ie>
//
/////////////////////////////////////////////////////////////////////////////

class Champ_couplage_transport_ie : public Champ_Fonc_P0_base
{

  Declare_instanciable( Champ_couplage_transport_ie ) ;

public :
  virtual void  mettre_a_jour(double un_temps);
  const REF(Domaine)& domaine() const;
  REF(Domaine)& domaine();
  const REF(Sous_Zone)& sous_zone() const;
  REF(Sous_Zone)& sous_zone();
  virtual void associer_zone_dis_base(const Zone_dis_base&) ;
  virtual const Zone_dis_base& zone_dis_base() const ;
  const DoubleTab& get_derivee() const
  {
    return derivee_;
  } ;
protected :

  REF(Zone_dis_base) zdb_;
  REF(Probleme_base) mon_pb_;
  REF(Domaine) mon_dom_;
  REF(Sous_Zone) mon_ssz_;
  //DoubleTab mon_barycentre_;
  int sens_;
  DoubleTab derivee_;

  REF(Probleme_base) pb_electrique_;
  REF(Domaine) dom_electrique_;
  REF(Sous_Zone) szz_electrique_;

  REF(Probleme_base) pb_ionique_;
  REF(Domaine) dom_ionique_;
  REF(Sous_Zone) szz_ionique_;

  DoubleTab bary_dom_electrique_, bary_dom_ionique_; // pour convertir les champs donnees vers champ P0
};

inline const REF(Domaine)& Champ_couplage_transport_ie::domaine() const
{
  return mon_dom_;
}

inline REF(Domaine)& Champ_couplage_transport_ie::domaine()
{
  return mon_dom_;
}


inline const REF(Sous_Zone)& Champ_couplage_transport_ie::sous_zone() const
{
  return mon_ssz_;
}

inline REF(Sous_Zone)& Champ_couplage_transport_ie::sous_zone()
{
  return mon_ssz_;
}

inline void Champ_couplage_transport_ie::associer_zone_dis_base(const Zone_dis_base& la_zone_dis_base)
{
  zdb_=la_zone_dis_base;
}

inline const Zone_dis_base& Champ_couplage_transport_ie::zone_dis_base() const
{
  return zdb_.valeur();
}

#endif /* Champ_couplage_transport_ie_included */
