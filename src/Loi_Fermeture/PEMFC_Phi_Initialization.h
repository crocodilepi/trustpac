/////////////////////////////////////////////////////////////////////////////
//
// File      : PEMFC_Phi_Initialization.h
// Directory : $PEMFC_ROOT/src/Loi_Fermeture
//
/////////////////////////////////////////////////////////////////////////////

#ifndef PEMFC_Phi_Initialization_included
#define PEMFC_Phi_Initialization_included

#include <Loi_Fermeture_base.h>
#include <Param.h>
#include <Discretisation_base.h>
#include <Champ_Fonc.h>
#include <Ref_Champ_base.h>
#include <Ref_Champ_Inc.h>
#include <Ref_Operateur_base.h>
#include <DoubleTab.h>
#include <Ref_Sous_Zone.h>


/////////////////////////////////////////////////////////////////////////////
//
// .DESCRIPTION : class PEMFC_Phi_Initialization
// cette clase permet generer le champ d'initialisation de phi
// la relation suivante doit etre verifiee :
// phi = Erev - psi
// <Description of class PEMFC_Phi_Initialization>
//
/////////////////////////////////////////////////////////////////////////////

class PEMFC_Phi_Initialization : public Loi_Fermeture_base
{

  Declare_instanciable( PEMFC_Phi_Initialization ) ;

public :
  void set_param( Param& param );
  void compute_phi_init( void ) ;

protected :

  Nom nom_ssz_CLa_;
  Nom nom_ssz_CLc_;
  Nom nom_pb_psi_ ;
  Nom nom_loi_ ;
  Nom nom_pb_phi_;

  Champ_Fonc ch_phi_init_ ;

  REF(Sous_Zone) CL_a_;
  REF(Sous_Zone) CL_c_;

};

#endif /* PEMFC_Phi_Initialization_included */
