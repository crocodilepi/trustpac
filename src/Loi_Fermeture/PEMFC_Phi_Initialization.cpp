/////////////////////////////////////////////////////////////////////////////
//
// File      : PEMFC_Phi_Initialization.cpp
// Directory : $PEMFC_ROOT/src/Loi_Fermeture
//
/////////////////////////////////////////////////////////////////////////////

#include <PEMFC_Phi_Initialization.h>
#include <Probleme_base.h>
#include <Interprete.h>
#include <Domaine.h>
#include <Discret_Thermique.h>
#include <Equation_base.h>
#include <Zone_VF.h>

Implemente_instanciable( PEMFC_Phi_Initialization, "PEMFC_Phi_Initialization", Loi_Fermeture_base ) ;

Sortie& PEMFC_Phi_Initialization::printOn( Sortie& os ) const
{
  Loi_Fermeture_base::printOn( os );
  return os;
}

Entree& PEMFC_Phi_Initialization::readOn( Entree& is )
{
  Loi_Fermeture_base::readOn( is );
  compute_phi_init( );
  return is;
}

void PEMFC_Phi_Initialization::set_param(Param& param)
{
  param.ajouter("nom_pb_psi", &nom_pb_psi_, Param::REQUIRED);
  param.ajouter("nom_loi", &nom_loi_, Param::REQUIRED);
  param.ajouter("nom_pb_phi", &nom_pb_phi_, Param::REQUIRED);
  param.ajouter("nom_ssz_CLc", &nom_ssz_CLc_, Param::REQUIRED);
  param.ajouter("nom_ssz_CLa", &nom_ssz_CLa_, Param::REQUIRED);
}

void PEMFC_Phi_Initialization::compute_phi_init(  )
{
  const Probleme_base& pb_phi =  ref_cast( Probleme_base, interprete( ).objet( nom_pb_phi_ ) );
  const Probleme_base& pb_psi =  ref_cast( Probleme_base, interprete( ).objet( nom_pb_psi_ ) );
  const Equation_base& equation = pb_phi.get_equation_by_name( "Conduction" );
  const Discretisation_base& dis = equation.discretisation( );

  dis.discretiser_champ("champ_elem", equation.zone_dis( ).valeur( ), "phi_init", "unit", 1 ,0. , ch_phi_init_ );
  champs_compris_.ajoute_champ( ch_phi_init_ );

  const Domaine& dom = pb_phi.domaine();
  CL_a_ = dom.ss_zone( nom_ssz_CLa_ );
  CL_c_ = dom.ss_zone( nom_ssz_CLc_ );

  Loi_Fermeture_base& loi = ref_cast( Loi_Fermeture_base, interprete( ).objet( nom_loi_ ) );
  loi.mettre_a_jour( 0.0 );

  const DoubleTab& Erev = loi.get_champ( "Erev" ).valeurs( );
  // const DoubleTab& PsiP1NC  = pb_psi.get_champ("temperature").valeurs( );

  const Zone_VF& zvf = ref_cast( Zone_VF, equation.zone_dis().valeur( ) );
  const DoubleTab& xp= zvf.xp( );

  DoubleTab PsiP0;
  equation.zone_dis( ).zone( ).creer_tableau_elements( PsiP0 );
  pb_psi.get_champ("temperature").valeur_aux( xp, PsiP0 );

  DoubleTab& phi_init = ch_phi_init_.valeurs( );
  for (int poly = 0; poly < CL_a_.valeur().nb_elem_tot(); poly++)
    {
      const int& elem = CL_a_.valeur( )( poly );
      phi_init( elem ) = Erev( elem ) - PsiP0( elem ) ;
    }
  for (int poly = 0; poly < CL_c_.valeur().nb_elem_tot(); poly++)
    {
      const int& elem = CL_c_.valeur( )( poly );
      phi_init( elem ) = Erev( elem ) - PsiP0( elem ) ;
    }
}


