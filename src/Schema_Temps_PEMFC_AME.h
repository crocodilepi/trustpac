#ifndef  Schema_Temps_PEMFC_AME_included
#define Schema_Temps_PEMFC_AME_included
#include <Schema_Temps_base.h>
#include <Matrice.h>
#include <Noms.h>
#include <MD_Vector_composite.h>
#include <SolveurSys.h>

class Probleme_Couple;

class Schema_Temps_PEMFC_AME : Schema_Temps_base
{
  Declare_instanciable(Schema_Temps_PEMFC_AME);
public:
  void set_param(Param& param);
  int nb_valeurs_temporelles() const { return 3; }
  int nb_valeurs_futures() const { return 1; }
  double temps_futur(int) const { return temps_courant()+pas_de_temps(); }
  double temps_defaut() const { return temps_courant()+pas_de_temps(); }
  bool iterateTimeStep(bool &converged);
  void completer() {}
  void ComputePhiPsi( Probleme_Couple& pbc );
  int faire_un_pas_de_temps_eqn_base(Equation_base&) {exit(-1);return 0;}

protected:
  bool initialization( Probleme_Couple& pbc );

  Nom name_pb_phipsi_;
  SolveurSys solverPhiPsi_; // matrix solver for the coupled phi/psi problem
  bool already_initialized_;
  Matrice the_matrix_;

  Noms pb_names_;

  MD_Vector_composite mdc_;
  MD_Vector md_;

  DoubleVect rhs_global_ ;
  DoubleVect unknown_global_ ;

  int nb_pbs_ ;
  int dumpMatrix_;
};
#endif
