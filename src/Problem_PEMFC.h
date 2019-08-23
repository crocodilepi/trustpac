#ifndef Problem_PEMFC_included
#define Problem_PEMFC_included

#include <Probleme_Couple.h>
class Problem_PEMFC : Probleme_Couple
{
  Declare_instanciable(Problem_PEMFC);
public:
  void associer_sch_tps_base (Schema_Temps_base &);
};
#endif
