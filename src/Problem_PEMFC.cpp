#include <Problem_PEMFC.h>
#include <Schema_Temps_PEMFC_AME.h>

Implemente_instanciable(Problem_PEMFC,"Problem_PEMFC",Probleme_Couple);

Sortie& Problem_PEMFC::printOn(Sortie& s) const
{
  return  Probleme_Couple::printOn(s);
}

Entree& Problem_PEMFC::readOn(Entree& s)
{
  Probleme_Couple::readOn(s);
  return s;
}

void Problem_PEMFC::associer_sch_tps_base (Schema_Temps_base & sch)
{
  Cerr << "Error in Problem_PEMFC::associer_sch_tps_base:" << endl
      << " Problem_PEMFC must not be associated with a time scheme" << endl
      << " Each subproblem must have its own associated time scheme" << endl;
  exit();
}
