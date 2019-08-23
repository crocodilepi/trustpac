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

#include <Equation_AME.h>
#if 0
Implemente_instanciable_sans_constructeur(Equation_AME,"Equation_AME",Equation_base);

Conduction::Conduction()
{
  // todo: declare fields that can be accessed from outside
  // champs_compris_.ajoute_nom_compris("temperature_paroi");
}

Sortie& Conduction::printOn(Sortie& s ) const
{
  return s << que_suis_je() << "\n";
}

// This method is called before the readOn method, we declare here
//  the entries to be read in the .data file, which add to the default entries
void Conduction::set_param(Param& param)
{
  // Call the base class method to define base parameters
  Equation_base::set_param(param);
}

// This method is called during the .data file parsing. It reads what is inside
//  equation_AME { } block .
Entree& Conduction::readOn(Entree& is )
{
  // call the default reader
  Equation_base::readOn(is);
  return is;
}

int Conduction::nombre_d_operateurs() const
{
  return 1;
}

const Operateur& Conduction::operateur(int i) const
{
  assert (i==0);
  return operator_ame;
}

Operateur& Conduction::operateur(int i)
{
  assert (i==0);
  return operator_ame;
}


void Conduction::discretiser()
{
  const dis & dis = discretisation();
  Cerr << "Discretization of Equation AME" << finl;

  // Create the metadata for a composite vector composed of all unknowns:
  // This will finally depend on the discretisation type (vef or vdf)
  // Currently : hardcoded for vdf:
  {
    // Create a composite metadata for discretized fields:
    MD_Vector_composite md;
    const MD_Vector & at_nodes = zone().md_vector_elements();
    const MD_Vector & at_elements = zone().domaine().md_vector_sommets();
    const MD_Vector & at_faces = zone().domaine().md_vector_sommets();

#define ADD_VAR(name,where) md.add_part(where); unknowns_names_.ajouter(name)
    // Questions: where is N2 ?
    // "C_gas" = "C_N2 + C_H2_O2"
    // We are ready to declare all the unknowns of the problem:
    ADD_VAR("C_H2",at_nodes);
    ADD_VAR("C_O2",at_nodes);
    ADD_VAR("C_N2",at_nodes);
    ADD_VAR("C_vap",at_nodes); // vapour in gdl, dissolved in membrane
    ADD_VAR("P_gas",at_nodes);
    ADD_VAR("psi",at_nodes);  // electric potential
    ADD_VAR("T", at_nodes); // temperature !
    ADD_VAR("phi", at_nodes);
    ADD_VAR("psi", at_nodes);
    // Create a clean copy:
    md_vector_unknowns_.copy(md);
  }

  dis.temperature(schema_temps(), zone_dis(), la_temperature);
  champs_compris_.ajoute_champ(la_temperature);
  Equation_base::discretiser();
}
#endif
