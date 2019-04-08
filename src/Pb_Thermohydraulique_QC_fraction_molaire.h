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
/////////////////////////////////////////////////////////////////////////////
//
// File      : Pb_Thermohydraulique_QC_fraction_molaire.h
// Directory : $PEMFC_ROOT/src/Cas4
//
/////////////////////////////////////////////////////////////////////////////

#ifndef Pb_Thermohydraulique_QC_fraction_molaire_included
#define Pb_Thermohydraulique_QC_fraction_molaire_included

#include <Pb_Thermohydraulique_QC.h>
#include <Convection_Diffusion_Fraction_Molaire_QC.h>
#include <Convection_Concentration_Gaz.h>

/////////////////////////////////////////////////////////////////////////////
//
// .DESCRIPTION : class Pb_Thermohydraulique_QC_fraction_molaire
//
// <Description of class Pb_Thermohydraulique_QC_fraction_molaire>
//
/////////////////////////////////////////////////////////////////////////////

class Pb_Thermohydraulique_QC_fraction_molaire : public Pb_Thermohydraulique_QC
{

  Declare_instanciable( Pb_Thermohydraulique_QC_fraction_molaire ) ;

public :
  const Equation_base& equation(int) const ;
  Equation_base& equation(int);
  int nombre_d_equations() const;
  int verifier();
protected :
  Convection_Concentration_Gaz eq_cg_;
  Convection_Diffusion_Fraction_Molaire_QC eq_fraction_molaire_;
};

#endif /* Pb_Thermohydraulique_QC_fraction_molaire_included */
