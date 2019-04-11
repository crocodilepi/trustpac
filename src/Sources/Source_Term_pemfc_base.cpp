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
// File      : Source_Term_pemfc_base.cpp
// Directory : $PEMFC_ROOT/src/Sources
//
/////////////////////////////////////////////////////////////////////////////

#include <Source_Term_pemfc_base.h>
#include <Probleme_base.h>

Implemente_base( Source_Term_pemfc_base, "Source_Term_pemfc_base", Source_base ) ;

Sortie& Source_Term_pemfc_base::printOn( Sortie& os ) const
{
  Source_base::printOn( os );
  return os;
}

Entree& Source_Term_pemfc_base::readOn( Entree& is )
{
	Source_base::readOn( is );
	Param param(que_suis_je());
	Nom nom_champ_diffu, nom_pb_T, nom_pb_c, nom_champ_T, nom_champ_c;
	param.ajouter("nom_espece",&nom_espece_,Param::REQUIRED);
	param.ajouter("epsilon", &epsilon_, Param::REQUIRED);
	param.ajouter("epsilon_ionomer", &epsilon_ionomer_, Param::REQUIRED);
	param.ajouter("gamma_CL", &gamma_CL_, Param::REQUIRED);
	param.ajouter("gamma_CL", &gamma_CL_, Param::REQUIRED);
	param.ajouter("champ_D", &nom_champ_diffu, Param::REQUIRED);
	param.ajouter("champ_T", &nom_champ_T, Param::REQUIRED);
	param.ajouter("pb_c", &nom_pb_c, Param::REQUIRED);
	param.ajouter("champ_c", &nom_champ_T, Param::REQUIRED);
	param.lire_avec_accolades(is);
	REF(Probleme_base) pb_T = ref_cast(Probleme_base,interprete().objet(nom_pb_T));
	REF(Probleme_base) pb_c = ref_cast(Probleme_base,interprete().objet(nom_pb_c));
	T_ = pb_T.valeur().get_champ(nom_champ_T);	// champ P0 if VDF, champ P1NC if VEF
	c_ = pb_c.valeur().get_champ(nom_champ_c);  // TO-DO verify
	Da_ = equation().probleme().get_champ(nom_champ_diffu);	// champ P0 if VDF, champ P1NC if VEF
	C_ = equation().inconnue();					// // champ P0 if VDF, champ P1NC if VEF
	thickness_ionomer_ = (1-epsilon_)*epsilon_ionomer_ / gamma_CL_;
	return is;
}

double Source_Term_pemfc_base::eval_f(double diffu, double Ci, double ci, double T) const
{
	double R = 8.314;
	double H; 			// H Henri constant (mol.m3/Pa)
	if (nom_espece_ == "H2") {
		H = 1. / 4.5e3;
	} else if (nom_espece_ == "02") {
		H = 1. / 1.33e5 * exp(666/T);
	} else if (nom_espece_ == "N2") {
		H = 6.4e-6 * exp(1300*(1./T - 1/298.15));
	}
	double Ceq;
	if (nom_espece_ == "H2O" || nom_espece_ == "vap") {
		double Psat = 1.5;	// Pa								// TO-DO:vefify this value!!!!!
		double activ = ci * R * T / Psat;
		double lambda_eq = 0.043+17.81*activ-39.85*activ*activ+36*activ*activ*activ;
		double CSO3 = 2036;
		Ceq = lambda_eq * CSO3;
	} else {
		Ceq = ci * R * T * H;
	}
	return diffu * gamma_CL_ / thickness_ionomer_ * (Ceq - Ci);
}

double Source_Term_pemfc_base::eval_derivee_f(double diffu) const
{
	// expression_derivee_par_rapport_inconnue
	return (- diffu * gamma_CL_ / thickness_ionomer_);
}

DoubleTab& Source_Term_pemfc_base::ajouter(DoubleTab& resu) const
{
	assert(resu.dimension(0)==volumes_.size());
	int size=resu.dimension(0);

	const DoubleTab& T=T_.valeur().valeurs();
	const DoubleTab& C=C_.valeur().valeurs();
	const DoubleTab& D=Da_.valeur().valeurs();
	const DoubleTab& c=c_.valeur().valeurs();

	// check
	assert(resu.dimension(0)==T.size());
	assert(resu.dimension(0)==C.size());
	assert(resu.dimension(0)==D.size());
	assert(resu.dimension(0)==c.size());

	for (int i=0; i<size; i++)
	{
	  resu(i)+= volumes_(i) * eval_f(D(i,0), C(i,0), c(i,0), T(i,0));
	}
	return resu;
}

DoubleTab& Source_Term_pemfc_base::calculer(DoubleTab& resu) const
{
	resu = 0;
	ajouter(resu);
	return resu;
}

void Source_Term_pemfc_base::mettre_a_jour(double temps)
{
	// update the involving field
	T_.valeur().mettre_a_jour(temps);
	C_.valeur().mettre_a_jour(temps);
	c_.valeur().mettre_a_jour(temps);
	Da_.valeur().mettre_a_jour(temps);
}

void Source_Term_pemfc_base::contribuer_a_avec(const DoubleTab& inco, Matrice_Morse& mat) const
{
	int size=inco.dimension(0);
	const DoubleTab& D=Da_.valeur().valeurs();
	for (int i=0; i<size; i++)
	{
	  mat.coef(i,i)+=volumes_(i) * eval_derivee_f(D(i,0));
	}
}

void Source_Term_pemfc_base::associer_pb(const Probleme_base& pb)
{

}
