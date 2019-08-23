#include <Schema_Temps_PEMFC_AME.h>
#include <Matrice_Bloc.h>
#include <Matrice_Morse.h>
#include <Param.h>
#include <Probleme_Couple.h>
#include <Noms.h>
#include <Probleme_base.h>
#include <Equation_base.h>
#include <Matrice_Dense.h>
#include <ConstDoubleTab_parts.h>
#include <Scatter.h>
#include <Zone_VEF.h>
#include <Solv_GCP.h>
#include <Loi_Fermeture_transport_ionique.h>
#include <EcrFicCollecte.h>
Implemente_instanciable(Schema_Temps_PEMFC_AME,"Schema_Temps_PEMFC_AME",Schema_Temps_base);

Sortie& Schema_Temps_PEMFC_AME::printOn(Sortie& s) const
{
  return  Schema_Temps_base::printOn(s);
}

Entree& Schema_Temps_PEMFC_AME::readOn(Entree& s)
{
  already_initialized_ = false ;
  Schema_Temps_base::readOn(s);
  return s;
}

void Schema_Temps_PEMFC_AME::set_param(Param& param)
{
  // ATTENTION, le solveur GCP de trust ne fonctionne pas car la structures
  // parallele de l'inconnue est en plusieurs blocs.
  // En revanche les solveurs petsc semblent fonctionner.
  param.ajouter("solverPhiPsi", & solverPhiPsi_, Param::REQUIRED);
  param.ajouter("name_coupled_problem", & name_pb_phipsi_, Param::REQUIRED);
  param.ajouter_flag("dump_matrix", & dumpMatrix_);
  Schema_Temps_base::set_param(param);
}

void copierMatrice(const Matrice_Morse & m, Matrice_Dense & m2)
{
  int ni = m.nb_lignes();
  int nj = m.nb_colonnes();
  m2.dimensionner(ni,nj);
  for (int i = 0; i < ni; i++)
    for (int j = 0; j < nj; j++) {
      m2.set_coefficient(i,j, m(i,j));
    }
}

// initialization step to do "typer" only once and not at each time step
// also : avoid to resize objects at each time step
bool Schema_Temps_PEMFC_AME::initialization( Probleme_Couple& pbc )
{

  if ( already_initialized_ )
    {
      return true;
    }
  else
    {

      int i;
      nb_pbs_ = 2;
      already_initialized_ = true ;
      the_matrix_.typer("Matrice_Bloc");

      Matrice_Bloc & matrix = ref_cast(Matrice_Bloc, the_matrix_.valeur());
      matrix.dimensionner( nb_pbs_ , nb_pbs_ );

      pb_names_.dimensionner( nb_pbs_ );
      pb_names_[0] = "pb_psi";
      pb_names_[1] = "pb_phi";


      for (i = 0; i < nb_pbs_; i++)
	{
	  Probleme_base & pb = ref_cast(Probleme_base, pbc.probleme(pb_names_[i]));
	  Equation_base & eq = pb.equation(0);
	  Matrice & mat = matrix.get_bloc(i,i);
	  mat.typer("Matrice_Morse");
	  Matrice_Morse & matmorse = ref_cast(Matrice_Morse, mat.valeur());
	  eq.dimensionner_matrice(matmorse);
	}


      // Create the parallel data structure of the entire rhs and unknowns
      // concatenation of the two array structures:
      for (i = 0; i < nb_pbs_; i++)
	{
	  Probleme_base & pb = ref_cast(Probleme_base, pbc.probleme(pb_names_[i]));
	  Equation_base & eq = pb.equation(0);
	  mdc_.add_part(eq.inconnue().valeur().valeurs().get_md_vector());
	}

      // In order to use this structure descriptor to build vectors,
      // we have to create a "constant" MD_Vector that cannot anymore be changed
      // this is done by creating a copy of it that is accessible via the MD_Vector pointer structure
      md_.copy(mdc_);

      // Create the global vectors with the created structure:
      rhs_global_.resize( md_.valeur().get_nb_items_tot() );
      unknown_global_.resize( md_.valeur().get_nb_items_tot() );
      rhs_global_.set_md_vector( md_ );
      unknown_global_.set_md_vector( md_ );

      const int& npsi = matrix.get_bloc(0,0).valeur().nb_lignes();
      const int& nphi = matrix.get_bloc(1,1).valeur().nb_lignes();

      for (int ibloc = 0; ibloc < nb_pbs_; ibloc++)
	{
	  const int jbloc = 1-ibloc;
	  matrix.get_bloc(ibloc,jbloc).typer("Matrice_Morse");
	  Matrice_Morse &m = ref_cast(Matrice_Morse,matrix.get_bloc(ibloc,jbloc).valeur());
	  if (ibloc==0)
	    m.dimensionner(npsi,nphi,npsi);
	  else
	    m.dimensionner(nphi,npsi,nphi);
	}

    }
  return true;
}


void Schema_Temps_PEMFC_AME::ComputePhiPsi( Probleme_Couple& pbc )
{

  Matrice_Bloc & matrix = ref_cast(Matrice_Bloc, the_matrix_.valeur());
  int i,j;
  // Create a "pointer" array that points to the 2 parts of unknown_global:
  // (no storage allocated, just pointers inside the unknown_global array):
  DoubleTab_parts rhs_parts( rhs_global_ );
  DoubleTab_parts unknown_parts( unknown_global_ );

  for (i = 0; i < nb_pbs_; i++) {
    Probleme_base & pb = ref_cast(Probleme_base, pbc.probleme(pb_names_[i]));
    Equation_base & eq = pb.equation(0);
    unknown_parts[i] = eq.inconnue(); // copy the structure of the unknown
  }

  // Get access to the closure law:
  Objet_U & obj = interprete().objet("loi_trans_ion");
  Loi_Fermeture_transport_ionique & clos_law = ref_cast(Loi_Fermeture_transport_ionique,obj);

  const ArrOfInt & face_index_psi = clos_law.face_index_psi();
  const ArrOfInt & face_index_phi = clos_law.face_index_phi();
  //const ArrOfInt & face_index_T = clos_law.face_index_T();
  //const ArrOfInt & face_is_anode = clos_law.face_is_anode();
  const int nPhiPsiFaces = face_index_psi.size_array();
  // Compute the volume associated with common faces.
  DoubleVect volumeFaces(nPhiPsiFaces);
  DoubleVect coefCoupling(nPhiPsiFaces);
  DoubleVect rhsCoupling(nPhiPsiFaces);
  DoubleTab invmasse_psi = unknown_parts[0]; // create a copy
  invmasse_psi = 1.;
  ref_cast(Probleme_base, pbc.probleme(pb_names_[0])).equation(0).solv_masse().appliquer(invmasse_psi);
  DoubleTab invmasse_phi = unknown_parts[1]; // create a copy
  invmasse_phi = 1.;
  ref_cast(Probleme_base, pbc.probleme(pb_names_[1])).equation(0).solv_masse().appliquer(invmasse_phi);

  DoubleTab ir, DirDpsi;

  // ****************************************
  // Loop to converge
  int iteration = 0;
  while (true) {
    ArrOfDouble ir_old = ir;
    clos_law.compute_ir_DirDpsi(unknown_parts[0], unknown_parts[1], ir, DirDpsi);
    if (ir_old.size_array() == ir.size_array()) {
      double max_delta = 0;
      double max_val = 0.;
      for (i = 0; i < nPhiPsiFaces; i++) {
        max_delta = max(max_delta, fabs(ir[i]-ir_old[i]));
        max_val = max(max_val, fabs(ir[i]));
      }
      max_delta = mp_max(max_delta);
      max_val = mp_max(max_val);
      Cerr << "Phi_psi convergence iteration " << iteration;
      Cerr << " delta=" << max_delta << " / max=" << max_val << endl;
      iteration++;
      double convergence_threshold = 1e-5;
      if (max_val == 0. or max_delta/max_val < convergence_threshold)
        break;
    }

    for (i = 0; i < nb_pbs_; i++) {
      Probleme_base & pb = ref_cast(Probleme_base, pbc.probleme(pb_names_[i]));
      Equation_base & eq = pb.equation(0);
      Matrice & mat = matrix.get_bloc(i,i);
      Matrice_Morse & matmorse = ref_cast(Matrice_Morse, mat.valeur());
      matmorse.get_set_coeff() = 0.;
      rhs_parts[i] = 0.;
      eq.assembler(matmorse,eq.inconnue(),rhs_parts[i]);
      // For each operator: call contribuer_a_avec(inco, matrix)
      // and operator.ajouter(rhs)
      // then sources.ajouter(rhs)
      // and sources().contribuer_a_avec(inco,matrice);
      // and finally rhs = rhs + matrix *inco
      //   eq.assembler(matrix(i,i), eq.inconnue(), rhs[i]);
      // so
      //   (matrix*x)[i] = - "integral on element i "(-div(k grad x))
      //   rhs = matrix * x(time n) + ??????
    }


    for (i = 0; i < nPhiPsiFaces; i++) {
      int ipsi = face_index_psi[i];
      int iphi = face_index_phi[i];
      // if face is on the boundary of one domain, invmass is higher,
      // take the highest of the 2 values
      volumeFaces[i] = 1./(max(invmasse_psi[ipsi], invmasse_phi[iphi]));
      coefCoupling[i] = DirDpsi[i] * volumeFaces[i]; // value in equations with laplacian of psi
      rhsCoupling[i] = (- ir[i] + DirDpsi[i] * (unknown_parts[0][ipsi] - unknown_parts[1][iphi])) * volumeFaces[i];
    }

    //Cerr << "volumeFaces" << endl << volumeFaces << endl;
    // Add coupling via d/dphi:
    for (int ibloc = 0; ibloc < nb_pbs_; ibloc++) {
      int jbloc = 1-ibloc;
      Matrice_Morse &M = ref_cast(Matrice_Morse,matrix.get_bloc(ibloc,ibloc).valeur());
      Matrice_Morse &m = ref_cast(Matrice_Morse,matrix.get_bloc(ibloc,jbloc).valeur());
      IntVect & tab1 = m.get_set_tab1();
      IntVect & tab2 = m.get_set_tab2();
      DoubleVect & coeff = m.get_set_coeff();
      coeff = 0.;
      int nlines = m.nb_lignes();
      for (i = 0; i < nlines+1; i++)
        tab1[i] = i+1; // fortran index in tab2
      for (i = 0; i < nlines; i++)
        tab2[i] = 1; // fortran column index
      for (i = 0; i < nPhiPsiFaces; i++) {
        int ipsi = face_index_psi[i];
        int iphi = face_index_phi[i];
        if (ibloc == 0) {
          tab2[ipsi] = iphi+1; // fortran index :(
          coeff[ipsi] = - coefCoupling[i];
          M(ipsi,ipsi) += coefCoupling[i];
          rhs_parts[0][ipsi] += rhsCoupling[i];
        } else {
          tab2[iphi] = ipsi+1; // fortran index :(
          coeff[iphi] = - coefCoupling[i];
          M(iphi,iphi) += coefCoupling[i];
          rhs_parts[1][iphi] -= rhsCoupling[i];
        }
      }
    }
    // Normalize matrix and RHS:
    double norm = rhs_global_.mp_norme_vect();
    Cerr << "Normalize matrix and rhs, norm=" << norm << endl;
    if (norm > 1e-12) {
            rhs_global_ *= (1./norm);
            for (j = 0; j < 2; j++) {
                    for (i = 0; i < 2; i++) {
                            Matrice_Morse &m = ref_cast(Matrice_Morse,matrix.get_bloc(i,j).valeur());
                            m.get_set_coeff() *= (1./norm);
                    }
            }
    }

    if (dumpMatrix_) {
      for (j = 0; j < 2; j++) {
        for (i = 0; i < 2; i++) {
          Matrice_Dense mm;
          Matrice_Morse &m = ref_cast(Matrice_Morse,matrix.get_bloc(i,j).valeur());
          copierMatrice(m,mm);
          EcrFicCollecte f(Nom("Mat")+Nom(i)+Nom(j)+Nom(".txt"));
          f << mm;
        }
      }
    }


    solverPhiPsi_.valeur().fixer_nouvelle_matrice(1);
    solverPhiPsi_.resoudre_systeme( matrix, rhs_global_, unknown_global_ );
  }

  // Copy solution back to the equation unknowns:
  for (i = 0; i < nb_pbs_; i++) {
    Probleme_base & pb = ref_cast(Probleme_base, pbc.probleme(pb_names_[i]));
    Equation_base & eq = pb.equation(0);
    eq.inconnue().futur() = unknown_parts[i]; // Note : inconnue().valeurs() donne ce que serait inconnue().present()
  }
  // Add coupling via d/dphi:
  //matrix(0,1).typer("Matrice_Morse");
  // Loop on elements of the CL anodes and cathodes

}

bool Schema_Temps_PEMFC_AME::iterateTimeStep(bool &converged)
{
  // for(i=0; i<pbc.nb_problemes(); i++)
  // // recopie le present dans le passe ????
  // Initialiser_Champs(ref_cast(Probleme_base,pbc.probleme(i)));
  if (name_pb_phipsi_ == "disabled")
    return true;

  Probleme_Couple & pbc = ref_cast(Probleme_Couple, interprete().objet(name_pb_phipsi_));
  if( ! initialization( pbc ) )
    {
      Cerr << "Error in Schema_Temps_PEMFC_AME::faire_un_pas_de_temps_pb_couple " << finl;
      Cerr << "A problem occured during the initialization step"<<finl;
      Cerr << "Abording..."<<finl;
      Process::abort( );
    }

  // use newton algorithm to compute coupled psi/phi
  ComputePhiPsi( pbc );

  converged = true;
  return true;
}
