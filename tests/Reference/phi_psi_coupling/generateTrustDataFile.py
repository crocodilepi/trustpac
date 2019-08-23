# This script takes the output of the step1 script (trustReadGeometry.data)
# and generates the phi_psi_coupling.data file.
# Takes no arguments, works in current directory.
class Generator:
    def __init__(self, datafilename, params):
        self.P = params
        self.f = open(datafilename,"w")
    def close(self):
        self.f.close()
    def __del__(self):
        self.close()
    def pr(self,string):
        newstring = string
        for name,value in self.P.items():
            s = '$'+name
            newstring = newstring.replace(s,str(value))
        self.f.write(newstring)
    def insertFile(self,file):
        with open(file) as f:
            lines = f.readlines()
            self.f.writelines(lines)
    def probes(self,s):
        pass

# Define parameters for this simulation
def Diag33(a00,a11=None,a22=None):
    if a11==None:
        a11=a00
    if a22==None:
        a22=a11
    return "%f 0 0  0 %f 0  0 0 %f"%(a00,a11,a22)

P={}
P["sigma_x_GDLa"]=100
P["sigma_y_GDLa"]=100
P["sigma_x_MPLa"]=100
P["sigma_y_MPLa"]=100
P["sigma_CLa"]=20
P["CathodeVoltage"]=0.3
P["AnodeVoltage"]=0.0
def AnodeCathodeFunc(anodeVal,cathodeVal):
    # generates a champ_fonc_xyzt with two different values in anode and cathode regions
    return " champ_fonc_txyz dom_phi 1 (z>(-220e-6))*(%s)+(z<(-220e-6))*(%s)"%(str(anodeVal),str(cathodeVal))
P["i0_field_override"] = "i0_field_override"+AnodeCathodeFunc(anodeVal=30,cathodeVal=5e-6)
P["Erev_field_override"] = "Erev_field_override"+AnodeCathodeFunc(anodeVal=0.032,cathodeVal=1.14)

P["MatrixSigma_CLa"]=Diag33(P["sigma_CLa"])
P["MatrixSigma_MPLa"]=Diag33(P["sigma_x_MPLa"],P["sigma_x_MPLa"],P["sigma_y_MPLa"])
P["MatrixSigma_GDLa"]=Diag33(P["sigma_x_GDLa"],P["sigma_x_GDLa"],P["sigma_y_GDLa"])

g = Generator("phi_psi_coupling.data",params=P)
g.pr("dimension 3\n")
g.pr("# PARALLEL OK 2 #\n")
g.insertFile("trustReadGeometry.data")
g.pr("""
Pb_Conduction pb_phi
Pb_Conduction pb_psi
Pb_Conduction pb_T
VEFPrep1B dis

Schema_Temps_PEMFC_AME sch_pemfc_ame
Read sch_pemfc_ame
{
    tinit 0.
    tmax 10
    nb_pas_dt_max 1
    dt_min 1e-15
    dt_max 1
    dt_impr 0.001
    dt_sauv 400.
    # solver gmres { diag seuil 1e-8 nb_it_max 50 impr } # # does not seem to work well on this system #
    solverPhiPsi petsc cholesky { } # Cholesky is the fastest for sequential simulation and small meshes, 1s @ 60000 ddl #
    # solverPhiPsi petsc gcp { precond spai { } seuil 1e-9 impr } # # 3.9s@60000ddl #
    # solverPhiPsi petsc gcp { precond diag { } seuil 1e-9 impr } # # 3.3s@60000ddl #
    # solverPhiPsi petsc gcp { precond boomeramg { } seuil 1e-9 impr } # # 3.2s@60000ddl #
    facsec 1e4
    name_coupled_problem pbc
}
Schema_Temps_PEMFC_AME sch_disabled
Read sch_disabled
{
    tinit 0.
    tmax 10
    nb_pas_dt_max 1
    dt_min 1e-15
    dt_max 1
    dt_impr 0.001
    dt_sauv 400.
    solverPhiPsi petsc cholesky { }
    facsec 1e4
    name_coupled_problem disabled
}
# Empty time scheme for T #
Schema_Euler_Explicite sch_T { }

# proprietaire pour le transport ionique #
Solide prop_phi
Read prop_phi
{
    rho     champ_uniforme 1 1.
    Cp      champ_uniforme 1 1.
      # this is a uniform isotropic conduction field. #
      # In 2D: diffusion matrix is [[ a b ] [b a ]] => encoded [a b b a] #
      # In 3D: diffusion matrix is
        [[ a b c ]
         [ d e f ]
         [ g h i  ]] => encoded [a b c d e f g h i ] #
    lambda  champ_uniforme_morceaux dom_phi 9 {
        defaut
            6. 0  0
            0  6 0
            0  0  6
        dom_phi_CLa
            0.8 0  0
            0  0.8 0
            0  0  0.8
        dom_phi_CLc
            0.8  0   0
            0   0.8  0
            0    0   0.8  }
}

# couplage de conductivite kappa, courant Ii = kappa.gradT, laplacien Di = div(kappa.gradT) avec T,  C_H2O #
Loi_Fermeture_transport_ionique loi_trans_ion
associer pb_phi loi_trans_ion

# Electronic transport properties #
Solide prop_psi
Read prop_psi
{
    rho     champ_uniforme 1 1.
    Cp      champ_uniforme 1 1.
    lambda  champ_uniforme_morceaux dom_psi 9 {
        defaut      $MatrixSigma_CLa
        dom_psi_CLc $MatrixSigma_CLa
        dom_psi_GDLa $MatrixSigma_GDLa
        dom_psi_GDLc $MatrixSigma_GDLa
        dom_psi_MPLa $MatrixSigma_MPLa
        dom_psi_MPLc $MatrixSigma_MPLa
    }
}

# dummy properties #
Solide prop_T
Read prop_T
{
    rho     champ_uniforme 1 1.
    Cp      champ_uniforme 1 1.
    lambda  champ_uniforme 1 1.
}

Associate pb_phi dom_phi
Associate pb_phi prop_phi

Associate pb_psi dom_psi
Associate pb_psi prop_psi

Associate pb_T dom_T
Associate pb_T prop_T

Problem_PEMFC pbc
Associate pbc pb_phi
Associate pbc pb_psi
Associate pbc pb_T

Associate pb_phi sch_pemfc_ame
Associate pb_psi sch_disabled
Associate pb_T   sch_T

Discretize pbc dis

Read loi_trans_ion
{
   T_0 		353.15		# default value #
   CSO3 	2036.
   por_naf champ_uniforme_morceaux dom_phi 1 { defaut 0.0 dom_phi_CLa 0.47 dom_phi_CLc 0.47 }	# porosite de Nafion : CL 0.47, MB 0 #
   eps_naf champ_uniforme_morceaux dom_phi 1 { defaut 1. dom_phi_CLa 0.2 dom_phi_CLc 0.2	}	# ionomer proportionnel de Nafion: CA 0.2, MB 1 #
   tor_naf champ_uniforme_morceaux dom_phi 1 { defaut 1. dom_phi_CLa 1.5 dom_phi_CLc 1.3	}		# tortuosite de Nafion #
   nom_ssz_CLc dom_phi_CLc
   nom_ssz_CLa dom_phi_CLa
   nom_pb_T 	pb_T
   nom_champ_T 	temperature
   nom_pb_C_H2 none
   nom_champ_C_H2 none
   nom_pb_C_O2 none
   nom_champ_C_O2 none
   nom_pb_C_H2O none
   nom_champ_C_H2O none
   nom_pb_psi pb_psi
   nom_champ_psi temperature
   nom_pb_phi pb_phi
   nom_champ_phi temperature
   $i0_field_override
   $Erev_field_override
}
# Temperature problem is only here to provide a temperature field for the other equ. #
Read pb_T
{
	Conduction
	{
		equation_non_resolue 1
		diffusion { }
		initial_conditions { temperature Champ_Uniforme 1 357.0 }
		boundary_conditions {
            I_Periodic periodique
			defaultBC paroi_adiabatique
			I_GDLa_CHa paroi_adiabatique
			I_GDLc_CHc paroi_adiabatique
			I_CHa_Bipo paroi_adiabatique
			I_CHc_Bipo paroi_adiabatique
            I_Bipo_Cooling paroi_adiabatique
		}
	}
}

Read pb_phi
{
	Conduction
	{
	    equation_non_resolue 0
		diffusion { Tensor }
		initial_conditions { temperature Champ_Uniforme 1 -0.2 }
		boundary_conditions { defaultBC paroi_adiabatique }
    }
	Post_processing { fields dt_post 1e-8 { TEMPERATURE ELEM IR_ ELEM } }
    # Post_processing  { format lata fields dt_post 1e-8 {
      TEMPERATURE ELEM TEMPERATURE FACES
      ir_ FACES
    } } #
}

Read pb_psi
{
	Conduction {
        equation_non_resolue 0
        diffusion { Tensor }
        initial_conditions {
        	temperature champ_uniforme_morceaux dom_psi 1 {
                defaut $AnodeVoltage
                dom_psi_CLc $CathodeVoltage
                dom_psi_MPLc $CathodeVoltage
                dom_psi_GDLc $CathodeVoltage }
        }
        boundary_conditions {
            I_Periodic periodique
            defaultBC paroi_adiabatique
            I_GDLa_Bipo  paroi_echange_externe_impose
               t_ext champ_front_uniforme 1 $AnodeVoltage
               h_imp champ_front_uniforme 1 1e9
            I_GDLc_Bipo  paroi_echange_externe_impose
               t_ext champ_front_uniforme 1 $CathodeVoltage
               h_imp champ_front_uniforme 1 1e9
        }
	}
	Post_processing	{ fields dt_post 1e-8 { TEMPERATURE ELEM } }
    # Post_processing  { format lata fields dt_post 1e-8 { TEMPERATURE ELEM TEMPERATURE FACES } } #
}

Solve pbc
End
""")
