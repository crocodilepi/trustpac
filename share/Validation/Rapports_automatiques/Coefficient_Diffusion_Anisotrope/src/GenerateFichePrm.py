import subprocess

#
# Creation de la fiche prm
#
with open( "fiche.prm", "w" ) as fiche :
    fiche.write("parameters \n")
    fiche.write("{\n")
    fiche.write("  title \"Anisotropic diffusion coefficient in VEF discretization\" \n")
    fiche.write("  description \"In this test case, we want to solve a 2D thermal conduction problem with \latex_($\\rho$=$C_p$=1 \latex_), in other words :\" \n")
    fiche.write("  description \"\latex_( $$ - \\nabla \cdot \left( K ~ \\nabla T \\right) = P $$ \latex_) with the anisotropic coefficients matrix given by\"\n")
    fiche.write("  description \"\latex_($$ K = \left( \\begin{matrix} K_{00} & K_{01} \\\\ K_{10} & K_{11} \end{matrix} \\right) $$\latex_)\"\n")
    fiche.write("  description \"We are looking for an analytic solution expressed as : \latex_( $$T_{A} = Ax^2+By^2+Cxy+Dx+Ey $$ \latex_) so the thermal power is given by \latex_( $$ P = -2~A~K_{00} - 2~B~K_{11} - C~( K_{01} + K_{10} ) $$ \latex_)\"\n")
    fiche.write("  author \"Stephane Veys\"\n")


cmd="./AddTestCases.sh"
subprocess.call( cmd , shell=True )

with open( "fiche.prm", "a" ) as fiche :
    fiche.write("}\n")



