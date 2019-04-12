import sys
import math

def ExtraCoeffs( A, B, C, D, E, K00, K01, K10, K11, case_number, facsec ) :

    #
    # definition de la solution analytique
    #

    analytical_solution=""
    eps = 1e-16
    A_printed=0
    B_printed=0
    C_printed=0
    D_printed=0
    E_printed=0
    if float( A ) > eps :
        A_printed=1
        analytical_solution = analytical_solution +  A +" x^2 "
    if float( B ) > eps  :
        B_printed=1
        if A_printed==1 :
            analytical_solution = analytical_solution+" + "
        analytical_solution = analytical_solution +  B +" y^2 "
    if float( C ) > eps :
        C_printed=1
        if A_printed==1 or B_printed==1 :
            analytical_solution = analytical_solution+" + "
        analytical_solution = analytical_solution +  C +" xy  "
    if float( D ) > eps :
        D_printed=1
        if A_printed==1 or B_printed==1 or C_printed==1 :
            analytical_solution = analytical_solution+" + "
        analytical_solution = analytical_solution +  D +" x  "
    if float( E ) > eps :
        E_printed=1
        if A_printed==1 or B_printed==1 or C_printed==1 or D_printed==1 :
            analytical_solution = analytical_solution+" + "
        analytical_solution = analytical_solution + E +" y  "


    # 
    # Analyse de la solution analytique
    # 

    solution_P1=""
    if float( A ) < eps and float( B ) < eps and float( C ) < eps :
        solution_P1="The analytical solution should be exactly capted by TRUST"
    else :
        solution_P1="The analytical solution can't be exactly capted by TRUST because an interpolation error is expected"

    #
    # calcul des ratios
    #
    fK00 = float( K00 )
    fK01 = float( K01 )
    fK10 = float( K10 )
    fK11 = float( K11 )

    ratio_extra_extra = 0
    min_ratio_extra_diag = 0
    max_ratio_extra_diag = 0
    max_coeff_extra = fK01 
    min_coeff_extra = fK10 
    max_coeff_diag = fK00 
    min_coeff_diag = fK11 

    if fK10 > fK01 :
        max_coeff_extra = fK10
        min_coeff_extra = fK01
    if fK11 > fK00 :
        max_coeff_diag = fK11
        min_coeff_diag = fK00

    ratio_extra_extra = max_coeff_extra / min_coeff_extra
    if min_coeff_diag > max_coeff_extra :
        min_ratio_extra_diag = min_coeff_diag / max_coeff_extra
    else :
        min_ratio_extra_diag = max_coeff_extra / min_coeff_diag
    if max_coeff_diag > min_coeff_extra :
        max_ratio_extra_diag = max_coeff_diag / min_coeff_extra
    else :
        max_ratio_extra_diag = min_coeff_extra / max_coeff_diag

    #
    # Calcul de la norme (simple) de la matrice
    # 
    matrix_norm = math.sqrt( fK00*fK00 + fK11*fK11 + fK01*fK01 + fK10*fK10 )
    
    #
    # Ecriture dans le fichier prm
    #

    with open( "../fiche.prm", "a" ) as fiche :
        fiche.write("chapter \n")
        fiche.write("{\n")
        fiche.write("  title \"Impact of extra-diagonals coefficients\" \n")
        fiche.write("  description \"Coefficients used : \n")
        fiche.write("  description \"\latex_($$ K = \left( \\begin{matrix} "+str(K00)+" & "+str(K01)+" \\\\ "+str(K10)+" & "+str(K11)+" \end{matrix} \\right) $$\latex_)\"\n")
        fiche.write("  description \"Ratio between extra diagonal coefficients : %.3e\"\n" %ratio_extra_extra )
        fiche.write("  description \"Minimal ratio between extra diagonal and diagonal coefficients : %.3e\"\n" %min_ratio_extra_diag )
        fiche.write("  description \"Maximal ratio between extra diagonal and diagonal coefficients : %.3e\"\n" %max_ratio_extra_diag )
        fiche.write("  description \"Matrix norm : %.3e \"\n" %matrix_norm )
        fiche.write("  description \"Analytical solution : "+analytical_solution+"\" \n")
        fiche.write("  description \""+str( solution_P1 )+"\"\n")
        fiche.write("  description \" facsec : "+facsec+"\" \n")
        fiche.write("  Figure {\n")
        fiche.write("    dimension 2 \n")
        fiche.write("    include_description_curves 0 \n")
        fiche.write("    labely Error \n")
        fiche.write("    labelx time\n")
        fiche.write("    Curve {\n")
        fiche.write("      file "+str( case_number )+"/mydata_ERREUR_ABSOLUE.son \n")
        fiche.write("      legend Absolute \n")
        fiche.write("      style lines \n")
        fiche.write("    }\n")
        fiche.write("    Curve {\n")
        fiche.write("      file "+str( case_number )+"/mydata_ERREUR_RELATIVE.son \n")
        fiche.write("      legend Relative \n")
        fiche.write("      style lines \n")
        fiche.write("    }\n")
        fiche.write("  }\n")
        fiche.write("}\n")


def PlotErrorInFunctionOfParameter( parameters_name, error_file ) :
    nb_parameters = len( parameters_name )
    the_title = "Error in function of "
    for p in range( 0, nb_parameters ) :
        parameter_name = parameters_name[ p ]
        parameter_name = parameter_name[1:-1] #suppression des '@'
        the_title = the_title + parameter_name

    already_written = 0
    with open( "../fiche.prm", "r" ) as fiche :
        for line in fiche :
            if line == the_title :
                already_written=1

    if already_written == 0 :
        with open( "../fiche.prm", "a" ) as fiche :
            fiche.write("chapter \n")
            fiche.write("{\n")
            fiche.write("  title \""+the_title+"\" \n")
            fiche.write("  Figure {\n")
            fiche.write("    dimension 2 \n")
            fiche.write("    include_description_curves 0 \n")
            fiche.write("    labely Error \n")
            fiche.write("    labelx "+parameters_name[ 0 ]+"\n")
            fiche.write("    Curve {\n")
            fiche.write("      file "+str( error_file )+" \n")
            fiche.write("      columns $1 $2 \n")
            fiche.write("      legend Absolute \n")
            fiche.write("      style lines \n")
            fiche.write("    }\n")
            fiche.write("    Curve {\n")
            fiche.write("      file "+str( error_file )+" \n")
            fiche.write("      columns $1 $3 \n")
            fiche.write("      legend Relative \n")
            fiche.write("      style lines \n")
            fiche.write("    }\n")
            fiche.write("  }\n")
            fiche.write("}\n")

