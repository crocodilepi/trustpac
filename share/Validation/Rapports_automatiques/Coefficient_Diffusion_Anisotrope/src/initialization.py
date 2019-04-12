import sys
import ToolBox
import AddChapter

case_number = sys.argv[ 1 ]
print "case_number "+case_number

# Lecture du fichier parameters qui va donner les valeurs des parametres mais aussi si on souhaite ecrire un chapitre dans la prm
parameters_values, chapter = ToolBox.ReadParametersFile( )

# si on decide de faire varier un parametre alors parameters_values contiendra un string (le nom du parametre qui varie)
# il ne faut donc pas essayer d'evaluer le terme de puissance car on a pas encore de valeur a donner a ce parametre
compute_power_term = ToolBox.PossibilityToComputePowerTerm( parameters_values )

if compute_power_term == 1 :
    ToolBox.generate_2D_power( parameters_values )
    
ToolBox.create_mydata_file( parameters_values, compute_power_term )

print "initialization chapter "+chapter
if chapter == "ExtraCoeffs" : 
    A=parameters_values[ 8 ]
    B=parameters_values[ 9 ]
    C=parameters_values[ 10 ]
    D=parameters_values[ 16 ]
    E=parameters_values[ 17 ]
    K00=parameters_values[ 4 ]
    K01=parameters_values[ 5 ]
    K10=parameters_values[ 6 ]
    K11=parameters_values[ 7 ]
    facsec=parameters_values[ 3 ]
    AddChapter.ExtraCoeffs( A, B, C, D, E, K00, K01, K10, K11, case_number, facsec )



