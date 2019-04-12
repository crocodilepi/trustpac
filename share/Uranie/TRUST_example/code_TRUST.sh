# -*- coding: latin-1 -*-

##
# Script permettant de lancer le calcul
# puis de preparer la sortie attendue par Uranie
# ici la quantité d'intérêt est la norme du champ de vitesse en un point
## 

if [ "$TRUST_ROOT" = "" ]
then
    echo "TRUST environment not set !!"
    exit 1
fi
   
# options  
nb_procs=1             # nombre de procs
cp ../../upwind.geo .  # il est necessaire de copier le .geo
cp ../../generate_quantity_of_interest.py . # script permettant de construire la norme de la vitesse en un point

if [ "$nb_procs" = 1 ]
then
    # Lancement du calcul en séquentiel
    trust upwind &> .code_out
    tail -1 upwind_SONDE_VITESSE.son > vitesse_last_line
    python generate_quantity_of_interest.py vitesse_last_line
else
    # Préparation du calcul parallel
    make_PAR.data $trust_case $nb_procs &> .make_par
    # Lancement du calcul
    trust PAR_$trust_case $nb_procs &> .code_out
    tail -1 PAR_upwind_SONDE_VITESSE.son > vitesse_last_line
    python generate_quantity_of_interest.py vitesse_last_line
fi
