

Avant tout chose, ne pas oublier de sourcer l'environement TRUST

Le lancement d'un script MyScript.C se fait via la commande : root -l MyScript.C
La sortie de root se fait via la commande : .q


Cas d'étude : Upwind ( cas test de TRUST )
=============

Il s'agit d'un cas 2D, qui tourne très vite, qui simule l'écoulement d'un fluide autour d'un obstacle.
Une condition de non glissement (paroi fixe) est appliquée sur les contours de l'obstacle, tandis que le champ de vitesse est imposé en entrée et la pression est imposée en sortie.
Une condition de symétrie est appliquée partout ailleurs.

Il a été décidé de faire porter l'incertitude sur les paramètres suivants :
- le "seuil_statio" utilisé dans le schéma en temps
- le seuil de convergence du solveur pression ( GCP ) utilisé pour résoudre Navier Stokes
- la valeur de "Omega" utilisé par le préconditionneur SSOR du solveur pression


Organisation des scripts
=======================

Le script rootlogon.C est un script ESSENTIEL au bon fonctionnement de Uranie 1.4.0 
Le script sobol.C permet de calculer les indices de Sobol
Le script uncertainties.C permet de calculer les quantiles et de générer un affichage de type Cobweb
Le script DrawCobweb.C permet de calculer les quantiles et de générer un affichage de type Cobweb à partir des résultats générés par le script uncertainties.C


==== Résultats ====

Nombre de runs : 5 000


quantiles[  2.5% ] = 1.17347
quantiles[ 97.5% ] = 1.17719
=> on est sûr à 95% que la norme du champ de vitesse sera comprise entre 1.17347 et 1.17719

 ** Input att [seuil_statio] First [0] Total Order[0.001]
 ** Input att [seuil_NS] First [0.622923] Total Order[0.804627]
 ** Input att [Omega] First [0.184073] Total Order[0.415222]
=> C'est bien seuil_NS qui est le plus influent 



