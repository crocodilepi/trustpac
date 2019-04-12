Statistiques d'initialisation du calcul

Temps total                       0.517159

Statistiques de resolution du probleme

Temps total                       2.61585

Timesteps                         7
Secondes / pas de temps           0.373648
Dont solveurs Ax=B                0.0871894 23% (1 appel/pas de temps)
Dont assemblage matrice_implicite 0.00248614 0% (1 appel/pas de temps)
Dont operateurs diffusion         0.00128657 0% (1 appel/pas de temps)
Dont operateurs gradient          0.00509057 1% (21 appels/pas de temps)
Dont operations postraitement     0.2724 72% (1 appel/pas de temps)
Dont calcul dt                    5.67143e-05 0% (1 appel/pas de temps)
Dont operations sauvegarde        0.000614429 0% (1 appel/pas de temps)
Dont calcul divers                0.004523 1%
Nb echange_espace_virtuel / pas de temps 94
Nb solveur / pas de temps         1
Secondes / solveur                0.0871894
Iterations / solveur              1
I/O:
Nb sauvegardes         : 7
Secondes / sauvegarde  : 0.000614429
Donnees ecrites [Mo]   : 0.57
Debit           [Mo/s] : 132

Timesteps = nombre de pas de temps
Nb solveur = nombre de resolutions de systeme lineaire
Nb assemblage implicite = nombre d'assemblage de matrice pour le schema implicite
Iterations = nombre moyen d'iterations du solveur
Communications = fraction du temps passe dans les
                 communications entre processeurs (hors io fichiers)
Network latency = temps d'un mpsum mesure par un bench interne sur 0.1s
Network bandwidth = maximum sur l'ensemble des processeurs
                    de la bande passante moyenne des operations send_recv
Waiting time = estimation du temps d'attente des differents processeurs

Max_waiting_time grand  => probablement mauvais decoupage
Communications > 30%    => trop de processeurs ou reseau trop lent

Statistiques de post resolution

Temps total                       0.090899

