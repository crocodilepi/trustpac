Statistiques d'initialisation du calcul

Temps total                       1.66573

Statistiques de resolution du probleme

Temps total                       136.661

Timesteps                         807
Secondes / pas de temps           0.169307
Dont solveurs Ax=B                0.0855878 50% (1 appel/pas de temps)
Dont assemblage matrice_implicite 0.0309639 18% (1 appel/pas de temps)
Dont operateurs diffusion         0.0084013 4% (1 appel/pas de temps)
Dont operateurs gradient          8.77113e-05 0% (0.0210657 appels/pas de temps)
Dont operations postraitement     0.0114373 6% (1 appel/pas de temps)
Dont calcul dt                    0.000633971 0% (1 appel/pas de temps)
Dont operations sauvegarde        0.00139379 0% (0.972739 appels/pas de temps)
Dont calcul divers                0.0308011 18%
Nb echange_espace_virtuel / pas de temps 12.1685
Nb solveur / pas de temps         1
Secondes / solveur                0.0855878
Iterations / solveur              1
I/O:
Nb sauvegardes         : 785
Secondes / sauvegarde  : 0.00143285
Donnees ecrites [Mo]   : 63.80
Debit           [Mo/s] : 56

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

Temps total                       1.38002

