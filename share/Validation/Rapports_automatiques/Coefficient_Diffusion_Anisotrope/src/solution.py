# Coefficient de diffusion anisotrope K ( Matrice 2x2 )
#
# K00  K01
# K10  K11
#
# Solution T de la forme  Ax*x + B*y*y + C*x*y
#
# - div( K grad T ) = - 2 A K00 - 2 B K11 - C ( K01 + K10 )
#

import sys
from sympy import *

def deux( A, B, C, D, E, K00, K01, K10, K11, K, Tensoriel ) : 

    print "le script solution 2D is called with "

    xs = Symbol("x")
    ys = Symbol("y")

    # Remplissage du coefficient de diffusion
    K00 = float( K00 )
    K01 = float( K01 )
    K10 = float( K10 )
    K11 = float( K11 )
    K = float( K )

    # Solution analytique 
    A = float( A )
    B = float( B )
    C = float( C )
    D = float( D )
    E = float( E )
    
    F = 1e60

    if Tensoriel == "true" :
        # div K grad T
        GradTx =  2*A*xs + C * ys
        GradTy =  2*B*ys + C * xs
        KgradTx = K00 * GradTx + K01 * GradTy
        KgradTy = K10 * GradTx + K11 * GradTy

        # dGradTx_x = - 2*A
        # dGradTx_y = - C
        # dGradTy_y = - 2*B
        # dGradTy_x = - C
        # div_K_grad_T = K00 * dGradTx_x + K01 * dGradTy_x + K10 * dGradTx_y + K11 * dGradTy_y  
 
        # Calcul du second membre 
        F = -2*A*K00 - 2*B*K11 - C*( K01 + K10 )
    else :
        # Div K Grad T = K Delta T avec Delta T = 2A + 2B
        F = - K * ( 2*A + 2*B )


    with open( "puissance" , "w" ) as f :
        f.write( str(F) )







def debug( ) :

    # # Solutions aux bords
    # x = 0.0
    # nx =-1
    # ny = 0
    # SolutionGauche = A*x*x + B*ys*ys + C*x*ys
    # FluxGauche = KgradTx*nx

    # x = 1.0
    # nx = 1
    # SolutionDroit = A*x*x + B*ys*ys + C*x*ys
    # FluxDroit = KgradTx*nx

    # y = 0.0
    # nx = 0
    # ny = -1
    # SolutionBas = A*xs*xs + B*y*y + C*xs*y
    # FluxBas = KgradTy*ny

    # y = 1.0
    # ny = 1
    # SolutionHaut = A*xs*xs + B*y*y + C*xs*y
    # FluxHaut = KgradTy*ny

    # print "Solution sur le bord gauche : "+str( SolutionGauche )
    # print "Solution sur le bord droit  : "+str( SolutionDroit )
    # print "Solution sur le bord haut   : "+str( SolutionHaut )
    # print "Solution sur le bord bas    : "+str( SolutionBas )
    # print "Terme de puissance          : "+str( F )

    # print ""
    # print "Flux gauche : "+str( FluxGauche )
    # print "Flux droit : "+str( FluxDroit )
    # print "Flux haut : "+str( FluxHaut )
    # print "Flux bas : "+str( FluxBas )
    # print ""

    ### Codage
    # print "\nCODAGE\n"
    nb_comp=2
    Six = 1.0
    Siy = 0.0
    Sjx = 0.5
    Sjy = 0.5

    Diffu_S = K00 * Sjx + K01 * Sjy
    # print "Diffu_S : "+str( Diffu_S )
    DSiSj = Diffu_S * Six
    Diffu_S = K10 * Sjx + K11 * Sjy
    # print "Diffu_S : "+str( Diffu_S )
    DSiSj = DSiSj + Diffu_S * Siy
    # print "DSiSj : "+str( DSiSj )


    erreur_DSiSj= abs( DSiSj - ( 9.00000000000000e+01 )  )
    # print "erreur "+str( erreur_DSiSj )
