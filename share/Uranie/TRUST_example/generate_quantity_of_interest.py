# -*- coding: latin-1 -*-

# Le but de ce script est de calculer la norme de la vitesse en un point ( en 2D )

import sys
from math import sqrt


nb_args = len( sys.argv )
if nb_args != 2 :
    print "Error in generate_quantity_of_interest.py"
    print "You must specify the name of the file containing data to manipulate"
    sys.exit( 255 )

filename = sys.argv[ 1 ]

with open( filename, 'r') as input_file: 
    with open( "quantity_of_interest", 'wb') as output_file:
        line_idx = 0
        for line in input_file:
            if line_idx > 0 :
                print "Error while reading "+str( filename )+" because there are more than one line !"
                sys.exit( 255 )
            line_idx = line_idx + 1
            all_cols = line.split( ) # récupération de toutes les colonne 
            nb_cols = len( all_cols )
            if( nb_cols != 3 ) :
                print "Error while reading "+str( filename )+" because there are "+str( nb_cols )+" columns instead of 3 as expected"
                sys.exit( 255 )
            Vx = float( all_cols[ 1 ] )
            Vy = float( all_cols[ 2 ] )
            V = sqrt( Vx*Vx + Vy*Vy )
            output_file.write( "%.15e" %V ) # Ecriture de la norme de la vitesse avec une precision de 15
                

