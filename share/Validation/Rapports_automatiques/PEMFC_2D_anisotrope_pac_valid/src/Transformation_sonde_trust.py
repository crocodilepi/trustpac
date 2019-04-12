# -*- coding: latin-1 -*-

TRUST_file = "PEMFC_2D_anisotrope_pac_SONDE_TEMPERATURE.son" 
output_file = "trust_sonde.txt"


x1 = 0.0
x2 = 10.0

# d'abord on va compter le nombre de valeurs a notre disposition
# et ensuite on va faire la conversion
with open(TRUST_file, 'r') as infile :
    with open(output_file, 'w') as outfile :
        outfile.write("# Position \t Temperature \n")
        last_line = infile.readlines( )[ -1 ]
        nb_values = len( last_line.split( ) )
        if nb_values == 0 :
            print "Error during the conversion of the TRUST probe file"
        else :
            delta = (x2 - x1) / ( nb_values - 1 )            
            all_values = last_line.split( )
            for i in range(1, nb_values) : # time is ignored
                probe = all_values[ i ]
                position = x1 + (i-1)*delta
                str_pos = "%.3e"%position
                outfile.write(str_pos+"\t"+probe+"\n")
