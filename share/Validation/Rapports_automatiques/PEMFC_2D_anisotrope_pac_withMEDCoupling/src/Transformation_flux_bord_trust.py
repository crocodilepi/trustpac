# -*- coding: latin-1 -*-

TRUST_file = "PEMFC_2D_anisotrope_pac_withMEDCoupling_pb_Diffusion_chaleur.face" 
output_file = "trust_flux_bord_withMEDCoupling.txt"

with open(TRUST_file, 'r') as infile :
    with open(output_file, 'w') as outfile :
        outfile.write("# x \t y \t surface \t flux/surface \t flux \n")
        lines = infile.readlines();
        # find the portion of last time of flux au bord
        line0 = 0;
        lineEnd = len(lines);
        t = 0;
        for i, line in enumerate(lines):
            if "Flux par face sur Haut" in line:
                line0 = i;
        #print line0
        for i, line in enumerate(lines[line0+1:]):
			if "Flux par face" in line:
				lineEnd = i; 
				break;
        #print lineEnd
        
        for line in lines[lineEnd:line0:-1]:
			linestr = line.split()
			x = linestr[4]
			y = linestr[6]
			s = linestr[8]
			fs = linestr[10]
			f = linestr[12]
			outfile.write(x+"\t"+y+"\t"+s+"\t"+fs+"\t"+f+"\n")
