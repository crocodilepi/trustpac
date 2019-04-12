import sys
import subprocess
import math
import solution
import AddChapter

def IsFloat( value ):
  try:
    float( value )
    return 1
  except ValueError:
    return 0

def sampling ( LowerBound, UpperBound, Method, size ) :
  print "ToolBox.sampling"
  sample = []
  if float( size ) < 2 :
    print "Error in ToolBox.sampling, the sample size must be greater than 1 ! "
    sys.exit( 255 )
  if Method == "equidistributed" :
    total =  float( UpperBound ) - float( LowerBound ) 
    step  = total / ( int( size ) - 1 )
    for idx in range( 0, int( size ) ) :
      value = LowerBound + idx*step
      sample.append( value )
  elif Method == "log" :
    total =  math.log10( UpperBound ) - math.log10( LowerBound ) 
    step = total / ( int( size ) - 1 )
    for idx in range( 0, int( size ) ) :
      exp = idx*step
      value = LowerBound + pow( 10, exp )
      sample.append( value )
  else :
    print "Error in ToolBox.sampling the method "+str( Method )+" is not recognized"

  return sample




# Si on utilise un echantillon predefini et qu'on a 2 parametres qui varient
# alors on a d'abords l'echantillonnage associe au premier parametre, puis ensuite celui associe au second
def CreateTestCases( ParameterNames, LowerBounds, UpperBounds, Methods, SampleSize, PredefinedSample ) :
  print "ToolBox.CreateTestCases"
  nb_parameters = len( ParameterNames )

  nb_elements = int( SampleSize )
  use_predefined_sample = 0
  if len( PredefinedSample ) > 0 :
    use_predefined_sample = 1
    nb_elements = len( PredefinedSample )

  # Creation de tous les repertoires
  cmd = "ini=$PWD;"
  cmd = cmd + "cd .. ;"
  cmd = cmd + "mkdir variability ;"
  cmd = cmd + "cd variability ;"
  for e in range(0, nb_elements ) :
    rep_name = "v"+str( e )
    cmd = cmd + "mkdir "+rep_name+" ;"
    cmd = cmd + "cd "+rep_name+" ;"
    cmd = cmd + "cp $ini/generic.data . ;"
    cmd = cmd + "cp $ini/parameters . ;"
    cmd = cmd + "cd ../ ;"
  print "commande "+str( cmd )
  subprocess.call( cmd , shell=True )
    
  # Echantillonnage
  the_sample = []
  if use_predefined_sample == 1 :
    the_sample = PredefinedSample
  else :
    for p in range(0, nb_parameters ) :
      s = sampling( LowerBounds[ p ], UpperBounds[ p ], Methods[ p ], SampleSize )
      the_sample = the_sample + s

  # Mise a jour du fichier parameters en fonction de l'echantillonnage    
  cmd = "ini=$PWD;"
  cmd = cmd + "cd .. ;"
  cmd = cmd + "cd variability ;"
  idx = 0
  for e in range( 0, nb_elements ) :
    rep_name = "v"+str( e )
    cmd = cmd + "cd "+rep_name+" ;"
    cmd = cmd + "sed \" \n"
    for p in range(0, nb_parameters ) :
      cmd=cmd+"s/"+ParameterNames[ p ]+"/"+str( the_sample[ e ] )+"/g \n" # l'ajout de g permet de traiter toutes les occurences d'une ligne
    cmd=cmd+"\""
    cmd=cmd+" generic.data > mydata.data ;"
  subprocess.call([cmd], shell=True)




def AddNewTestCasesToPrm(  ) :
  print "ToolBox.AddNewTestCasesToPrm"
  nb_args = len( sys.argv )
  with open( "fiche.prm", "a" ) as fiche :
    for idx in range( 2, nb_args ) :
      arg = sys.argv[ idx ]
      test_case = arg[2:] # l'argument est sous la forme ./case4 donc on enleve les 2 premiers caracteres
      fiche.write("  testcase "+test_case+" mydata.data \n")


def AddNewTestCaseToPrm( test_case ) :
  print "ToolBox.AddNewTestCaseToPrm"
  with open( "fiche.prm", "a" ) as fiche :
    fiche.write("  testcase "+test_case+" mydata.data \n")


def read_last_line_of_file( filename ) :
  nb_lines = 0
  line_idx=0
  error=-1e38
  # d'abord compter le nombre de ligne
  with open( filename ,"r" ) as f :
    for line in f :
      nb_lines = nb_lines + 1
  with open( filename ,"r" ) as f :
    # puis recuperer la derniere
    for line in f :
      line_idx = line_idx + 1
      if line_idx == nb_lines :
        all_values = line.split( )
        time = all_values[ 0 ]
        error = all_values[ 1 ]
  return error

# Cette fonction fait appel a des commandes shell pour creer le fichier error_file si besoin ainsi que le post_run
def error_in_function_of( variables_name, variables_value, error_file ) :  
  nb_variables = len( variables_name.split( ) )
  cmd = "./error_in_function_of.sh "+str( nb_variables )+" "+variables_name+" "+variables_value+" "+error_file
  subprocess.call( cmd , shell=True )

# Le post_run fait appel a cette fonction
def add_error_in_function_of( ) :
  nb_args=len( sys.argv )
  print "add_error_in_fonction_of arguments : "+str( sys.argv )
  nb_values = nb_args - 1
  # les valeurs commencent a sys.argv[ 2 ]
  error_file = sys.argv[ nb_values ]
  all_values=""
  for idx in range( 2, nb_values ) :
    all_values = all_values +sys.argv[ idx ]+"\t"
    print "all_values = "+str( all_values )
  filename = "../"+error_file
  absolute_error = read_last_line_of_file( "mydata_ERREUR_ABSOLUE.son" )
  relative_error = read_last_line_of_file( "mydata_ERREUR_RELATIVE.son" )
  print "add_error variables_value "+str( all_values )
  print "add_error filename "+str( filename )
  print "add_error absolute_error "+str( absolute_error )
  print "add_error relative_error "+str( relative_error )
  with open( filename, "a" ) as err :
    err.write( str(all_values)+"\t"+str(absolute_error)+"\t"+str(relative_error)+"\n" )

  

def GetKeywords( ) :
  keywords=["@puissance@",          # 0
            "@initial_conditions@", # 1
            "@nb_pas_dt_max@",      # 2 
            "@facsec@",             # 3
            "@K00@",                # 4
            "@K01@",                # 5
            "@K10@",                # 6
            "@K11@",                # 7
            "@A@",                  # 8
            "@B@",                  # 9
            "@C@",                  # 10
            "@erreur_periode@",     # 11
            "@Maillage@",           # 12
            "@K@",                  # 13
            "@utilisation_operateur_tensoriel@", # 14
            "@lambda@",             # 15
            "@D@",                  # 16
            "@E@"]                  # 17

  return keywords


def ReadParametersFile( ) :
  print "ToolBox.ReadParametersFile"
  parameters_values=[ ]

  nb_variables = 0
  variables_name=""
  variables_value=""

  chapter = "no_chapter"
  post_process_error = 0

  # initialisation
  keywords = GetKeywords( )
  size=len( keywords )
  values=[]
  for idx in range(0, size ):
    parameters_values.append("default_value" )

  # Lecture des parametres
  with open("parameters","r") as f :
    for line in f :
      all_words = line.split( )
      nb_words = len( all_words )
      if nb_words > 0 :
        first_word = all_words[ 0 ]
        if first_word in keywords :
          idx = keywords.index( first_word )
          if idx == 1  : # initial conditions 
            parameters_values[ idx ] = ""
            for value_idx in range( 1, nb_words ) :
              parameters_values[ idx ] = parameters_values[ idx ] +" "+ all_words[ value_idx ]
          else :
            parameters_values[ idx ] = all_words[ 1 ]
        else : # else associe au first_word
          if first_word == "chapter" or first_word == "Chapter" :
            chapter = all_words[ 1 ]
          if first_word == "error_in_function_of" :
            print "line : "+str( line )
            nb_variables = int( all_words[ 1 ] )
            idx = 2
            for p in range( 0, nb_variables ) :
              i = idx+p
              print "on va taper dans la valeur "+str( i )+" car nb_variables = "+str( nb_variables )+" et idx = "+str( idx )+" et p "+str( p )
              variables_name = variables_name + " "+all_words[ i ] 
            idx = 2 + nb_variables
            error_file = all_words[ idx ]
            post_process_error = 1
          if first_word == "variability" :
            print "The variability keyword is not yet available !"
            print "Aborting..."
            cmd = "exit 255"
            subprocess.call( cmd , shell=True )
            use_variability = 1
            nb_params = int( all_words[ 1 ] )
            parameter_names = [ ]
            lower_bounds = [ ]
            upper_bounds = [ ]
            methods = [ ]
            for p in range( 0, nb_params ) :
              idx= 2 + p * 3
              parameter_names.append( all_words[ idx + 0 ] )
              lower_bounds.append( float( all_words[ idx + 1 ] ) )
              upper_bounds.append( float( all_words[ idx + 2 ] ) )
              methods.append( "log" )
            nb_runs = all_words[ nb_words - 1 ]
            CreateTestCases( parameter_names, lower_bounds, upper_bounds, methods, nb_runs, [] )


  Tensoriel=parameters_values[ 14 ]
  K00=parameters_values[ 4 ]
  K01=parameters_values[ 5 ]
  K10=parameters_values[ 6 ]
  K11=parameters_values[ 7 ]
  K=parameters_values[ 13 ]
  if Tensoriel == "true" :
    parameters_values[ 15 ] = "champ_fonc_xyz dom 4 "+str( K00 )+" "+str( K01 )+" "+str( K10 )+" "+str( K11 )
  else :
    parameters_values[ 15 ] = "champ_uniforme 1 "+str( K )

  if post_process_error == 1 :
    split_variables_name = variables_name.split( )
    nb_names = len( split_variables_name )
    for idx in range( 0, nb_names ) :
      variable_name = split_variables_name[ idx ]
      if variable_name in keywords :
        idx = keywords.index( variable_name )
        variables_value = variables_value+" "+str( parameters_values[ idx ] )
      else :
        print "erreur, la variable "+variable_name+" n'est pas dans la liste des variables !"
        sys.exit( 255 )
    error_in_function_of( variables_name, variables_value, error_file )
    # AddChapter.PlotErrorInFunctionOfParameter( variables_name, error_file )

  return ( parameters_values, chapter )





def generate_2D_power( parameters_values ) :

  print "ToolBox.generate_2D_power"
  A=parameters_values[ 8 ]
  B=parameters_values[ 9 ]
  C=parameters_values[ 10 ]
  D=parameters_values[ 16 ]
  E=parameters_values[ 17 ]
  Tensoriel=parameters_values[ 14 ]
  K00=parameters_values[ 4 ]
  K01=parameters_values[ 5 ]
  K10=parameters_values[ 6 ]
  K11=parameters_values[ 7 ]
  K=parameters_values[ 13 ]
  solution.deux( A, B, C, D, E, K00, K01, K10, K11, K, Tensoriel )



    
def create_mydata_file( parameters_values, compute_power_term ) :

  print "ToolBox.create_mydata_file with compute_power_term "+str( compute_power_term )
  keywords = GetKeywords( )
  size = len( keywords )

  # Lecture du terme de puissance si ce dernier a ete calcule
  # sinon on le remplace par son nom
  if compute_power_term == 1 : 
    with open("puissance","r") as f : 
      for line in f :
        all_words = line.split( )
        if len( all_words ) != 1 :
          print "Erreur pendant la lecture du fichier de puissance"
          cmd="exit 255 "
          subprocess.call([cmd], shell=True)
        parameters_values[ 0 ] = all_words[ 0 ]
  else :
    parameters_values[ 0 ] = keywords[ 0 ]
        
  print "Lancement de la commande sed"
  # Lancement de la commande sed
  Tensoriel=parameters_values[ 14 ]
  cmd="sed \" \n"
  for idx in range( 0, size ):
    cmd=cmd+"s/"+keywords[ idx ]+"/"+parameters_values[ idx ]+"/g \n" # l'ajout de g permet de traiter toutes les occurences d'une ligne
  if Tensoriel != "true" :
    cmd=cmd+"s/diffusion { Tensor }/diffusion { }/ \n"
  cmd=cmd+"\""
  cmd=cmd+" generic.data > mydata.data"
  subprocess.call([cmd], shell=True)

 

def PossibilityToComputePowerTerm( parameters_values ) :
  print "ToolBox.PossibilityToComputePowerTerm"
  indices = [ 8, 9, 10, 16, 17, 4, 5, 6, 7, 13 ]
  for idx in indices :
    v = parameters_values[ idx ]
    isfloat = IsFloat( v )
    if isfloat == 0 :
      return 0
  return 1



if __name__ == "__main__" :
  cmd="exit 255"
  nb_args = len( sys.argv )
  if nb_args == 1 :
    print "il faut au moins un argument (le nom de la fonction) lorsqu'on fait appel a ToolBox !"
    subprocess.call( cmd , shell=True )
  else :
    if sys.argv[ 1 ] == "AddNewTestCasesToPrm" :
      AddNewTestCasesToPrm( )
    elif sys.argv[ 1 ] == "add_error_in_function_of" :
      add_error_in_function_of( )
    else :
      print "Erreur, le nom de la fonction "+sys.argv[ 1 ]+" n'est pas reconnu "
      subprocess.call( cmd , shell=True )
