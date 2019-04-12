{
  // Define the DataServer
  TDataServer *tds = new TDataServer("tds", "Data server");

  // Add the study attributes 
  tds->addAttribute( new TUniformDistribution("seuil_statio" , 1e-8  , 1e-3 ) );
  tds->addAttribute( new TUniformDistribution("seuil_NS"     , 1e-8  , 1e-3 ) );
  tds->addAttribute( new TUniformDistribution("Omega"        , 1     , 2    ) );

  TString sJDDReferenceFlag = TString( "upwind.data" );
  tds->getAttribute( "seuil_statio" )->setFileFlag( sJDDReferenceFlag, "@seuil_statio@" );
  tds->getAttribute( "seuil_NS" )->setFileFlag( sJDDReferenceFlag, "@seuil_NS@" );
  tds->getAttribute( "Omega" )->setFileFlag( sJDDReferenceFlag, "@Omega@" );

  // The output file of the code is written in the file 'quantity_of_interest'
  TOutputFileRow *fout = new TOutputFileRow( "quantity_of_interest" ); 

  // The quantity of interest is the velocity norm 
  TAttribute * velocity_norm = new TAttribute( "velocity_norm" );
  velocity_norm->setDefaultValue( 0.0 );
  fout->addAttribute( velocity_norm, 1 ) ; // 1ere (et unique) colonne du fichier

  TCode *mycode = new TCode( tds, "../../code_TRUST.sh upwind" );
  mycode->addOutputFile( fout );

  Int_t ns = 5000;
  TSobol * tsobol = new TSobol( tds, mycode, ns );
  tsobol->generateSample( );
  tsobol->computeIndexes( );
}
