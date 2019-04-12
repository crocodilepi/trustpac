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
  
  //define the sampling
  Int_t sampling_size = 5000;
  TSampling *sampling = new TSampling( tds, "lhs", sampling_size );
  sampling->generateSample();

  // The output file of the code is written in the file 'quantity_of_interest'
  TOutputFileRow *fout = new TOutputFileRow( "quantity_of_interest" ); 

  // The quantity of interest is the velocity norm 
  TAttribute * velocity_norm = new TAttribute( "velocity_norm" );
  velocity_norm->setDefaultValue( 0.0 );
  fout->addAttribute( velocity_norm, 1 ) ; // 1ere (et unique) colonne du fichier

  TCode *mycode = new TCode( tds, "../../code_TRUST.sh upwind" );
  mycode->addOutputFile( fout );


  TLauncher *tl = new TLauncher( tds, mycode );
  tl->run( );

  // save results in order to be able to draw with DrawCobweb.C
  tds->exportData( "data_save.dat" );
  
  //setup for quantiles computation
  
  nProba = 2; //nProba = number of quantiles we want

  proba = new Double_t[ nProba ]; //will contain percentage of quantile (for example 2.5% and 97.5%)
  proba[0] = 0.025; // 2.5%
  proba[1] = 0.975; // 97.5%
  quantiles = new Double_t[nProba];

  //quantiles will be written in a file
  ofstream quantile_file("quantiles.txt", ios::out | ios::trunc);

  //now compute quantiles
  tds->computeQuantile( "velocity_norm", nProba, proba, quantiles );
  quantile_file << " lower_bound : "<< quantiles[ 0 ] << " upper_bound" << quantiles[ 1 ] << endl;
  cout << "quantiles[  2.5% ] = "<<quantiles[ 0 ] << endl;
  cout << "quantiles[ 97.5% ] = "<<quantiles[ 1 ] << endl;
  cout << " ====== "<< endl;
  
  quantile_file.close( );

  // Draw the cobweb plot
  tds->getTuple()->Draw("seuil_statio:seuil_NS:Omega:velocity_norm","","para"); 
  TParallelCoord* para = ( TParallelCoord* )gPad->GetListOfPrimitives( )->FindObject( "ParaCoord" );
  TParallelCoordVar* axis = ( TParallelCoordVar* )para->GetVarList( )->FindObject( "seuil_NS" );
  axis->AddRange( new TParallelCoordRange( axis, 1e-8, 5e-5 ) ); 
  para->AddSelection("blue");
  para->GetCurrentSelection()->SetLineColor(kBlue);

}
