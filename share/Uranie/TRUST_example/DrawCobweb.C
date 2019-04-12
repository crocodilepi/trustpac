{
  // Define the DataServer
  TDataServer *tds = new TDataServer( "tds", "Data server" );
  tds->fileDataRead( "data_save.dat" );
  
  //setup for quantiles computation
  
  nProba = 2; //nProba = number of quantiles we want

  proba = new Double_t[ nProba ]; //will contain percentage of quantile (for example 2.5% and 97.5%)
  proba[ 0 ] = 0.025; // 2.5%
  proba[ 1 ] = 0.975; // 97.5%
  quantiles = new Double_t[ nProba ];

  //quantiles will be written in a file
  ofstream quantile_file( "quantiles.txt", ios::out | ios::trunc );

  //now compute quantiles
  tds->computeQuantile( "velocity_norm", nProba, proba, quantiles );
  quantile_file << " lower_bound : "<< quantiles[ 0 ] << " upper_bound" << quantiles[ 1 ] << endl;
  quantile_file.close( );

  // display some values
  cout << "quantiles[  2.5% ] = "<<quantiles[ 0 ] << endl;
  cout << "quantiles[ 97.5% ] = "<<quantiles[ 1 ] << endl;
  cout << " ====== "<< endl;
  
  // Draw the cobweb plot
  tds->getTuple( )->Draw( "seuil_statio:seuil_NS:Omega:velocity_norm", "", "para" ); 
  TParallelCoord* para = ( TParallelCoord* )gPad->GetListOfPrimitives( )->FindObject( "ParaCoord" );
  TParallelCoordVar* axis = ( TParallelCoordVar* )para->GetVarList( )->FindObject( "seuil_NS" );
  axis->AddRange( new TParallelCoordRange( axis, 1e-8, 5e-5 ) ); 
  para->AddSelection( "blue" );
  para->GetCurrentSelection( )->SetLineColor( kBlue );

  // historique des seuils utilises pour les graph de cobweb
  // axis->AddRange( new TParallelCoordRange( axis, 1.4, 1.5 ) ); // Omega  
  // axis->AddRange( new TParallelCoordRange( axis, 1.176, 1.177 ) ); // velocity_norm  
}
