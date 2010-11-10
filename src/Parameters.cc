template<typename T>
Parameters<T>::Parameters() :
    n_surfs(1),   
    sh_order(12),
    filter_freq(8),
    bending_modulus(1e-2),
    position_solver_iter(15),
    tension_solver_iter(15),
    position_solver_restart(1),
    tension_solver_restart(1),
    position_solver_tol((typeid(T) == typeid(float)) ? 1e-4: 1e-8),
    tension_solver_tol((typeid(T) == typeid(float)) ? 5e-4: 1e-8),
    n_steps(-1),
    time_horizon(1),
    ts(1),
    scheme(Explicit),
    singular_stokes(ViaSpHarm),
    rep_maxit(10),
    rep_up_freq(24),
    rep_filter_freq(4),
    rep_ts(1),
    rep_tol((typeid(T) == typeid(float)) ? 1e-3 : 1e-4),
    bg_flow_param(1e-1),
    save_data(false),
    save_stride(-1),
    error_factor(1)
{}

template<typename T>
Parameters<T>::Parameters(int argc, char** argv)
{
  this->parseInput(argc, argv);
}

template<typename T>
Parameters<T>::~Parameters()
{}

template<typename T>
Error_t Parameters<T>::parseInput(int argc, char** argv)
{
  // 1. CREATE AN OBJECT
  AnyOption opt;

  // 2. SET PREFERENCES
  // opt.noPOSIX(); // do not check for POSIX style character options
  opt.setVerbose(); // print warnings about unknown options
  opt.autoUsagePrint(true); // print usage for bad options 

  // 3. SET THE USAGE/HELP   
  this->setUsage(&opt);
     
  // 4. SET THE OPTIONS
  this->setOptions(&opt);
 
  // 5. PROCESS THE COMMANDLINE AND RESOURCE FILE, COMMAND LINE
  // OPTIONS OVERRIDE THE FILE
  opt.processCommandArgs( argc, argv );

  if( opt.getValue( 'f' ) != NULL  || opt.getValue( "optFile" ) )
    opt.processFile( opt.getValue('f') );  
  
  opt.processCommandArgs( argc, argv );

  //print usage if no options
  if( ! opt.hasOptions() ) { 
    opt.printUsage();
    return InvalidParameter;
  }

  // 6. GET THE VALUES 
  this->getOptionValues(&opt);

  return Success;
}

template<typename T>
void Parameters<T>::adjustFreqs()
{
  this->filter_freq = 2 * this->sh_order / 3;
  this->rep_up_freq = 2 * this->sh_order;
  this->rep_filter_freq = this->sh_order / 3;
}

template<typename T>
void Parameters<T>::setUsage(AnyOption *opt)
{
  opt->addUsage( "" );
  opt->addUsage( " Usage: " );
  opt->addUsage( "" );
  opt->addUsage( "  -f  --optFile              The name of the options file to be parsed.");
  opt->addUsage( "  -h  --help                 Prints this help " );
  opt->addUsage( "  -i  --initFile             The file containing the initial shape of vesicles");
  opt->addUsage( "  -o  --outFile              The output file");
  opt->addUsage( "  -s  --saveData             Flag to save data to file" );
  opt->addUsage( "" );
  opt->addUsage( "      --bendingModulus       The bending modulus of the interfaces" );
  opt->addUsage( "      --bgFlowParam          Single parameter passed to the background flow class" );
  opt->addUsage( "      --centFile             The file containing the initial center points *NOT IMPEMENTED*");
  opt->addUsage( "      --filterFreq           The differentiation filter frequency" );
  opt->addUsage( "      --nSurfs               The number of surfaces" );
  opt->addUsage( "      --positionIterMax      Maximum number of iterations for the position solver" );
  opt->addUsage( "      --positionRestart      Maximum number of restarts for the position solver" );
  opt->addUsage( "      --positionTol          The tolerence for the position solver" );
  opt->addUsage( "      --repFilterFreq        The filter freq for the reparametrization" );
  opt->addUsage( "      --repMaxIter           Maximum number of reparametrization steps" );
  opt->addUsage( "      --repTimeStep          The Time step for the reparametrization" );
  opt->addUsage( "      --repTol               The absolute value tol on the velocity of reparametrization" );
  opt->addUsage( "      --repUpFreq            The upsampling frequency for the reparametrization" );
  opt->addUsage( "      --saveStride           The frequency of saving to file (in time scale)" );
  opt->addUsage( "      --scheme               The time stepping scheme {Explicit | SemiImplicit}" );
  opt->addUsage( "      --shOrder              The spherical harmonics order" );
  opt->addUsage( "      --singularStokes       The scheme for the singular stokes evaluation" );
  opt->addUsage( "      --tensionIterMax       Maximum number of iterations for the tension solver" );
  opt->addUsage( "      --tensionRestart       Maximum number of restarts for the tension solver" );
  opt->addUsage( "      --tensionTol           The tolerence for the tension solver" );
  opt->addUsage( "      --timeHorizon          The time horizon of the simulation" );
  opt->addUsage( "      --timeStep             The time step size" );
}

template<typename T>
void Parameters<T>::setOptions(AnyOption *opt)
{
  // by default all options will be checked on the command line and
  // from option/resource file

  // a flag (takes no argument), supporting long and short forms
  opt->setFlag( "help", 'h' );   
  opt->setFlag( "saveData", 's' );   

  //an option (takes an argument), supporting long and short forms
  opt->setOption( "optFile", 'f');
  opt->setOption( "initFile", 'i');
  opt->setOption( "outFile", 'o');
  
  //an option (takes an argument), supporting only long form
  opt->setOption( "bendingModulus" );
  opt->setOption( "bgFlowParam" );
  opt->setOption( "centFile");
  opt->setOption( "filterFreq" );
  opt->setOption( "nSurfs" );
  opt->setOption( "positionIterMax" );
  opt->setOption( "positionRestart" );
  opt->setOption( "positionTol" );
  opt->setOption( "repFilterFreq" );
  opt->setOption( "repMaxIter" );
  opt->setOption( "repTimeStep" );
  opt->setOption( "repTol" );
  opt->setOption( "repUpFreq" );
  opt->setOption( "saveStride" );
  opt->setOption( "scheme" );
  opt->setOption( "shOrder" );
  opt->setOption( "singularStokes" );
  opt->setOption( "tensionIterMax" );
  opt->setOption( "tensionRestart" );
  opt->setOption( "tensionTol" );
  opt->setOption( "timeHorizon" );
  opt->setOption( "timeStep" );

  //for options that will be checked only on the command and line not
  //in option/resource file 
  //opt.setCommandFlag( ... ); 

  // for options that will be checked only from the option/resource
  // file
  //opt.setFileOption(  ... ); 
}

template<typename T>
void Parameters<T>::getOptionValues(AnyOption *opt)
{
  if( opt->getFlag( "help" ) || opt->getFlag( 'h' ) ) 
    opt->printUsage();
  
  if( opt->getFlag( "saveData" ) || opt->getFlag( 's' ) )
    this->save_data = true;

  //an option (takes an argument), supporting long and short forms
  if( opt->getValue( "initFile" ) != NULL || opt->getValue( 'i' ) !=NULL )
    this->init_file_name = opt->getValue( 'i' );
  
  if( opt->getValue( 'o' ) !=NULL || opt->getValue( "outFile") !=NULL )
    this->save_file_name = opt->getValue( 'o' );
  
  
  //shOrder set first so that the adjusted frequencies can be overridden
  if( opt->getValue( "shOrder" ) != NULL  )
  {
    this->sh_order =  atoi(opt->getValue( "shOrder" ));
    this->adjustFreqs();
  }

  if( opt->getValue( "bendingModulus" ) != NULL  )
    this->bending_modulus=  atof(opt->getValue( "bendingModulus" ));

  if( opt->getValue( "bgFlowParam" ) != NULL  )
    this->bg_flow_param =  atof(opt->getValue( "bgFlowParam" ));

  if( opt->getValue( "centFile") !=NULL )
    this->cntrs_file_name = opt->getValue( "centFile" );

  if( opt->getValue( "filterFreq" ) != NULL  )
    this->filter_freq =  atoi(opt->getValue( "filterFreq" ));

  if( opt->getValue( "nSurfs" ) != NULL  )
    this->n_surfs =  atoi(opt->getValue( "nSurfs" ));

  if( opt->getValue( "positionIterMax" ) != NULL  )
    this->position_solver_iter =  atoi(opt->getValue( "positionIterMax" ));
  
  if( opt->getValue( "positionRestart" ) != NULL  )
    this->position_solver_restart =  atoi(opt->getValue( "positionRestart" ));

  if( opt->getValue( "positionTol" ) != NULL  )
    this->position_solver_tol =  atof(opt->getValue( "positionTol" ));

  if( opt->getValue( "repFilterFreq"  ) != NULL  )
    this->rep_filter_freq =  atoi(opt->getValue( "repFilterFreq"  ));

  if( opt->getValue( "repMaxIter" ) != NULL  )
    this->rep_maxit =  atoi(opt->getValue( "repMaxIter" ));

  if( opt->getValue( "repTimeStep" ) != NULL  )
    this->rep_ts =  atof(opt->getValue( "repTimeStep" ));

  if( opt->getValue( "repTol" ) != NULL  )
    this->rep_tol =  atof(opt->getValue( "repTol" ));

  if( opt->getValue( "repUpFreq" ) != NULL  )
    this->rep_up_freq =  atoi(opt->getValue( "repUpFreq" ));

  if( opt->getValue( "saveStride" ) != NULL  )
    this->save_stride =  atof(opt->getValue( "saveStride" ));

  if( opt->getValue( "scheme" ) != NULL  )
    this->scheme = EnumifyScheme(opt->getValue( "scheme" ));

  if( opt->getValue( "singularStokes" ) != NULL  )
    this->singular_stokes = EnumifyStokesRot(opt->getValue( "singularStokes" ));

  if( opt->getValue( "tensionIterMax" ) != NULL  )
    this->tension_solver_iter =  atoi(opt->getValue( "tensionIterMax" ));

  if( opt->getValue( "tensionRestart" ) != NULL  )
    this->tension_solver_restart =  atoi(opt->getValue( "tensionRestart" ));
  
  if( opt->getValue( "tensionTol" ) != NULL  )
    this->tension_solver_tol =  atof(opt->getValue( "tensionTol" ));

  if( opt->getValue( "timeHorizon" ) != NULL  )
    this->time_horizon =  atof(opt->getValue( "timeHorizon" ));

  if( opt->getValue( "timeStep" ) != NULL  )
    this->ts =  atof(opt->getValue( "timeStep" ));

//   other methods: (bool) opt.getFlag( ... long or short ... )
}

template<typename T>
std::ostream& operator<<(std::ostream& output, const Parameters<T>& par)
{
    output<<"\n ===================================="<<std::endl;
    output<<"  Simulator parameters"<<std::endl;
    output<<" ===================================="<<std::endl;
    output<<"  Surface:"<<std::endl;
    output<<"    Number of surfaces       : "<<par.n_surfs<<std::endl;
    output<<"    SH order                 : "<<par.sh_order<<std::endl;
    output<<"    Filter freq              : "<<par.filter_freq<<std::endl;
    output<<"    Bending modulus          : "<<par.bending_modulus<<std::endl;
    output<<" ------------------------------------"<<std::endl;
    output<<"  Solver:"<<std::endl;
    output<<"    Position solver iter     : "<<par.position_solver_iter<<std::endl;
    output<<"    tension solver iter      : "<<par.tension_solver_iter<<std::endl;
    output<<"    Position solver restart  : "<<par.position_solver_restart<<std::endl;
    output<<"    Tension solver restart   : "<<par.tension_solver_restart<<std::endl;
    output<<"    Position solver tol      : "<<par.position_solver_tol<<std::endl;
    output<<"    Tension solver tol       : "<<par.tension_solver_tol<<std::endl;

    output<<" ------------------------------------"<<std::endl;
    output<<"  Time stepper:"<<std::endl;
    output<<"    Number of time steps     : "<<par.n_steps<<std::endl;
    output<<"    Time horizon             : "<<par.time_horizon<<std::endl;
    output<<"    Step size                : "<<par.ts<<std::endl;
    output<<"    Scheme                   : "<<par.scheme<<std::endl;
    output<<"    Singular Stokes          : "<<par.singular_stokes<<std::endl;
    output<<" ------------------------------------"<<std::endl;
    output<<"  Reparametrization:"<<std::endl;
    output<<"    Rep maxit                : "<<par.rep_maxit<<std::endl;
    output<<"    Rep upsample freq        : "<<par.rep_up_freq<<std::endl;
    output<<"    Rep filter freq          : "<<par.rep_filter_freq<<std::endl;
    output<<"    Rep step size            : "<<par.rep_ts<<std::endl;
    output<<"    Rep tol                  : "<<par.rep_tol<<std::endl;
    output<<" ------------------------------------"<<std::endl;
    output<<"  Initialization:"<<std::endl;
    output<<"    Init file name           : "<<par.init_file_name<<std::endl;
    output<<"    Centers file name        : "<<par.cntrs_file_name<<std::endl;
    output<<" ------------------------------------"<<std::endl;
    output<<"  Saving:"<<std::endl;
    output<<"    Save data                : "<<std::boolalpha<<par.save_data<<std::endl;
    output<<"    Save file name           : "<<par.save_file_name<<std::endl;
    output<<"    Save stride              : "<<par.save_stride<<std::endl;
    output<<"    Error Factor             : "<<par.error_factor<<std::endl;
    output<<" ------------------------------------"<<std::endl;
    output<<"  Background flow:"<<std::endl;
    output<<"    Background flow parameter: "<<par.bg_flow_param<<std::endl;
    output<<" ===================================="<<std::endl<<std::endl;

    return output;
}
