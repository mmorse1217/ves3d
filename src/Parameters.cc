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
    upsample_interaction(false),
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

  // 5. PROCESS THE COMMANDLINE AND RESOURCE FILE
  opt.processCommandArgs( argc, argv );
  if( opt.getValue( 'h' ) != NULL  || opt.getValue( "help" ) ){
      opt.printUsage();
      exit(0);
  }

  if( opt.getValue( 'f' ) != NULL  || opt.getValue( "opt-file" ) ){
      INFO("Parsing input file "<<opt.getValue('f'));

      if (!opt.processFile(opt.getValue('f')))
          return ErrorEvent::InvalidParameter;
  }
  // reporcess time commandline to override the file content
  opt.processCommandArgs( argc, argv );

  //print usage if no options
  if( ! opt.hasOptions() ) {
      CERR_LOC("You need to pass the simulation options"
          " either through command-line or file.",
          "",
          opt.printUsage());
      opt.printUsage();

      return ErrorEvent::InvalidParameter;
  }

  // 6. GET THE VALUES
  this->getOptionValues(&opt);
  return ErrorEvent::Success;
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
  opt->addUsage( "List of options: " );
  opt->addUsage( "" );
  opt->addUsage( " -f  --opt-file              The name of the options file to be parsed.");
  opt->addUsage( " -h  --help                  Prints this help " );
  opt->addUsage( " -i  --init-file             The file containing the initial shape of vesicles");
  opt->addUsage( " -o  --out-file              The output file");
  opt->addUsage( " -s  --save-data             Flag to save data to file" );
  opt->addUsage( "" );
  opt->addUsage( "     --bending-modulus       The bending modulus of the interfaces" );
  opt->addUsage( "     --bg-flow-param         Single parameter passed to the background flow class" );
  opt->addUsage( "     --cent-file             The file containing the initial center points");
  opt->addUsage( "     --error-factor          The permissible increase factor in area and volume error");
  opt->addUsage( "     --filter-freq           The differentiation filter frequency" );
  opt->addUsage( "     --n-surfs               The number of surfaces" );
  opt->addUsage( "     --position-iter-max     Maximum number of iterations for the position solver" );
  opt->addUsage( "     --position-restart      Maximum number of restarts for the position solver" );
  opt->addUsage( "     --position-tol          The tolerence for the position solver" );
  opt->addUsage( "     --rep-filter-freq       The filter freq for the reparametrization" );
  opt->addUsage( "     --rep-max-iter          Maximum number of reparametrization steps" );
  opt->addUsage( "     --rep-timestep          The Time step for the reparametrization" );
  opt->addUsage( "     --rep-tol               The absolute value tol on the velocity of reparametrization" );
  opt->addUsage( "     --rep-up-freq           The upsampling frequency for the reparametrization" );
  opt->addUsage( "     --save-stride           The frequency of saving to file (in time scale)" );
  opt->addUsage( "     --time-scheme           The time stepping scheme {Explicit | BlockImplicit}" );
  opt->addUsage( "     --sh-order              The spherical harmonics order" );
  opt->addUsage( "     --singular-stokes       The scheme for the singular stokes evaluation" );
  opt->addUsage( "     --tension-iter-max      Maximum number of iterations for the tension solver" );
  opt->addUsage( "     --tension-restart       Maximum number of restarts for the tension solver" );
  opt->addUsage( "     --tension-tol           The tolerence for the tension solver" );
  opt->addUsage( "     --time-horizon          The time horizon of the simulation" );
  opt->addUsage( "     --timestep              The time step size" );
  opt->addUsage( "     --upsample-interaction  Flag to whether upsample (and filter) the interaction force" );
}

template<typename T>
void Parameters<T>::setOptions(AnyOption *opt)
{
  // by default all options will be checked on the command line and
  // from option/resource file

  // a flag (takes no argument), supporting long and short forms
  opt->setCommandFlag( "help", 'h');
  opt->setFlag( "save-data", 's' );
  opt->setFlag( "upsample-interaction");


  //an option (takes an argument), supporting long and short forms
  opt->setOption( "init-file", 'i');
  opt->setOption( "out-file", 'o');

  //an option (takes an argument), supporting only long form
  opt->setOption( "bending-modulus" );
  opt->setOption( "bg-flow-param" );
  opt->setOption( "cent-file");
  opt->setOption( "error-factor" );
  opt->setOption( "filter-freq" );
  opt->setOption( "n-surfs" );
  opt->setOption( "position-iter-max" );
  opt->setOption( "position-restart" );
  opt->setOption( "position-tol" );
  opt->setOption( "rep-filter-freq" );
  opt->setOption( "rep-max-iter" );
  opt->setOption( "rep-timestep" );
  opt->setOption( "rep-tol" );
  opt->setOption( "rep-up-freq" );
  opt->setOption( "save-stride" );
  opt->setOption( "time-scheme" );
  opt->setOption( "sh-order" );
  opt->setOption( "singular-stokes" );
  opt->setOption( "tension-iter-max" );
  opt->setOption( "tension-restart" );
  opt->setOption( "tension-tol" );
  opt->setOption( "time-horizon" );
  opt->setOption( "timestep" );

  //for options that will be checked only on the command and line not
  //in option/resource file
  opt->setCommandOption( "opt-file", 'f');

  // for options that will be checked only from the option/resource
  // file
  //opt.setFileOption(  ... );
}

template<typename T>
void Parameters<T>::getOptionValues(AnyOption *opt)
{
  if( opt->getFlag( "save-data" ) || opt->getFlag( 's' ) )
      this->save_data = true;
  else
      this->save_data = false;

  if( opt->getFlag( "upsample-interaction" ) )
      this->upsample_interaction = true;
  else
      this->upsample_interaction = false;

  //an option (takes an argument), supporting long and short forms
  if( opt->getValue( "init-file" ) != NULL || opt->getValue( 'i' ) !=NULL )
    this->init_file_name = opt->getValue( 'i' );

  if( opt->getValue( 'o' ) !=NULL || opt->getValue( "out-file") !=NULL )
    this->save_file_name = opt->getValue( 'o' );


  //shOrder set first so that the adjusted frequencies can be overridden
  if( opt->getValue( "sh-order" ) != NULL  )
  {
    this->sh_order =  atoi(opt->getValue( "sh-order" ));
    this->adjustFreqs();
  }

  if( opt->getValue( "bending-modulus" ) != NULL  )
    this->bending_modulus=  atof(opt->getValue( "bending-modulus" ));

  if( opt->getValue( "bg-flow-param" ) != NULL  )
    this->bg_flow_param =  atof(opt->getValue( "bg-flow-param" ));

  if( opt->getValue( "cent-file") !=NULL )
    this->cntrs_file_name = opt->getValue( "cent-file" );

  if( opt->getValue( "filter-freq" ) != NULL  )
    this->filter_freq =  atoi(opt->getValue( "filter-freq" ));

  if( opt->getValue( "error-factor" ) != NULL  )
    this->error_factor =  atof(opt->getValue( "error-factor" ));

  if( opt->getValue( "n-surfs" ) != NULL  )
    this->n_surfs =  atoi(opt->getValue( "n-surfs" ));

  if( opt->getValue( "position-iter-max" ) != NULL  )
    this->position_solver_iter =  atoi(opt->getValue( "position-iter-max" ));

  if( opt->getValue( "position-restart" ) != NULL  )
    this->position_solver_restart =  atoi(opt->getValue( "position-restart" ));

  if( opt->getValue( "position-tol" ) != NULL  )
    this->position_solver_tol =  atof(opt->getValue( "position-tol" ));

  if( opt->getValue( "rep-filter-freq"  ) != NULL  )
    this->rep_filter_freq =  atoi(opt->getValue( "rep-filter-freq"  ));

  if( opt->getValue( "rep-max-iter" ) != NULL  )
    this->rep_maxit =  atoi(opt->getValue( "rep-max-iter" ));

  if( opt->getValue( "rep-timestep" ) != NULL  )
    this->rep_ts =  atof(opt->getValue( "rep-timestep" ));

  if( opt->getValue( "rep-tol" ) != NULL  )
    this->rep_tol =  atof(opt->getValue( "rep-tol" ));

  if( opt->getValue( "rep-up-freq" ) != NULL  )
    this->rep_up_freq =  atoi(opt->getValue( "rep-up-freq" ));

  if( opt->getValue( "save-stride" ) != NULL  )
    this->save_stride =  atof(opt->getValue( "save-stride" ));

  if( opt->getValue( "time-scheme" ) != NULL  )
    this->scheme = EnumifyScheme(opt->getValue( "time-scheme" ));

  if( opt->getValue( "singular-stokes" ) != NULL  )
    this->singular_stokes = EnumifyStokesRot(opt->getValue( "singular-stokes" ));

  if( opt->getValue( "tension-iter-max" ) != NULL  )
    this->tension_solver_iter =  atoi(opt->getValue( "tension-iter-max" ));

  if( opt->getValue( "tension-restart" ) != NULL  )
    this->tension_solver_restart =  atoi(opt->getValue( "tension-restart" ));

  if( opt->getValue( "tension-tol" ) != NULL  )
    this->tension_solver_tol =  atof(opt->getValue( "tension-tol" ));

  if( opt->getValue( "time-horizon" ) != NULL  )
    this->time_horizon =  atof(opt->getValue( "time-horizon" ));

  if( opt->getValue( "timestep" ) != NULL  )
    this->ts =  atof(opt->getValue( "timestep" ));

//   other methods: (bool) opt.getFlag( ... long or short ... )
}

template<typename T>
std::ostream& operator<<(std::ostream& output, const Parameters<T>& par)
{
    output<<"===================================="<<std::endl;
    output<<" Simulator parameters"<<std::endl;
    output<<"===================================="<<std::endl;
    output<<" Surface:"<<std::endl;
    output<<"   Number of surfaces       : "<<par.n_surfs<<std::endl;
    output<<"   SH order                 : "<<par.sh_order<<std::endl;
    output<<"   Filter freq              : "<<par.filter_freq<<std::endl;
    output<<"   Bending modulus          : "<<par.bending_modulus<<std::endl;
    output<<"------------------------------------"<<std::endl;
    output<<" Solver:"<<std::endl;
    output<<"   Position solver iter     : "<<par.position_solver_iter<<std::endl;
    output<<"   tension solver iter      : "<<par.tension_solver_iter<<std::endl;
    output<<"   Position solver restart  : "<<par.position_solver_restart<<std::endl;
    output<<"   Tension solver restart   : "<<par.tension_solver_restart<<std::endl;
    output<<"   Position solver tol      : "<<par.position_solver_tol<<std::endl;
    output<<"   Tension solver tol       : "<<par.tension_solver_tol<<std::endl;

    output<<"------------------------------------"<<std::endl;
    output<<" Time stepper:"<<std::endl;
    output<<"   Time horizon             : "<<par.time_horizon<<std::endl;
    output<<"   Step size                : "<<par.ts<<std::endl;
    output<<"   Scheme                   : "<<par.scheme<<std::endl;
    output<<"   Singular Stokes          : "<<par.singular_stokes<<std::endl;
    output<<"------------------------------------"<<std::endl;
    output<<" Reparametrization:"<<std::endl;
    output<<"   Rep maxit                : "<<par.rep_maxit<<std::endl;
    output<<"   Rep upsample freq        : "<<par.rep_up_freq<<std::endl;
    output<<"   Rep filter freq          : "<<par.rep_filter_freq<<std::endl;
    output<<"   Rep step size            : "<<par.rep_ts<<std::endl;
    output<<"   Rep tol                  : "<<par.rep_tol<<std::endl;
    output<<"------------------------------------"<<std::endl;
    output<<" Initialization:"<<std::endl;
    output<<"   Init file name           : "<<par.init_file_name<<std::endl;
    output<<"   Centers file name        : "<<par.cntrs_file_name<<std::endl;
    output<<"------------------------------------"<<std::endl;
    output<<" Saving:"<<std::endl;
    output<<"   Save data                : "<<std::boolalpha<<par.save_data<<std::endl;
    output<<"   Save file name           : "<<par.save_file_name<<std::endl;
    output<<"   Save stride              : "<<par.save_stride<<std::endl;
    output<<"   Error Factor             : "<<par.error_factor<<std::endl;
    output<<"------------------------------------"<<std::endl;
    output<<" Background flow:"<<std::endl;
    output<<"   Background flow parameter: "<<par.bg_flow_param<<std::endl;
    output<<"   Upsample Interaction     : "<<std::boolalpha<<par.upsample_interaction<<std::endl;
    output<<"====================================";

    return output;
}
