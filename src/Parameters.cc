template<typename T>
Parameters<T>::Parameters(){
    init();
}

template<typename T>
Parameters<T>::Parameters(int argc, char** argv)
{
    init();
    parseInput(argc, argv);
}

template<typename T>
Parameters<T>::Parameters(std::istream &is, Format format){
    unpack(is, format);
}

template<typename T>
void Parameters<T>::init()
{
    // ordered alphabetically
    bending_modulus	    = 1e-2;
    bg_flow                 = ShearFlow;
    bg_flow_param	    = 1e-1;
    error_factor	    = 1;
    filter_freq		    = 8;
    n_surfs		    = 1;
    num_threads             = 4;
    position_solver_iter    = 15;
    position_solver_restart = 1;
    position_solver_tol	    = (typeid(T) == typeid(float)) ? 1e-4: 1e-8;
    rep_filter_freq	    = 4;
    rep_maxit		    = 10;
    rep_tol		    = (typeid(T) == typeid(float)) ? 1e-3 : 1e-4;
    rep_ts		    = 1;
    rep_up_freq		    = 24;
    save_data		    = false;
    save_stride		    = -1;
    scheme		    = JacobiBlockImplicit;
    sh_order		    = 12;
    singular_stokes	    = ViaSpHarm;
    tension_solver_iter	    = 15;
    tension_solver_restart  = 1;
    tension_solver_tol	    = (typeid(T) == typeid(float)) ? 5e-4: 1e-8;
    time_horizon	    = 1;
    time_iter_max           = 100;
    time_precond            = NoPrecond;
    time_tol                = 1e-6;
    ts			    = 1;
    upsample_interaction    = false;
    viscosity_contrast      = 1.0;
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
          return ErrorEvent::InvalidParameterError;
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

      return ErrorEvent::InvalidParameterError;
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

  // ordered alphabetically
  opt->addUsage( "     --bending-modulus       The bending modulus of the interfaces" );
  opt->addUsage( "     --bg-flow-param         Single parameter passed to the background flow class" );
  opt->addUsage( "     --bg-flow-type          Type of the background flow" );
  opt->addUsage( "     --cent-file             The file containing the initial center points");
  opt->addUsage( "     --error-factor          The permissible increase factor in area and volume error");
  opt->addUsage( "     --filter-freq           The differentiation filter frequency" );
  opt->addUsage( "     --n-surfs               The number of surfaces" );
  opt->addUsage( "     --num-threads           The number OpenMP threads" );
  opt->addUsage( "     --position-iter-max     Maximum number of iterations for the position solver" );
  opt->addUsage( "     --position-restart      Maximum number of restarts for the position solver" );
  opt->addUsage( "     --position-tol          The tolerence for the position solver" );
  opt->addUsage( "     --rep-filter-freq       The filter freq for the reparametrization" );
  opt->addUsage( "     --rep-max-iter          Maximum number of reparametrization steps" );
  opt->addUsage( "     --rep-timestep          The Time step for the reparametrization" );
  opt->addUsage( "     --rep-tol               The absolute value tol on the velocity of reparametrization" );
  opt->addUsage( "     --rep-up-freq           The upsampling frequency for the reparametrization" );
  opt->addUsage( "     --save-stride           The frequency of saving to file (in time scale)" );
  opt->addUsage( "     --sh-order              The spherical harmonics order" );
  opt->addUsage( "     --singular-stokes       The scheme for the singular stokes evaluation" );
  opt->addUsage( "     --tension-iter-max      Maximum number of iterations for the tension solver" );
  opt->addUsage( "     --tension-restart       Maximum number of restarts for the tension solver" );
  opt->addUsage( "     --tension-tol           The tolerence for the tension solver" );
  opt->addUsage( "     --time-horizon          The time horizon of the simulation" );
  opt->addUsage( "     --time-iter-max         Maximum number of iteration for the choice of time stepper" );
  opt->addUsage( "     --time-precond          The type of preconditioner to use" );
  opt->addUsage( "     --time-scheme           The time stepping scheme" );
  opt->addUsage( "     --time-tol              The desired error tolerance in the time stepping" );
  opt->addUsage( "     --timestep              The time step size" );
  opt->addUsage( "     --upsample-interaction  Flag to whether upsample (and filter) the interaction force" );
  opt->addUsage( "     --viscosity-contrast    The viscosity contrast of vesicles" );
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
  //ordered alphabetically
  opt->setOption( "bending-modulus" );
  opt->setOption( "bg-flow-param" );
  opt->setOption( "bg-flow-type" );
  opt->setOption( "cent-file");
  opt->setOption( "error-factor" );
  opt->setOption( "filter-freq" );
  opt->setOption( "n-surfs" );
  opt->setOption( "num-threads" );
  opt->setOption( "position-iter-max" );
  opt->setOption( "position-restart" );
  opt->setOption( "position-tol" );
  opt->setOption( "rep-filter-freq" );
  opt->setOption( "rep-max-iter" );
  opt->setOption( "rep-timestep" );
  opt->setOption( "rep-tol" );
  opt->setOption( "rep-up-freq" );
  opt->setOption( "save-stride" );
  opt->setOption( "sh-order" );
  opt->setOption( "singular-stokes" );
  opt->setOption( "tension-iter-max" );
  opt->setOption( "tension-restart" );
  opt->setOption( "tension-tol" );
  opt->setOption( "time-horizon" );
  opt->setOption( "time-iter-max" );
  opt->setOption( "time-precond" );
  opt->setOption( "time-scheme" );
  opt->setOption( "time-tol" );
  opt->setOption( "timestep" );
  opt->setOption( "viscosity-contrast" );

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
    this->bending_modulus = atof(opt->getValue( "bending-modulus" ));

  if( opt->getValue( "viscosity-contrast" ) != NULL  )
    this->viscosity_contrast = atof(opt->getValue( "viscosity-contrast" ));

  if( opt->getValue( "bg-flow-type" ) != NULL  )
    this->bg_flow = EnumifyBgFlow(opt->getValue( "bg-flow-type" ));
    ASSERT(this->bg_flow != UnknownFlow, "Failed to parse the background flow name");

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

  if( opt->getValue( "num-threads" ) != NULL  )
    this->num_threads =  atoi(opt->getValue( "num-threads" ));

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

  if( opt->getValue( "time-scheme" ) != NULL  ){
    this->scheme = EnumifyScheme(opt->getValue( "time-scheme" ));
    ASSERT(this->scheme != UnknownScheme, "Failed to parse the time scheme name");
  }

  if( opt->getValue( "time-precond" ) != NULL  ){
    this->time_precond = EnumifyPrecond(opt->getValue( "time-precond" ));
    ASSERT(this->time_precond != UnknownPrecond, "Failed to parse the preconditioner name" );
  }

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

  if( opt->getValue( "time-tol" ) != NULL  )
    this->time_tol =  atof(opt->getValue( "time-tol" ));

  if( opt->getValue( "time-iter-max" ) != NULL  )
    this->time_iter_max =  atof(opt->getValue( "time-iter-max" ));

//   other methods: (bool) opt.getFlag( ... long or short ... )
}

template<typename T>
Error_t Parameters<T>::pack(std::ostream &os, Format format) const
{
    ASSERT(format==Streamable::ASCII, "BIN is not supported yet");

    os<<"PARAMETERS\n";
    os<<"n_surfs: "<<n_surfs<<"\n";
    os<<"sh_order: "<<sh_order<<"\n";
    os<<"filter_freq: "<<filter_freq<<"\n";
    os<<"bending_modulus: "<<bending_modulus<<"\n";
    os<<"viscosity_contrast: "<<viscosity_contrast<<"\n";
    os<<"position_solver_iter: "<<position_solver_iter<<"\n";
    os<<"tension_solver_iter: "<<tension_solver_iter<<"\n";
    os<<"position_solver_restart: "<<position_solver_restart<<"\n";
    os<<"tension_solver_restart: "<<tension_solver_restart<<"\n";
    os<<"position_solver_tol: "<<position_solver_tol<<"\n";
    os<<"tension_solver_tol: "<<tension_solver_tol<<"\n";
    os<<"time_horizon: "<<time_horizon<<"\n";
    os<<"ts: "<<ts<<"\n";
    os<<"time_tol: "<<time_tol<<"\n";
    os<<"time_iter_max: "<<time_iter_max<<"\n";
    os<<"scheme: "<<scheme<<"\n";
    os<<"time_precond: "<<time_precond<<"\n";
    os<<"bg_flow: "<< bg_flow<<"\n";
    os<<"singular_stokes: "<< singular_stokes<<"\n";
    os<<"rep_maxit: "<<rep_maxit<<"\n";
    os<<"rep_up_freq: "<<rep_up_freq<<"\n";
    os<<"rep_filter_freq: "<<rep_filter_freq<<"\n";
    os<<"rep_ts: "<<rep_ts<<"\n";
    os<<"rep_tol: "<<rep_tol<<"\n";
    os<<"bg_flow_param: "<<bg_flow_param<<"\n";
    os<<"upsample_interaction: "<<upsample_interaction<<"\n";
    os<<"save_data: "<<save_data<<"\n";
    os<<"save_stride: "<<save_stride<<"\n";
    os<<"init_file_name: "<<init_file_name<<" |\n"; //added | to stop >> from consuming next line if string is empty
    os<<"cntrs_file_name: "<<cntrs_file_name<<" |\n"; //added | to stop >> from consuming next line if string is empty
    os<<"save_file_name: "<<save_file_name<<" |\n"; //added | to stop >> from consuming next line if string is empty
    os<<"error_factor: "<<error_factor<<"\n";
    os<<"num_threads: "<<num_threads<<"\n";
    os<<"/PARAMETERS\n";
    return ErrorEvent::Success;
}

template<typename T>
Error_t Parameters<T>::unpack(std::istream &is, Format format)
{
    ASSERT(format==Streamable::ASCII, "BIN is not supported yet");
    std::string s, key;
    is>>s;
    ASSERT(s=="PARAMETERS", "Bad input string (missing header).");

    is>>key>>n_surfs;			ASSERT(key=="n_surfs:", "Unexpected key (expeded n_surfs)");
    is>>key>>sh_order;			ASSERT(key=="sh_order:", "Unexpected key (expected sh_order)");
    is>>key>>filter_freq;		ASSERT(key=="filter_freq:", "Unexpected key (expected filter_freq)");
    is>>key>>bending_modulus;		ASSERT(key=="bending_modulus:", "Unexpected key (expected bending_modulus)");
    is>>key>>viscosity_contrast;	ASSERT(key=="viscosity_contrast:", "Unexpected key (expected viscosity_contrast)");
    is>>key>>position_solver_iter;	ASSERT(key=="position_solver_iter:", "Unexpected key (expected position_solver_iter)");
    is>>key>>tension_solver_iter;	ASSERT(key=="tension_solver_iter:", "Unexpected key (expected tension_solver_iter)");
    is>>key>>position_solver_restart;	ASSERT(key=="position_solver_restart:", "Unexpected key (expected position_solver_restart)");
    is>>key>>tension_solver_restart;	ASSERT(key=="tension_solver_restart:", "Unexpected key (expected tension_solver_restart)");
    is>>key>>position_solver_tol;	ASSERT(key=="position_solver_tol:", "Unexpected key (expected position_solver_tol)");
    is>>key>>tension_solver_tol;	ASSERT(key=="tension_solver_tol:", "Unexpected key (expected tension_solver_tol)");
    is>>key>>time_horizon;		ASSERT(key=="time_horizon:", "Unexpected key (expected time_horizon)");
    is>>key>>ts;			ASSERT(key=="ts:", "Unexpected key (expected ts)");
    is>>key>>time_tol;			ASSERT(key=="time_tol:", "Unexpected key (expected time_tol)");
    is>>key>>time_iter_max;		ASSERT(key=="time_iter_max:", "Unexpected key (expected time_iter_max)");

    //enums
    is>>key>>s; 			ASSERT(key=="scheme:", "Unexpected key (expected scheme)");
    scheme=EnumifyScheme(s.c_str());
    is>>key>>s;          		ASSERT(key=="time_precond:", "Unexpected key (expected time_precond)");
    time_precond=EnumifyPrecond(s.c_str());
    is>>key>>s; 			ASSERT(key=="bg_flow:", "Unexpected key (expected  bg_flow)");
    bg_flow=EnumifyBgFlow(s.c_str());
    is>>key>>s;          		ASSERT(key=="singular_stokes:", "Unexpected key (expected  singular_stokes)");
    singular_stokes=EnumifyStokesRot(s.c_str());

    is>>key>>rep_maxit;			ASSERT(key=="rep_maxit:", "Unexpected key (expected rep_maxit)");
    is>>key>>rep_up_freq;		ASSERT(key=="rep_up_freq:", "Unexpected key (expected rep_up_freq)");
    is>>key>>rep_filter_freq;		ASSERT(key=="rep_filter_freq:", "Unexpected key (expected rep_filter_freq)");
    is>>key>>rep_ts;			ASSERT(key=="rep_ts:", "Unexpected key (expected rep_ts)");
    is>>key>>rep_tol;			ASSERT(key=="rep_tol:", "Unexpected key (expected rep_tol)");
    is>>key>>bg_flow_param;		ASSERT(key=="bg_flow_param:", "Unexpected key (expected bg_flow_param)");
    is>>key>>upsample_interaction;	ASSERT(key=="upsample_interaction:", "Unexpected key (expected upsample_interaction)");
    is>>key>>save_data;			ASSERT(key=="save_data:", "Unexpected key (expected save_data)");
    is>>key>>save_stride;		ASSERT(key=="save_stride:", "Unexpected key (expected save_stride)");

    is>>key>>s;		                ASSERT(key=="init_file_name:", "Unexpected key (expected init_file_name)");
    if (s!="|"){init_file_name=s; is>>s;}; //conume |
    is>>key>>s;          		ASSERT(key=="cntrs_file_name:", "Unexpected key (expected cntrs_file_name)");
    if (s!="|"){cntrs_file_name=s; is>>s;}; //conume |
    is>>key>>s;         		ASSERT(key=="save_file_name:", "Unexpected key (expected save_file_name)");
    if (s!="|"){save_file_name=s; is>>s;}; //conume |
    is>>key>>error_factor;		ASSERT(key=="error_factor:", "Unexpected key (expected error_factor)");
    is>>key>>num_threads;		ASSERT(key=="num_threads:", "Unexpected key (expected num_threads)");
    is>>s;
    ASSERT(s=="/PARAMETERS", "Bad input string (missing footer).");

    return ErrorEvent::Success;
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
    output<<"   viscosity contrast       : "<<par.viscosity_contrast<<std::endl;
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
    output<<"   Time step                : "<<par.ts<<std::endl;
    output<<"   Scheme                   : "<<par.scheme<<std::endl;
    output<<"   Time tol                 : "<<par.time_tol<<std::endl;
    output<<"   Time iter max            : "<<par.time_iter_max<<std::endl;
    output<<"   Precond                  : "<<par.time_precond<<std::endl;
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
    output<<"   Background flow type     : "<<par.bg_flow<<std::endl;
    output<<"   Background flow parameter: "<<par.bg_flow_param<<std::endl;
    output<<"   Upsample Interaction     : "<<std::boolalpha<<par.upsample_interaction<<std::endl;

    output<<"------------------------------------"<<std::endl;
    output<<" Misc:"<<std::endl;
    output<<"   OpenMP num threads       : "<<par.num_threads<<std::endl;
    output<<"====================================";

    return output;
}
