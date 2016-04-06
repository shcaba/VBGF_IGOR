DATA_SECTION

  !! ad_comm::change_datafile_name("VBGF_RE.nam");       //  get filenames
//DAT file name
  init_adstring datfilename
  init_adstring ctlfilename
  init_int readparfile
  init_adstring version_info   // read base string from file
  !! adstring(version_info)+=adstring("1.0 [Oct_10_2005]");

  !! ad_comm::change_datafile_name(datfilename);
  !!cout<<"START DAT file"<<endl;

  init_int read		//Number of Readers/Reads
  init_int n		//Number of Individuals


//Age and growth data
  init_matrix ag_data_T(1,n,1,read+2);
  matrix ag_data(1,read+2,1,n);
  !! ag_data = trans(ag_data_T);
   number eps
  !! eps = 1e-5;

  int i			//Counter for SF1
  int read_i		//counter for samples

 //!!cout<<column(ag_data,1)<<endl<<column(ag_data,n)<<endl;
 !!cout<<""<<endl;
 !!cout<<"END .DAT file"<<endl;
 !!cout<<""<<endl;

//CTL file inputs
  !! ad_comm::change_datafile_name(ctlfilename);
  !!cout<<"START CTL file"<<endl;

  init_int Linf_start	 //Starting estimate of Linf
  init_int Linf_phase	 //Phase estimation of Linf
  init_number k_start	 //Starting estimate of Linf
  init_int k_phase	 //Phase estimation of Linf
  init_number t0_start	 //Starting estimate of Linf
  init_int t0_phase	 //Phase estimation of Linf

  init_int read_use		//Number of reads to use
  init_number CV_e		//CV on the random effects (not estimated)
  init_int age_dist		//Age likelihood
  init_number age_re_phase	//Phase for random effects

  init_number CV_Lt_start	//Starting value for the variance on length (sigma)
  init_number CV_Lt_phase	//Phase for sigma estimation

  init_number alpha_start	//Starting value for mean random effects
  init_number alpha_phase	//Phase for alpha estimation
  init_number sigma_age_start	//Starting value for the variance on the random effects (sigma_a)
  init_number sigma_age_phase	//Phase for sigma_a estimation

  init_number beta_start	//Starting value for exponential parameter
  init_number beta_phase	//Phase for exponential parameter

  init_number gam_shape_start	//Starting value for exponential parameter
  init_number gam_shape_phase	//Phase for exponential parameter

  init_number gam_scale_start	//Starting value for exponential parameter
  init_number gam_scale_phase	//Phase for exponential parameter

  //  !!cout<<"Starting values"<<endl;
  //  !!cout<<"Linf  "<<Linf_start<<endl;
  //  !!cout<<"k  "<<k_start<<endl;
  //  !!cout<<"t0  "<<t0_start<<endl;
  //  !!cout<<""<<endl;
  //  
  //  !!cout<<"k_phase  "<<k_phase<<endl;
  //  !!cout<<"read_use  "<<read_use<<endl;
  //  !!cout<<"CV_e  "<<CV_e<<endl;
  //  !!cout<<"sigma_Lt_start  "<<sigma_Lt_start<<endl;
  //  !!cout<<"sigma_Lt_phase  "<<sigma_Lt_phase<<endl;
  //  !!cout<<"sigma_re_phase  "<<age_re_phase<<endl;
  //  !!cout<<"age_dist  "<<age_dist<<endl;

 // !!cout<<"END CTL file"<<endl;

  number Linf_lo

PARAMETER_SECTION
 LOCAL_CALCS
    if(readparfile>=1) 
    {cout<<"read par file"<<endl;
    ad_comm::change_pinfile_name("VBGF_RE.PAR");}
    Linf_lo =max(ag_data(2))*0.75;

 END_CALCS

//VBGF parameters
  init_bounded_number Linf(Linf_lo,Linf_lo*4,Linf_phase)		//Linf
  init_bounded_number k(0.0001,1,k_phase)				//k
  init_bounded_number t0(-15,0,t0_phase)				//t0
  init_bounded_number log_CV_Lt(-10,2,CV_Lt_phase)			//Process error in lengths
  init_bounded_number alpha(0,500,alpha_phase)					//Mean age
  init_bounded_number log_sigma_age(-100,100,sigma_age_phase)		//Observation error in ages
  init_bounded_number log_beta(-100,100,beta_phase)			//Exponential rate parameter
  init_bounded_number gam_shape(0,100,gam_shape_phase)			//Exponential rate parameter
  init_bounded_number gam_scale(0,100,gam_scale_phase)			//Exponential rate parameter
  vector start_re(1,n)

  random_effects_vector age_re(1,n,age_re_phase)	//Random effects age vector

  objective_function_value obj_fun;

PRELIMINARY_CALCS_SECTION
  if(Linf_start<max(ag_data(2)))
    {Linf=Linf_start;}
  else 
    {Linf=Linf_start;}
  k= k_start;
  t0= t0_start;
  log_CV_Lt = log(CV_Lt_start);
  alpha = alpha_start;
  log_sigma_age = log(sigma_age_start);
  log_beta = log(beta_start);
  gam_shape = gam_shape_start;
  gam_scale = gam_scale_start;
  if (CV_e < 0.05)
  {CV_e = 0.05;}
  cout<<"Starting values"<<endl;
  cout<<"Linf  "<<Linf<<endl;
  cout<<"k  "<<k<<endl;
  cout<<"t0  "<<t0<<endl;
  cout<<""<<endl;
  cout<<""<<endl;
  cout<<"CV_e"<<CV_e<<endl;
  cout <<"Parameters have been initialized-- END PARAMETER SECTION"<<endl;
  cout<<""<<endl;

PROCEDURE_SECTION

  obj_fun = 0;

  if(age_dist==1)
     {
     for (i =1;i<=n;i++)
     {sf1(Linf, k, age_re(i), t0, log_CV_Lt, log_sigma_age, alpha);}
     }

  if(age_dist==2)
  {
     for (i =1;i<=n;i++)
     {sf2(Linf, k, age_re(i), t0, log_CV_Lt, log_beta);}
  }

  if(age_dist==3)
  {
     for (i =1;i<=n;i++)
     {sf3(Linf, k, age_re(i), t0, log_CV_Lt, gam_shape, gam_scale);}
  }

SEPARABLE_FUNCTION void sf1(const dvariable& Linf_i, const dvariable& k_i, const dvariable& age_re_i, const dvariable& t0_i, const dvariable& log_CV_Lt_i, const dvariable& sigma_age_i, const dvariable& alpha_i)
       dvariable Lt_i = Linf_i*(1-mfexp(-k_i*(age_re_i-t0_i)));
       dvariable sigma_e_i;
       dvariable sigma_Lt_i;
       if(current_phase()==1)
       {
       sigma_e_i = CV_e*(ag_data(3,i)+age_re_i);
       }
       if(current_phase()>1)
       {
       sigma_e_i= (CV_e*age_re_i+eps);
       sigma_Lt_i = log_CV_Lt_i+log(Lt_i+eps);
       }
       obj_fun += sigma_Lt_i + 0.5*square((Lt_i-ag_data(2,i))/mfexp(sigma_Lt_i));
       obj_fun += sigma_age_i + 0.5*square((age_re_i-alpha_i)/mfexp(sigma_age_i));
       for (read_i=1;read_i<=read_use;read_i++)
       if (ag_data(2+read_i,i)>=0)
       {
         {obj_fun += log(sigma_e_i) + 0.5*square((ag_data(2+read_i,i)-age_re_i)/sigma_e_i);}
	}

SEPARABLE_FUNCTION void sf2(const dvariable& Linf_i, const dvariable& k_i, const dvariable& age_re_i, const dvariable& t0_i,const dvariable& log_CV_Lt_i, const dvariable& log_beta_i)
       dvariable Lt_i = Linf_i*(1-mfexp(-k_i*(age_re_i-t0_i)));
       dvariable sigma_e_i;
       dvariable sigma_Lt_i;
       if(current_phase()==1)
       {
       sigma_e_i = CV_e*(ag_data(3,i)+age_re_i);
       }
       if(current_phase()>1)
       {
       sigma_e_i= (CV_e*age_re_i+eps);
       sigma_Lt_i = log_CV_Lt_i+log(Lt_i+eps);
       }
       obj_fun += sigma_Lt_i + 0.5*square((Lt_i-ag_data(2,i))/mfexp(sigma_Lt_i));
       obj_fun += -log_beta_i+(mfexp(log_beta_i)*age_re_i);
       for (read_i=1;read_i<=read_use;read_i++)
       if (ag_data(2+read_i,i)>=0)
       {
         {obj_fun += log(sigma_e_i) + 0.5*square((ag_data(2+read_i,i)-age_re_i)/(sigma_e_i));}
	}

SEPARABLE_FUNCTION void sf3(const dvariable& Linf_i, const dvariable& k_i, const dvariable& age_re_i, const dvariable& t0_i,const dvariable& log_CV_Lt_i, const dvariable& gam_shape_i, const dvariable& gam_scale_i)
       dvariable Lt_i = Linf_i*(1-mfexp(-k_i*(age_re_i-t0_i)));
       dvariable sigma_e_i;
       dvariable sigma_Lt_i;
       if(current_phase()==1)
       {
       sigma_e_i = CV_e*(ag_data(3,i)+age_re_i);
       }
       if(current_phase()>1)
       {
       sigma_e_i= (CV_e*age_re_i+eps);
       sigma_Lt_i = log_CV_Lt_i+log(Lt_i+eps);
       }
       obj_fun += sigma_Lt_i + 0.5*square((Lt_i-ag_data(2,i))/mfexp(sigma_Lt_i));
       obj_fun += -((gam_shape_i-1)*log((age_re_i+eps)/gam_scale_i))+((age_re_i+eps)/gam_scale_i)+log(gam_scale_i)+gammln(gam_shape_i);
       for (read_i=1;read_i<=read_use;read_i++)
       if (ag_data(2+read_i,i)>=0)
       {
         {obj_fun += log(sigma_e_i) + 0.5*square((ag_data(2+read_i,i)-age_re_i)/(sigma_e_i));}
	}

REPORT_SECTION
  report<<"Version 1.0 [Oct_10_2005]"<<endl;
  report<< "Objective function" << endl;
  report<< obj_fun <<endl;
  report<< "Linf" <<endl;
  report<<Linf<<endl;
  report<< "k" <<endl;
  report<<k<<endl;
  report<< "t0" <<endl;
  report<<t0<<endl;
  report<<"CV_Lt"<<endl;
  report<<mfexp(log_CV_Lt)<<endl;
  if(age_dist==1)
  {report << "alpha" <<endl;
   report << alpha <<endl;
   report<< "log_sigma_age" <<endl;
   report<< log_sigma_age <<endl;}
  if(age_dist==2)
   {report<< "log_beta" <<endl;
   report<< log_beta <<endl;}
  if(age_dist==3)
   {report<< "gam_shape" <<endl;
   report<< gam_shape <<endl;
   report<< "gam_scale" <<endl;
   report<< gam_scale <<endl;}
  report<<"age_re"<<endl;
  report<<age_re<<endl;


//*======================================================================================*
TOP_OF_MAIN_SECTION
  arrmblsize = 50000000;
  gradient_structure::set_GRADSTACK_BUFFER_SIZE(1000000);
  gradient_structure::set_CMPDIF_BUFFER_SIZE(1000000);
  gradient_structure::set_MAX_NVAR_OFFSET(2394763);
  gradient_structure::set_NUM_DEPENDENT_VARIABLES(5000);

GLOBALS_SECTION
  #include <admodel.h>
  #include <time.h>
  time_t start,finish;
  long hour,minute,second;
  double elapsed_time;

FINAL_SECTION
  //Calculates how long is taking to run
  // this code is based on the Widow Rockfish model (from Erik H. Williams, NMFS-Santa Cruz) 
  time(&finish);
  elapsed_time = difftime(finish,start);
  hour = long(elapsed_time)/3600;
  minute = long(elapsed_time)%3600/60;
  second = (long(elapsed_time)%3600)%60;
  cout<<endl<<endl<<"starting time: "<<ctime(&start);
  cout<<"finishing time: "<<ctime(&finish);
  cout<<"This run took: ";
  cout<<hour<<" hours, "<<minute<<" minutes, "<<second<<" seconds."<<endl;
