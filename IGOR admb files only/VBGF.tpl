DATA_SECTION

  !! ad_comm::change_datafile_name("VBGF.nam");       //  get filenames
//DAT file name
  init_adstring datfilename
  init_adstring ctlfilename
  init_int readparfile
  init_adstring version_info   // read base string from file
  !! adstring(version_info)+=adstring("1.0 (Oct_10_2005)");

  !! ad_comm::change_datafile_name(datfilename);
  !!cout<<"START DAT file"<<endl;

  init_int read		//Number of Readers/Reads
  init_int n		//Number of Individuals

//Age and growth data
  init_matrix ag_data_T(1,n,1,read+2);
  matrix ag_data(1,read+2,1,n);
  !! ag_data = trans(ag_data_T);
  vector lengths(1,n)
 LOCAL_CALCS
  lengths =column(ag_data_T,2);
 END_CALCS
// !!cout<<column(ag_data,1)<<endl<<column(ag_data,n)<<endl;
// !!cout<<ag_data<<endl;
 !!cout<<"END .dat file"<<endl;
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
  init_int fit_choice		//Fit by standard regression, functional regression, or random effects
  init_int reg_type		//Regression type (1/2)
  init_number d_LB		//Lower bound for contamination value
  init_number d_UB		//Upper bound for contamination value
  init_number d_start		//Starting contamination value
  init_number d_phase		//Estimating phase for contamination value
  init_number CV_Lt_start	//Starting CV_Lt value
  init_number CV_Lt_phase	//Estimation phase of the length CV
  init_int fit_data		//Regression fit to lengths, ages, or both
  init_int age_t_phase		//Estimation phase of functional regression age vector
  init_number lambda_in		//Input ratio of observation to process error
//  init_number lam_phase	//Phase of lambda estimation
  init_number sd_all_ages		//Overall standard deviation of ageing (error)

 !!cout<<"Starting values"<<endl;
 !!cout<<"Linf  "<<Linf_start<<endl;
 !!cout<<"k  "<<k_start<<endl;
 !!cout<<"t0  "<<t0_start<<endl;
 !!cout<<""<<endl;

 !!cout<<"Fit choice "<<fit_choice<<endl<<"Regression type "<<reg_type<<endl<<"Regression fit "<<fit_data<<endl;
 !!cout<<"END CTL file"<<endl;

  int i			//counter
  number read_i		//counter
  number j		//counter
  number read_mod	//counter

PARAMETER_SECTION
 LOCAL_CALCS
    if(readparfile>=1) 
    {cout<<"read par file"<<endl;
    ad_comm::change_pinfile_name("VBGF.PAR");}
 END_CALCS

//VBGF parameters
  init_bounded_number Linf(0.000001,10000,Linf_phase)	//Linf
  init_bounded_number k(0.0000001,10,k_phase)	//k
  init_number t0(t0_phase)	//t0
  init_bounded_number d(d_LB,d_UB,d_phase)	//robust contamination parameter
  init_bounded_number d_age(d_LB,d_UB,d_phase)	//robust contamination parameter
  init_bounded_number log_CV_Lt(-100,100,CV_Lt_phase)	//Length (process) error
//  init_number lambda (lam_phase)  
  init_vector ages_FR(1,n,age_t_phase)		//Functional regression ages
	
  vector Lt(1,n)		//Length vector
  vector age_t(1,n)		//Age vector

//Vector used in the functional regression
  //vector mean_age(1,n)			//Mean age for each individual across reads
  vector temp_diff(1,read_use)		//For SD calculations
  matrix age_calc(1,read_use,1,n)	//Age vector for calculations
  vector Lt_FR(1,n)			//Length fits from ages_FR
  number diff_sq_Lt			//SS for lengths
  vector diff_sq_age(1,read_use)	//SS for ages
  number lambda				//Ratio of observation to process error
  number sigma_CV_Lt			//Length (process) error as CV in likelihood
  number sigma_Lt			//Length (process) error as SD for the fxn reg.
  //vector sd_ages(1,n)			//Standard deviation of ages by readers

  objective_function_value obj_fun;

 LOCAL_CALCS

  if(Linf_start<max(ag_data(2)))
    {Linf=1.1*max(ag_data(2));}
  else 
    {Linf=Linf_start;}
  k = k_start;
  t0= t0_start;
  d = d_start;
  d_age = d;
  ages_FR = ag_data(3);
  log_CV_Lt = log(CV_Lt_start);

//  if(read_use>1)
//  {
//  for (read_i=1;read_i<=read_use;read_i++)
//    {age_calc(read_i) = extract_row(ag_data,2+read_i);}
//    mean_age=colsum(age_calc)/read_use;
//  for (j=1;j<=n;j++)
//   {sd_ages(j)= sqrt(value(norm2(column(age_calc,j)- mean_age(j))/(read_use-1)));}
//  sd_all_ages=sum(sd_ages)/(n-1);
//  sd_all_ages=0.1;
//  }

  sigma_Lt = 20;


 END_CALCS
//  !!cout<<"age calc"<<endl;
//  !!cout<<age_calc<<endl;
//  !!cout<<"mean age"<<endl;
//  !!cout<<mean_age<<endl;
//  !!cout<<"sd ages 1"<<endl;
//  !!cout<<sd_ages(1)<<endl;
//  !!cout<<"sd all ages"<<endl;
//  !!cout<<sd_all_ages<<endl;
//  !!cout<<""<<endl;
  !!cout <<"Parameters have been initialized-- END PARAMETER SECTION"<<endl;

PROCEDURE_SECTION
//Basic VBGF equations
//  Lt=Linf*(1-mfexp(-k*(ag_data(3)-t0)));
//  age_t= -(log(1-(Lt/Linf_age))/k_age)+t0_age;
  obj_fun = 0;
  if(fit_choice==1)
    {
     if(reg_type==1)
       {
//	1) Fit to primary read
//	2) Fit to average read
//	3) Fit to median read
// Fit to likelihood instead of least squares
	   if(fit_data<4)
             for (j = 1; j<=n; j++)
                  {Lt(j)=Linf*(1-mfexp(-k*(ag_data(3,j)-t0)));
                   sigma_CV_Lt = log(Lt(j))+log_CV_Lt;
                   obj_fun += sigma_CV_Lt + 0.5*square((Lt(j)-ag_data(2,j))/mfexp(sigma_CV_Lt));}
//Fit to regression function
	   else if(fit_data==4)
             {Lt=Linf*(1-mfexp(-k*(ag_data(3)-t0)));
              obj_fun = regression(ag_data(2),Lt);}
///Fit to multiple ages
	   else if(fit_data==5)
             for (read_i=1;read_i<=read_use;read_i++)
              {Lt=Linf*(1-mfexp(-k*(ag_data(2+read_i)-t0)));
              obj_fun += regression(ag_data(2),Lt);}
// Fit to Ages
	   else if(fit_data==6)
             {age_t= -(log(1-(ag_data(2)/Linf))/k)+t0;
             obj_fun = regression(ag_data(3),age_t);}
//Fit to lengths and ages
	   else if(fit_data==7)
             {Lt=Linf*(1-mfexp(-k*(ag_data(3)-t0)));
             obj_fun = regression(ag_data(2),Lt);
             age_t= -(log(1-(ag_data(2)/Linf))/k)+t0; 
             obj_fun += regression(ag_data(3),age_t);}
           else
             {abort();}
       }
     else if(reg_type==2)
       {
           if(fit_data<4)
             {Lt=Linf*(1-mfexp(-k*(ag_data(3)-t0)));
             obj_fun = robust_regression(ag_data(2),Lt,d);}
           else if(fit_data==5)
             {age_t= -(log(1-(ag_data(2)/Linf))/k)+t0;
             obj_fun = robust_regression(ag_data(3),age_t);}
           else if(fit_data==6)
             for (read_i=1;read_i<=read_use;read_i++)
              {Lt=Linf*(1-mfexp(-k*(ag_data(2+read_i)-t0)));
              obj_fun += robust_regression(ag_data(2+read_i),Lt,d);}
           else if(fit_data==7)
             {Lt=Linf*(1-mfexp(-k*(ag_data(3)-t0)));
             obj_fun = robust_regression(ag_data(2),Lt,d);
             age_t= -(log(1-(ag_data(2)/Linf))/k)+t0; 
             obj_fun += robust_regression(ag_data(3),age_t);}
           else
             {abort();}
       }
     else
       {abort();}
    }
  else if(fit_choice==2)
   {
     Lt=Linf*(1-mfexp(-k*(ages_FR-t0)));
     diff_sq_Lt = norm2(ag_data(2)-Lt);
     if(lambda_in>0)
     {lambda = lambda_in;}
     else if (lambda_in<=0)
     sigma_Lt = sqrt(diff_sq_Lt/(n-1));
     {lambda=(sigma_Lt*sigma_Lt)/(sd_all_ages*sd_all_ages);}
     for (read_i=1;read_i<=read_use;read_i++)
       {diff_sq_age(read_i) = lambda*norm2(ag_data(2+read_i)-ages_FR);}
     obj_fun= diff_sq_Lt+sum(diff_sq_age);
   }

//  cout << obj_fun << endl;

REPORT_SECTION
  report<<"Version 1.0 (Oct_10_2005)"<<endl;
  report<< "Objective function" << endl;
  report<< obj_fun <<endl;
  report<< "Linf" <<endl;
  report<<Linf<<endl;
  report<< "k" <<endl;
  report<<k<<endl;
  report<< "t0" <<endl;
  report<<t0<<endl;
  if(fit_choice==1)
  {
    if(fit_data< 4)
      {
       report<< "CV_Lt"<<endl;
       report<<mfexp(log_CV_Lt)<<endl;
      }
    if(fit_data < 6 || fit_data==7)
      {
        report<<"Predicted Lengths (Lt)"<<endl;
        report<<Lt<<endl;
      }
    if(fit_data== 6 || fit_data==7)
      {
        report<<"Predicted Ages"<<endl;
        report<<age_t<<endl;
      }
  }
  if(fit_choice==1 && fit_data ==2)
  {
   report<< "d"<<endl;
   report<<d<<endl;
   report<< "d_age"<<endl;
   report<<d_age<<endl;
   }
  if(fit_choice==2)
  {     report<<"Lambda"<<endl;
        report<<lambda<<endl;
	report<<"sigma_Lt"<<endl;
	report<<sigma_Lt<<endl;
	report<<"sigma_age"<<endl;
	report<<sd_all_ages<<endl;}

//*======================================================================================*
TOP_OF_MAIN_SECTION
  arrmblsize = 50000000;
  gradient_structure::set_GRADSTACK_BUFFER_SIZE(10000000000);
  gradient_structure::set_CMPDIF_BUFFER_SIZE(100000000);
  gradient_structure::set_MAX_NVAR_OFFSET(5000);
  gradient_structure::set_NUM_DEPENDENT_VARIABLES(5000);

GLOBALS_SECTION
  #include <admodel.h>
  #include <time.h>
  time_t start,finish;
  long hour,minute,second;
  double elapsed_time;

FINAL_SECTION
  //Calculates model run time
  time(&finish);
  elapsed_time = difftime(finish,start);
  hour = long(elapsed_time)/3600;
  minute = long(elapsed_time)%3600/60;
  second = (long(elapsed_time)%3600)%60;
  cout<<endl<<endl<<"starting time: "<<ctime(&start);
  cout<<"finishing time: "<<ctime(&finish);
  cout<<"This run took: ";
  cout<<hour<<" hours, "<<minute<<" minutes, "<<second<<" seconds."<<endl<<endl<<endl;
