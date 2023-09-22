// Statistical catch-at-age model for 4T Spring Herring
// Acoustic index is in 1000t, expect q < 0.001
// model scale 1000 fish. input WAA in kg



// this version of the model allows projections for different M and recruitment scenarios







DATA_SECTION
  init_int dum;

  // Control parameters -- switch to estimation control file
  !! ad_comm::change_datafile_name("scaControl.ctl");
  // Total number of years for the model.
  init_int nT;
  // year rwM starts
  init_int rwT;
  // Number of fisheries
  init_int nFisheries;
  // First age class
  init_int firstAge;
  // Plus group age...last age class
  init_int plusGroupAge;
  vector age(firstAge,plusGroupAge);
  !! age.fill_seqadd(2,1);
  // Number of M groups
  init_int nSplitAge;
  // Lower ages of M groups
  init_vector splitAge(1,nSplitAge);
  // Initial natural mortality rates
  init_int ph_M;
  // M deviations
  init_int ph_mDevs;
  // Prior distributions
  // Natural mortality priors: initial Ms
  init_vector prior_initM(1,nSplitAge);    // means
  init_vector priorSD_initM(1,nSplitAge);  // std devs
  // Natural mortality priors: random walk Ms
  init_vector priorSD_rwM(1,nSplitAge);    // std devs


  // First age for avg'ing F
  init_int a1F;
  // Last age for avg'ing F
  init_int a2F;
  // Recruitment model type: 0=const avgR; 1=const recRate
  init_int recType;
  init_int ph_log_recRate;
  !! if( recType==1 ) ph_log_recRate = 1;
  // Baranov iterations
  init_int baranovIter;

  // Estimation phase numbers for model parameters
  // Unfished biomass
  init_int ph_log_AvgR;
  // Recruitment deviations
  init_int ph_initRecDevs;
  init_int ph_recDevs;
  // Number of selectivity time-blocks
  init_ivector nSelBlocks(1,nFisheries);
  // Selectivity time block boundaries
  init_matrix tBlock(1,nFisheries,1,nSelBlocks);
  // Selectivity logS50
  init_int ph_S50;
  // Selectivity logS95
  init_int ph_S95;
  // Recruitment autcorrelation
  init_int ph_gammaR;

  // Recruitment deviation prior std error
  init_number priorSD_R
  init_number QdevSD;
  init_int phzQdev;
  init_int QdevFirst;
  init_int QdevLast;

  // Catch levels for predictions
  init_vector Cpro(1,4);
  // end year for projection, pT
  init_int pT;
  // Simulation flag: 0=NO, 1=Yes
  init_int simFlag;
  // Simulation seed
  init_int simSeed;

  init_int eof1;


 LOC_CALCS
   if(eof1!=1111)
   {
     cout<<" data entry error control file " <<endl;
     cout<<" data entry error priorSD_R "<< priorSD_R <<endl;
     cout<<" data entry error baranovIter "<< baranovIter <<endl;
     exit(1);
   }
 END_CALCS



  // Catch biomass data -- switch to catch data file
  !! ad_comm::change_datafile_name("scaCatch.dat");
  // Matrix of landed catch biomass (tonnes) values.
  init_matrix landCatchMatrix(1,nT,1,nFisheries+2);
  // Total landed catch
  vector totLandCatch(1,nT);
  // Catch proportions by fishery
  vector propCatchFish(1,nFisheries);
  // Catch by fishery
  matrix  Ctg(1,nT,1,nFisheries);
  init_int eof2;


 LOC_CALCS
   if(eof2!=2222)
   {
     cout<<" data entry error catch file " <<endl;
     exit(1);
   }
 END_CALCS

  // Abundance index data -- switch to index data file
  !! ad_comm::change_datafile_name("scaIndex.dat");
  // Number of stock indices
  init_int nIndexSeries;
  //  Series index numbers.
  init_ivector idxIndex(1,nIndexSeries);
  // Is index relative? (0=false, 1=true).
  init_vector idxRelative(1,nIndexSeries);
  // Index weights in overall likelihood
  init_vector idxLikeWeight(1,nIndexSeries);
  // Vector indicating first year for each index.
  init_ivector idxFirstYear(1,nIndexSeries);
  // Vector indicating last year for each index.
  init_ivector idxLastYear(1,nIndexSeries);
  // Fraction of year when survey occurs
  init_vector fracYearSurvey(1,nIndexSeries);
  // Index series data
  init_matrix idxSeries(1,nIndexSeries,1,nT);
  // Create vector of max index series lengths.
  ivector nobs(1,nIndexSeries);
  !!nobs=idxLastYear-idxFirstYear+1;
  vector validObs(1,nIndexSeries);
  init_int eof3;


 LOC_CALCS
   if(eof3!=3333)
   {
     cout<<" data entry error index file " <<endl;
     cout<<" data entry error eof3"<< eof3 <<endl;
     cout<<" data entry error fracYearSurvey "<< fracYearSurvey <<endl;
     exit(1);
   }
 END_CALCS

  // Age-related data -- prop-age, weight-age, mat-age
  !! ad_comm::change_datafile_name("scaAges.dat");
  // Number of age series and their index numbers.
  init_int nAgeSeries;
  // Series index numbers
  init_ivector ageIndex(1,nAgeSeries);
  // First age class for observed paa
  init_ivector minAge(1,nAgeSeries);
  // Last age class
  init_ivector maxAge(1,nAgeSeries);
  // Vector indicating first year for each index.
  init_ivector ageFirstYear(1,nAgeSeries);
  // Vector indicating last year for each index.
  init_ivector ageLastYear(1,nAgeSeries);
  // matrix starting ages for paa
  init_matrix pst(1,nAgeSeries,1,nT);
  // matrix last ages for paa
  init_matrix pnd(1,nAgeSeries,1,nT);

  //Observed proportions-at-age
  init_3darray ageObsProp(1,nAgeSeries,firstAge,plusGroupAge,1,nT);

  // Beginning-of-year Weight-at-age by year
  init_matrix BwtAge(firstAge,plusGroupAge,1,nT);
  // Commercial fishery Weight-at-age by year
  init_matrix CwtAge(firstAge,plusGroupAge,1,nT);
  // cpue and acoustic Weight-at-age by year
  init_3darray IwtAge(1,nIndexSeries,firstAge,plusGroupAge,1,nT);
  // Maturity-at-age
  init_vector matAge(firstAge,plusGroupAge);
  // Initial values for age observation std errors
  init_vector init_tauAge(1,nAgeSeries);
  // M
  init_number M;
  // eof4
  init_int eof4;
  // tracking number of function evaluations
  int neval;


 LOC_CALCS
   if(eof4!=4444)
   {
     cout<<" data entry error age file " <<endl;
     cout<<" data entry error init_tauAge "<< init_tauAge <<endl;
     cout<<" data entry error BwtAge2 "<< BwtAge(2)(1,5) <<endl;
     cout<<" data entry error waaCPUE 11 "<< IwtAge(1)(11)(1,5) <<endl;
     cout<<" data entry error waaAcou 11 "<< IwtAge(2)(11)(1,5) <<endl;
     exit(1);
   }
   //cout<<" data read test " <<endl;

 END_CALCS

  //projection variables for output
    vector Mpro(firstAge,plusGroupAge);    // average of last 5 yr
    vector Sapro1(nT+1,pT);
    vector Sapro2(nT+1,pT);
    vector Sapro3(nT+1,pT);
    vector Sapro4(nT+1,pT);
    vector Bpro1(nT+1,pT);
    vector Bpro2(nT+1,pT);
    vector Bpro3(nT+1,pT);
    vector Bpro4(nT+1,pT);
    matrix Fpro(1,4,nT+1,pT); //  Ftg
    matrix Npro1(firstAge,plusGroupAge,nT+1,pT);
    matrix Npro2(firstAge,plusGroupAge,nT+1,pT);
    matrix Npro3(firstAge,plusGroupAge,nT+1,pT);
    matrix Npro4(firstAge,plusGroupAge,nT+1,pT);

  //!!cout<<Npro3;
  //!!exit(99);


//*******************************************************************/
PARAMETER_SECTION
  objective_function_value objFunc;
  // Average recruitment: initialized at 16 based on gbcod1fx2.par
  init_number log_AvgR(ph_log_AvgR);
  // Recruitment deviations: initial abundance 1978
  init_bounded_vector init_recDevs(firstAge,plusGroupAge,-5.,5.,ph_initRecDevs);
  // Recruitment deviations: 1978-2011
  init_bounded_vector recDevs(2,nT,-5.,5.,ph_recDevs);
  // Selectivity parameters
  init_bounded_matrix log_S50_gi(1,nFisheries,1,nSelBlocks,0.6,2.,ph_S50);
  init_bounded_matrix log_S95_step_gi(1,nFisheries,1,nSelBlocks,0.0,2.0,ph_S95);
  // Initial fishing mortality prior to 1978
  init_number Finit;
  // Recruitment autocorrelation
  init_number logit_gamma_R(ph_gammaR);
  // deviations in CPUE q
  init_bounded_vector Qdev(QdevFirst,QdevLast,-5,5,phzQdev);
  init_bounded_vector log_q(1,nIndexSeries,-12,1.5,1);
  init_bounded_vector log_sig(1,nIndexSeries,-3.0,1.5,4);
  // Initial natural mortality rates
  init_vector log_M(1,nSplitAge,ph_M);
  // Natural mortality ranWalk deviations
  init_bounded_matrix mDevs(1,nSplitAge,rwT,nT,-2.,2.,ph_mDevs);

  // Parameters: natural scale
  number avgR;
  number gamma_R;
  matrix Mt(1,nSplitAge,1,nT);
  matrix Mta(1,nT,firstAge,plusGroupAge);
  // Selectivity parameters: age-at-50%
  matrix S50_gt(1,nFisheries,1,nT);
  // Selectivity parameters: age-at-95%
  matrix S95_gt(1,nFisheries,1,nT);

  // Posterior quantities
  number recPrior;
  number qPrior;
  number mPrior;
  vector indexLikelihood(1,nIndexSeries);
  vector ageLikelihood(1,nAgeSeries);
  vector tauSquareAges(1,nAgeSeries);
  vector validObs(1,nIndexSeries);
  matrix q(1,nIndexSeries,1,nT);
  vector ss(1,nIndexSeries);
  vector sigIndex(1,nIndexSeries);

  // Derived variables
  // Unfished SSB-per-recruit
  number phiSSB;
  vector rw_recDevs(1,nT);
  // change recRate to a number parameter if recType is 1
  vector recRate(3,nT);
  vector SSBt(1,nT);
  vector SSBa(1,nT);
  matrix predNdx(1,nIndexSeries,1,nT);
  vector N4p(1,nT);
  vector B510(1,nT);
  matrix Nta(1,nT,firstAge,plusGroupAge);
  matrix Bta(1,nT,firstAge,plusGroupAge);
  matrix Zta(1,nT,firstAge,plusGroupAge);
  matrix Ftg(1,nT,1,nFisheries);
  matrix expBtg(1,nIndexSeries,1,nT);

  sdreport_vector avgF(1,nT);

  matrix  Bprime(1,nFisheries,firstAge,plusGroupAge);
  3darray uCgat(1,nAgeSeries,firstAge,plusGroupAge,1,nT);

  // Selectivity by age
  //matrix sel(1,nFisheries,1,plusGroupAge);
  3darray sel_gta(1,nFisheries,1,nT,firstAge,plusGroupAge);

 //*******************************************************************/
GLOBALS_SECTION
  void solveBaranov( const int& t, const int& iter);

  // Flag to control whether header written to mcmc output file.
  int mcmcHeader = 0;

  // Flag to control which parameters are done in MCMC phase.
  int mcmcFlag = 1;

  // Uncomment this line for 32-bit compiler.
  // #include <fstream.h>
  #include <admodel.h>
  ofstream mcoutSSB("mcoutSSB.dat");
  ofstream mcoutSSBa("mcoutSSBa.dat");
  ofstream mcoutF("mcoutF.dat");
  ofstream mcoutN2("mcoutN2.dat");
  ofstream mcoutN4("mcoutN4.dat");
  ofstream mcoutN4p("mcoutN4p.dat");
  ofstream mcoutPredC("mcoutPredC.dat");
  ofstream mcoutPredA("mcoutPredA.dat");
  ofstream mcoutq("mcoutq.dat");
  ofstream mcoutM1("mcoutM1.dat");
  ofstream mcoutM2("mcoutM2.dat");
  ofstream mcoutRp1("mcoutRp1.dat");
  ofstream mcoutRp2("mcoutRp2.dat");
  ofstream mcoutRp3("mcoutRp3.dat");
  ofstream mcoutRp4("mcoutRp4.dat");
  ofstream mcoutSpa1("mcoutSpa1.dat");
  ofstream mcoutSpa2("mcoutSpa2.dat");
  ofstream mcoutSpa3("mcoutSpa3.dat");
  ofstream mcoutSpa4("mcoutSpa4.dat");
  ofstream mcoutBp1("mcoutBp1.dat");
  ofstream mcoutBp2("mcoutBp2.dat");
  ofstream mcoutBp3("mcoutBp3.dat");
  ofstream mcoutBp4("mcoutBp4.dat");
  ofstream mcoutFp1("mcoutFp1.dat");
  ofstream mcoutFp2("mcoutFp2.dat");
  ofstream mcoutFp3("mcoutFp3.dat");
  ofstream mcoutFp4("mcoutFp4.dat");
  ofstream sim("scaSim.rep");

  int iseed=1337;
  random_number_generator rng(iseed);

TOP_OF_MAIN_SECTION
  arrmblsize=20000000;
  gradient_structure::set_CMPDIF_BUFFER_SIZE(25000000);
  gradient_structure::set_GRADSTACK_BUFFER_SIZE(1000000);

PRELIMINARY_CALCS_SECTION
  for ( int t=1;t<=landCatchMatrix.rowmax();t++ )
  {
    // total catch by fishery: Rows are different fisheries
    // and surveys. Columns are tSteps.
    for ( int i=3;i<=landCatchMatrix.colmax();i++ )
    {
        Ctg(t,i-2) = landCatchMatrix(landCatchMatrix(t,2),i);
    }
  }
  neval = 0;
  tauSquareAges  = square(init_tauAge);


PROCEDURE_SECTION
  //cout << "start procedure section..." << endl;
  initModelParameters();
  //cout << "initModelParameters done..." << endl;
  popInit();
  //cout << "popInit done..." << endl;
  popDynamics();
  //cout << "popDynamics done..." << endl;
  calc_index_likelihood();
  //cout << "index_likelihood done..." << endl;

  calc_age_likelihood();
  //cout << "age_likelihood done..." << endl;
  calc_rec_prior();
  //cout << "rec_prior done..." << endl;
  calc_M_prior();
  //cout << "M_prior done..." << endl;
  calc_q_prior();
  //cout << "M_prior done..." << endl;

  objFunc  = 0.;
  objFunc  = sum( elem_prod(idxLikeWeight, indexLikelihood) ) + sum( ageLikelihood );
  objFunc += recPrior;
  objFunc += qPrior;
  objFunc += mPrior;


  if ( mceval_phase() )
  {
    projection();
    mcoutSSB << SSBt    << endl;
    mcoutSSBa << SSBa    << endl;
    mcoutF   << avgF   << endl;
    mcoutN2 << column(Nta,2)    << endl;
    mcoutN4 << column(Nta,4)    << endl;
    mcoutN4p << N4p    << endl;
    mcoutPredC << predNdx(1)(13,45)    << endl;
    mcoutPredA << predNdx(2)(17,45)    << endl;
    mcoutq << q(1)(13,45) << endl;
    mcoutM1   << mfexp(Mt(1)) << endl;
    mcoutM2   << mfexp(Mt(2)) << endl;
    mcoutRp1   << Npro1(2)   << endl;
    mcoutRp2   << Npro2(2)   << endl;
    mcoutRp3   << Npro3(2)   << endl;
    mcoutRp4   << Npro4(2)   << endl;
    mcoutSpa1   << Sapro1   << endl;
    mcoutSpa2   << Sapro2   << endl;
    mcoutSpa3   << Sapro3   << endl;
    mcoutSpa4   << Sapro4   << endl;
    mcoutBp1   << Bpro1   << endl;
    mcoutBp2   << Bpro2   << endl;
    mcoutBp3   << Bpro3   << endl;
    mcoutBp4   << Bpro4   << endl;
    mcoutFp1   << Fpro(1)   << endl;
    mcoutFp2   << Fpro(2)   << endl;
    mcoutFp3   << Fpro(3)   << endl;
    mcoutFp4   << Fpro(4)   << endl;
    mcmcHeader = 1;
  }

FUNCTION initModelParameters
  {
  avgR    = mfexp( log_AvgR );

  gamma_R = exp( logit_gamma_R )/(1.+exp(logit_gamma_R));
  //cout << "avgR = " << avgR << endl;
  Mta.initialize();
  Mt.initialize();
  calc_Mta_1Split();
  //cout << "M_calc done" << endl;
  calc_sel_gta();
  //cout << "sel_calc done" << endl;
  get_q();
  //cout << "q_calc done" << endl;
  calcPhi();
  //cout << "Phi_calc done" << endl;
  sigIndex = mfexp(log_sig);
  }

FUNCTION get_q
  int i, t;
  for (i=1; i<=2; i++) q(i) = mfexp(log_q(i));
    if(active(Qdev))
  {
      for (t=QdevFirst; t<=QdevLast; t++) q(1,t)=q(1,t-1)*mfexp(Qdev(t));
  }


FUNCTION calcPhi
  {
  // SPC: calculating unfished SSB-per-recruit. Borrowing
  // Nta and Bta for now. Initializing with 1 recruit
  // to get through calcs.
  Nta.initialize(); SSBt.initialize();
  Nta(1,firstAge) = 1.;
  for( int a=firstAge+1; a<=(plusGroupAge-1); a++ )
  {
    Nta(1,a) = Nta(1,a-1)*exp( -M);
  }
  Nta(1,plusGroupAge) = Nta(1,plusGroupAge-1)*exp(-M)/(1.-exp(-M));
  Bta(1) = elem_prod( Nta(1), column(BwtAge,1) ); // units here are units of Wat
  phiSSB = sum( elem_prod( Bta(1),matAge ) );
  // Unfished spawning biomass.
  SSBt(1) = avgR*phiSSB;
  //cout << "phiSSB  = " << phiSSB << endl;
  //cout << "SSBt(1) = " << SSBt(1) << endl;
  }


FUNCTION calc_Mta_1Split
  {
  // Random walk for 2 age blocks for M
  // Build an Mta matrix that inclues age and time effects.
  // Set initial Mt's in years 1:13. All on log scale until the end.
  for( int t=1; t<rwT; t++ )
  {
    Mt(1)(t)      = log_M(1);
    Mt(2)(t)      = log_M(2);
    Mta(t)(firstAge,splitAge(2)-1)   = exp( Mt(1,t) );
    Mta(t)(splitAge(2),plusGroupAge) = exp( Mt(2,t) );
  }
  // Fill random walk.
  for( int t=rwT; t<=nT; t++ )
  {
    Mt(1,t) = Mt(1,t-1) + mDevs(1,t);
    Mt(2,t) = Mt(2,t-1) + mDevs(2,t);

    Mta(t)(firstAge,splitAge(2)-1)   = exp( Mt(1,t) );
    Mta(t)(splitAge(2),plusGroupAge) = exp( Mt(2,t) );
  }
  }


FUNCTION calc_sel_gta
  {
  //
  int lt, ut;
  sel_gta.initialize();
  S50_gt.initialize();
  S95_gt.initialize();
  for( int g=1; g<=nFisheries; g++ )
  {
    if( nSelBlocks(g) > 1 )
    {
      // time-varying selectivity in blocks.
      for( int i=1; i<=nSelBlocks(g)-1; i++ )
      {
        lt = tBlock(g,i); ut = tBlock(g,i+1)-1;
        S50_gt(g)(lt,ut) = mfexp( log_S50_gi(g,i) );
        S95_gt(g)(lt,ut) = S50_gt(g)(lt,ut) + mfexp( log_S95_step_gi(g,i) );
      }
      S50_gt(g)(ut+1,nT) = mfexp( log_S50_gi(g,nSelBlocks(g)) );
      S95_gt(g)(ut+1,nT) = S50_gt(g)(ut+1,nT) + mfexp( log_S95_step_gi(g,nSelBlocks(g)) );
    }
    else
    {
      // constant sel: i is vector 1,nT
      S50_gt(g)(1,nT) = mfexp( log_S50_gi(g,1) );
      S95_gt(g)(1,nT) = S50_gt(g) + mfexp( log_S95_step_gi(g,1) );
    }
    // time loop to fill
    for( int t=1; t<=nT; t++ )
    {
      dvariable tmp = log(19.)/( S95_gt(g)(t) - S50_gt(g)(t) );
      sel_gta(g)(t) = 1./( 1. + exp(-tmp*(age - S50_gt(g)(t) ) ) );
    }
  }
  }

FUNCTION popInit
  {
  Bta.initialize(); SSBt.initialize(); N4p.initialize();
  uCgat.initialize(); expBtg.initialize(); B510.initialize();
  SSBa.initialize();

  // Initialize abundance in first year assuming each cohort
  // (except plusGroup) had independent recruitment deviations
  // and age-dependent M.
  // age-1 does not involve M
  Nta(1)(firstAge) = avgR*mfexp( init_recDevs(firstAge) );
  for( int a=firstAge+1; a<plusGroupAge; a++ )
  {
    Nta(1)(a) = avgR*mfexp( init_recDevs(a) - sum( sel_gta(1)(1)(firstAge,a-1)*Finit + Mta(1)(firstAge,a-1) ) );
  }
  Nta(1)(plusGroupAge) = avgR*mfexp( init_recDevs(plusGroupAge) - sum( sel_gta(1)(1)(firstAge,plusGroupAge-1)*Finit +  Mta(1)(firstAge,plusGroupAge-1) ) );
  Nta(1)(plusGroupAge) *= 1./( 1.-mfexp(-( sel_gta(1)(1)(plusGroupAge)*Finit + Mta(1)(plusGroupAge) ) ) );

  // Biomass-at-age and SSB: SSBt(1) should match the value
  // at the end of calcPhi()
  Bta(1)   = elem_prod( Nta(1),column(BwtAge,1) );
  SSBt(1)  = sum( elem_prod( Bta(1),matAge ) );
  N4p(1) = sum(Nta(1)(4,plusGroupAge));
  B510(1) = sum(Bta(1)(5,10));
  //cout << "SSBt(1) = " << SSBt(1) << endl;

  // Solve catch equations, t=1: returns Zta(1) and Ftg(1)(1,nFisheries)
  int t=1;
  Zta.initialize(); Ftg.initialize();
  Zta(1) = sel_gta(1)(t)*Finit + Mta(1);
  //cout << "Zta(1) = " << Zta(1) << endl;
  solveBaranov(t,baranovIter);
  //cout << "Zta(1) = " << Zta(1) << endl;

  // predicted exploitable biomass by index - omitted - no indices start in year 1
  for( int i=1; i<=nIndexSeries; i++ )
  {
    int g = idxIndex(i);
    for( int k=minAge(g); k<=maxAge(g); k++ )
      expBtg(i,1) += Nta(1,k)*IwtAge(i,k,1)*sel_gta(g,1,k)*exp(-fracYearSurvey(i)*Zta(1,k));
  }

  // calculate April 1 SSB
  for( int k=4; k<=11; k++)
      SSBa(1) += Nta(1,k)*IwtAge(1,k,1)*exp(-0.2466*Mta(1,k));

  // predicted age-proportions uCgat in this year's catch
  // based on numbers-at-age
  for( int i=1; i<=nAgeSeries; i++ )
  {
    int g = ageIndex(i);
    dvar_vector N   = Nta(1)(minAge(g),maxAge(g));
    dvar_vector S   = sel_gta(g)(1)(minAge(g),maxAge(g));
    dvar_vector tmp = elem_prod( N, S );
    for( int k=minAge(g); k<=maxAge(g); k++ )
      uCgat(g)(k)(1) = tmp(k)/sum(tmp);
  }

  }

FUNCTION popDynamics
  {
  int t, g;
  for( int t=2; t<=nT; t++ )
  {
    // Age-2 recruitment
    if( recType==0 )
      Nta(t,firstAge) = mfexp( gamma_R*log(Nta(t-1,firstAge)) + (1.-gamma_R)*log_AvgR + recDevs(t) );

    //if( recType==1 )
    //{
    //  rw_recDevs(t) = gamma_R*rw_recDevs(t-1) + recDevs(t);
    //  Nta(t,1)      = mfexp(recRate+rw_recDevs(t))*SSBt(t-1);
    //}

    // age-3 to age-(A-1)
    for( int a=firstAge+1; a<=plusGroupAge-1; a++ )
    {
      Nta(t,a) = Nta(t-1,a-1)*exp( -Zta(t-1,a-1) );
    }
    Nta(t,plusGroupAge) = Nta(t-1,plusGroupAge-1)*exp( -Zta(t-1,plusGroupAge-1) ) +
                            Nta(t-1,plusGroupAge)*exp( -Zta(t-1,plusGroupAge) );

    // Biomass-at-age and SSB
    Bta(t)  = elem_prod( Nta(t),column(BwtAge,t));
    SSBt(t) = sum( elem_prod( Bta(t), matAge ) );
    N4p(t) = sum(Nta(t)(4,plusGroupAge));
    B510(t) = sum(Bta(t)(5,10));
    // Solve catch equations, t: returns Zta(t) and Ftg(t)(1,nFisheries)
    solveBaranov(t,baranovIter);

  // predicted exploitable biomass by index
  for( int i=1; i<=nIndexSeries; i++ )
  {
    int g = idxIndex(i);
    for( int k=minAge(g); k<=maxAge(g); k++ )
      expBtg(i,t) += Nta(t,k)*IwtAge(i,k,t)*sel_gta(g,t,k)*exp(-fracYearSurvey(i)*Zta(t,k));
  }

  // calculate April 1 SSB
  for( int k=4; k<=11; k++)
      SSBa(t) += Nta(t,k)*IwtAge(1,k,t)*exp(-0.2466*Mta(t,k));

    // predicted age-proportions uCgat in this year's catch
    for( int i=1; i<=nAgeSeries; i++ )
    {
      int g = ageIndex(i);
      dvar_vector N   = Nta(t)(minAge(g),maxAge(g));
      dvar_vector S   = sel_gta(g)(t)(minAge(g),maxAge(g));
      dvar_vector tmp = elem_prod( N, S );
      for( int k=minAge(g); k<=maxAge(g); k++ )
        uCgat(g)(k)(t) = tmp(k)/sum(tmp);
    }
  } // end t-loop

  if ( last_phase() ) neval+=1;
  else neval=-2.;  // start neval counter at -2 at start of last phase so it equals admb screen output
  }

FUNCTION calc_rec_prior
  {
  recPrior = 0.;
  if( active(init_recDevs) )
    recPrior =  0.5*norm2(init_recDevs)/square(priorSD_R);
  if( active(recDevs) )
    recPrior += 0.5*norm2(recDevs)/square(priorSD_R);
  // note: these could use separate prior SDs
  }

FUNCTION calc_q_prior
  {
  qPrior = 0.;
  if( active(Qdev))
    qPrior = 0.5*norm2(Qdev)/square(QdevSD);
  }

FUNCTION calc_M_prior
  {
  mPrior = 0.;
  if( active(mDevs) )
  {
    // Normal prior for M random walk residuals
    for( int j=1; j<=nSplitAge; j++ )
      mPrior +=  0.5*norm2(mDevs(j))/square(priorSD_rwM(j));
  }
  // Normal prior for initial M
  for( int j=1; j<=nSplitAge; j++ )
    mPrior += 0.5*square(mfexp(log_M(j)) - prior_initM(j))/square(priorSD_initM(j));

  }


FUNCTION calc_index_likelihood
  {
  // Initialize likelihood function terms.
  int i,t;
  indexLikelihood.initialize();
  predNdx.initialize();
  validObs.initialize();
  ss.initialize();
  dvar_vector z(1,nT);  // residuals

  for ( i=1;i<=nIndexSeries;i++ )
  {
    z.initialize();
    for ( t=1;t<=nT;t++ )
    {
      if ( idxSeries[i][t] > 0. )
      {
        predNdx(i,t) = q(i,t)*expBtg(i,t);
        z[t] = log(idxSeries[i][t]) - log( predNdx(i,t) );
        ss(i) += pow(z[t], 2);
        validObs[i] += 1;
      }
    }
    indexLikelihood(i) = validObs(i)*log_sig(i) + (ss(i)/(2.*square(sigIndex(i))));
  }
  }

FUNCTION calc_age_likelihood
  {
  dvariable meanDiff;
  dvariable nAges;
  dvariable etaSumSq;
  dvariable mnLL;
  int g;
  int astr;
  int aend;
  tauSquareAges.initialize();

  // Calculate predicted age-proportions for each index gear.
  for ( int i=1;i<=nAgeSeries;i++ )
  {
    nAges = 0.;
    etaSumSq   = 0.;
    mnLL       = 0.;
    g          = ageIndex(i);
    for ( int t=ageFirstYear(g); t<=ageLastYear(g); t++ )
    {
      dvar_vector res(1,maxAge(g));
      res.initialize();
      //nYearsAges += 1;
      dvar_vector obsPat  = column(ageObsProp(g),t)(pst(g,t),pnd(g,t));
      astr = minAge(g);
      aend = maxAge(g);
      if(pst(g,t) > minAge(g))
      {
        astr = pst(g,t);
        uCgat(g,astr,t) = sum(column(uCgat(g),t)(minAge(g),astr) );
      }
      if(pnd(g,t) < maxAge(g))
      {
        aend = pnd(g,t);
        uCgat(g,aend,t) = sum(column(uCgat(g),t)(aend,maxAge(g)));
      }
      dvar_vector predPat = ( column(uCgat(g),t)(astr,aend) )/sum( column(uCgat(g),t)(astr,aend) );

      // Column means for p_at residuals.
      int iRes = 0; dvariable sumRes = 0.;
      for( int a=astr; a<=aend; a++ )
      {
        if( obsPat(a) > 0 )
        {
         //cout << obsPat(a) << endl;
          nAges  += 1;
          //cout << "nAges" << nAges <<endl;
          iRes   += 1;
          //cout<< "iRes"<< iRes << endl;
          res(iRes) = log( obsPat(a) ) - log( predPat(a) );
          //cout<<"res"<< res<< endl;
          sumRes += res(iRes);
          //cout << "sumRes" << sumRes <<endl;
        }
        else
        {
          //cout << obsPat(a) << endl;
          mnLL += -50.*log(1.-predPat(a));
        }
      }
      if( iRes > 0 )
      {
        meanDiff = sumRes/iRes;
        // Proportion residual function.
        etaSumSq += norm2(res(1,iRes)-meanDiff);
      }
    }
    // MLE of variance in age-proportions.
    tauSquareAges(g) = etaSumSq/(nAges );
    ageLikelihood(g) = nAges*log(tauSquareAges(g)) + mnLL;
  }
  }

FUNCTION  void solveBaranov( int t, int nIter )
  {
  RETURN_ARRAYS_INCREMENT();
  dvariable f;
  dvariable J;
  dvar_vector Bprime(firstAge,plusGroupAge);
  dvar_vector tmp(firstAge,plusGroupAge);
  dvar_vector Za(firstAge,plusGroupAge);
  dvar_vector ZaNew(firstAge,plusGroupAge);

  // Initialize Z to current vector of Mta
  ZaNew.initialize();
  // Initial approximation of F...
  // Selected biomass
  Bprime = elem_prod(elem_prod( Nta(t),column(CwtAge,t)), sel_gta(1)(t) );
  Ftg(t,1) = Ctg(t,1)/sum( Bprime );
  ZaNew = Mta(t) + sel_gta(1)(t)*Ftg(t,1);
  // refine F for fisheries (surveys not critical)
  for( int i=1; i<=nIter; i++ )
  {
    // Total mortality
    Za=ZaNew; ZaNew=Mta(t);
    // Predicted catch given F
    tmp    = elem_div( elem_prod( Bprime*Ftg(t,1),1.-exp(-Za) ), Za );
    //cout << "predCatch = " << tmp << endl;
    // Function value: difference of pred - obs catch
    f =  sum(tmp) - Ctg(t,1); tmp=0;
    // Jacobian
    dvar_vector tmp1 = elem_div( Bprime, Za );
    dvar_vector tmp2 = elem_prod( sel_gta(1)(t), Za )*Ftg(t,1);
    dvar_vector tmp3 = 1. - mfexp( -Za );
    dvar_vector tmp4 = elem_prod( elem_div( sel_gta(1)(t)*Ftg(t,1), Za ), 1.-exp(-Za)  );

    tmp = elem_prod( tmp1, tmp2 + tmp3 - tmp4 );
    J = sum(tmp); tmp=0;
    Ftg(t,1) -= f/J;
    ZaNew += sel_gta(1)(t)*Ftg(t,1);
    //cout <<"t = " << t<< " iter = "<< i << " f = "<< f <<" J = "<< J << " f/J = " << f/J << endl;
    //cout <<"iter = "<< i << " Ftg = "<< Ftg(t, 1) << endl;
  }
  Zta(t) = ZaNew;
  avgF(t) = sum( elem_prod( Nta(t)(a1F,a2F),sel_gta(1)(t)(a1F,a2F)*Ftg(t,1) ) );
  avgF(t) /= sum(Nta(t)(a1F,a2F));
  //avgF2(t) = sum( elem_prod( Nta(t)(a1F2,a2F2),sel_gta(1)(t)(a1F2,a2F2)*Ftg(t,1) ) );
  //avgF2(t) /= sum(Nta(t)(a1F2,a2F2));
  RETURN_ARRAYS_DECREMENT();
  }

FUNCTION  double getFforward( const int& t, const int& p, const int& nIter, const dvector& Np);
  {
  double fp;
  double Jp;
  double Fest;
  dvector Bp(firstAge,plusGroupAge);
  dvector tmp(firstAge,plusGroupAge);
  dvector Zap(firstAge,plusGroupAge);
  dvector ZapNew(firstAge,plusGroupAge);
  dvector selAge(firstAge,plusGroupAge);

  selAge = value(sel_gta(1)(nT));

  // Initialize Z to current vector of Mta
  ZapNew.initialize();
  // Initial approximation of F...
  // Selected biomass
  Bp = elem_prod(elem_prod( Np,column(CwtAge,nT)), selAge );
  Fest = Cpro(p)/sum( Bp );
  ZapNew = Mpro + selAge*Fest;
  // refine F for fisheries (surveys not critical)
  for( int i=1; i<=nIter; i++ )
  {
    // Total mortality
    Zap=ZapNew; ZapNew=Mpro;
    // Predicted catch given F
    tmp    = elem_div( elem_prod( Bp*Fest,1.-exp(-Zap) ), Zap );
    //cout << "predCatch = " << tmp << endl;
    // Function value: difference of pred - obs catch
    fp =  sum(tmp) - Cpro(p); tmp=0;
    // Jacobian
    dvector tmp1 = elem_div( Bp, Zap );
    dvector tmp2 = elem_prod( selAge, Zap )*Fest;
    dvector tmp3 = 1. - mfexp( -Zap );
    dvector tmp4 = elem_prod( elem_div( selAge*Fest, Zap ), 1.-exp(-Zap)  );

    tmp = elem_prod( tmp1, tmp2 + tmp3 - tmp4 );
    Jp = sum(tmp); tmp=0;
    Fest -= fp/Jp;
    ZapNew += selAge*Fest;
    //cout <<"t = " << t<< " iter = "<< i << " f = "<< f <<" J = "<< J << " f/J = " << f/J << endl;
    //cout <<"iter = "<< i << " Ftg = "<< Ftg(t, 1) << endl;
  }
  return Fest;
  }



FUNCTION projection
 {
   int i; int j; int k; int p;
   dmatrix usewate(firstAge,plusGroupAge,nT+1,pT);
   dmatrix aprwate(firstAge,plusGroupAge,nT+1,pT);

  // randomly pick weight-at-age vectors from last 5 years
   for (j=nT+1; j<=pT; j++)
   {
     k=5.*randu(rng);
     for (i=firstAge; i<=plusGroupAge; i++) usewate(i,j)=BwtAge(i,nT-k);
     for (i=4; i<=9; i++) aprwate(i,j)=IwtAge(1)(i,nT-k);
     for (i=10; i<=11; i++) aprwate(i,j)=IwtAge(1)(9,nT-k);
   }


// ***************************** Natural mortality scenarios ************************


  // assign M to average of last 5 years
  //
  //    and reduce by multiplier
  //
  //     here a 75% reduction (multiply by 0.25 to bring it in to 25% of recent value)

  // Mpro = 0.0;
  // for (i=firstAge; i<=plusGroupAge; i++)
  //  {
  //   for (j=nT-4; j<=nT; j++) Mpro(i) += value(((Mta(j,i))/5)*0.01);
  //  }

  // assign M to average of set of years
  //
  //     here set to early time series when M was low (1978-1982)

  Mpro = 0.0;
  for (i=firstAge; i<=plusGroupAge; i++)
  {
   for (j=nT-41; j<=nT-37; j++) Mpro(i) += value(((Mta(j,i))/5));
  }

  // cout << Mpro << endl;

//*************************************************************************************

  // make data variables for required dvariables
  double log_AR = value(log_AvgR); //
  //double acorR = value(gamma_R); //
  double rnyr;
  dvector Rdev(3,pT);

  dvector MnT(firstAge,plusGroupAge);
  dvector NnT(firstAge,plusGroupAge);
  dvector ssb(1,nT);
  dvector selAge(firstAge,plusGroupAge);
  dvector Rdevuse(nT+1,pT);
  dvector Rdevuse1(nT+1,pT);
    dvector Rdevuse2(nT+1,pT);
      dvector Rdevuse3(nT+1,pT);
        dvector Rdevuse4(nT+1,pT);



  MnT = value(Mta(nT));
  NnT = value(Nta(nT));
  ssb = value(SSBt);
  selAge = value(sel_gta(1)(nT));



// *********************** Recruitment scenarios ************************




/////////////////////// high recruitment scenario

  // get rec residuals for year 2 to 16  of time series
  // (until last big recruitment event, 1994 = year 16)
  // you dont want recdev 0 so you need to put k+2
  //
  //rnyr = nT - firstAge - 1;
  rnyr = 17.0;
  for (j=nT+1; j<=pT; j++)
  {
   k = rnyr*randu(rng);
   //cout<< "k" << k <<endl;
   //cout<< k <<endl;
   Rdevuse(j)=value(recDevs(k+2));
   //cout<<Rdevuse<<endl;
   }

/////////////////////// recent recruitment scenario

  // assign Rdev1 to average of last 5 years
 // for (j=nT+1; j<=pT; j++)
  //{
  //   for (i=nT-4; i<=nT; i++) Rdevuse1(j) += value(recDevs(i))/5;
  //   cout << Rdevuse1(j) << endl;
  //}
  //cout << Rdevuse1(j) << endl;


  // assign recruitment to random recdevs in the last 5 years
  rnyr = 5.;
  //rnyr = 40.0;
  for (j=nT+1; j<=pT; j++)
  {
   k = rnyr*randu(rng);
   //cout<< "k" << k <<endl;
   //cout<< k <<endl;
   Rdevuse1(j)=value(recDevs(nT-k));
   //cout<<k<<endl;
   }

////////////////////// random recruitment scenario

  // get random rec residuals for whole time series

  rnyr = nT - firstAge - 1;
  //rnyr = 40.0;
  for (j=nT+1; j<=pT; j++)
  {
   k = rnyr*randu(rng);
   //cout<< "k" << k <<endl;
   //cout<< k <<endl;
   Rdevuse2(j)=value(recDevs(k+3));
   //cout<<k<<endl;
   }

// ***********************************************************************








// *********************************     FIRST PROJECTION     ************************************************

// this uses the high recruitment years (Rdevuse)

  Fpro = 0.0;
  // PROJECTION 1, F=0
  Sapro1 = 0.0;
  Bpro1 = 0.0;
  Npro1 = 0.0;

  // first projection year

  Npro1(firstAge,nT+1) = mfexp(log_AR+Rdevuse1(nT+1));
  //cout<<Npro1(firstAge,nT+1)<<endl;
  for (i=firstAge+1; i<=plusGroupAge-1; i++) Npro1(i,nT+1) = NnT(i-1)*mfexp(-(MnT(i-1)+selAge(i)*value(Ftg(nT,1))));
  Npro1(plusGroupAge,nT+1) =  NnT(plusGroupAge-1)*mfexp(-(MnT(plusGroupAge-1)+selAge(plusGroupAge-1)*value(Ftg(nT,1))));
  Npro1(plusGroupAge,nT+1) +=  NnT(plusGroupAge)*mfexp(-(MnT(plusGroupAge)+selAge(plusGroupAge)*value(Ftg(nT,1))));
  for (i=4; i<=plusGroupAge; i++) Bpro1(nT+1) += Npro1(i,nT+1)*usewate(i,nT+1);
  for (i=4; i<=plusGroupAge; i++) Sapro1(nT+1) += Npro1(i,nT+1)*aprwate(i,nT+1)*exp(-0.2466*Mpro(i));

    // second projection year

  Npro1(firstAge,nT+2) = mfexp(log_AR+Rdevuse1(nT+2));
  //cout<<Npro1(firstAge,nT+2)<<endl;
  for (i=firstAge+1; i<=plusGroupAge-1; i++) Npro1(i,nT+2) = Npro1(i-1,nT+1)*mfexp(-Mpro(i-1));
  Npro1(plusGroupAge,nT+2) =  Npro1(plusGroupAge-1,nT+1)*mfexp(-Mpro(plusGroupAge-1));
  Npro1(plusGroupAge,nT+2) += Npro1(plusGroupAge,nT+1)*mfexp(-Mpro(plusGroupAge));
  for (i=4; i<=plusGroupAge; i++) Bpro1(nT+2) += Npro1(i,nT+2)*usewate(i,nT+2);
  for (i=4; i<=plusGroupAge; i++) Sapro1(nT+2) += Npro1(i,nT+2)*aprwate(i,nT+2)*exp(-0.2466*Mpro(i));

  // remaining projection years
  for( j=nT+3; j<=pT; j++)
  {

    Npro1(firstAge,j) = mfexp(log_AR+Rdevuse1(j));
    //cout<< Npro1(firstAge,j)<<endl;
    for (i=firstAge+1; i<=plusGroupAge-1; i++) Npro1(i,j) = Npro1(i-1,j-1)*mfexp(-Mpro(i-1));
    Npro1(plusGroupAge,j) =  Npro1(plusGroupAge-1,j-1)*mfexp(-Mpro(plusGroupAge-1));
    Npro1(plusGroupAge,j) += Npro1(plusGroupAge,j-1)*mfexp(-Mpro(plusGroupAge));
    for (i=4; i<=plusGroupAge; i++) Bpro1(j) += Npro1(i,j)*usewate(i,j);
    for (i=4; i<=plusGroupAge; i++) Sapro1(j) += Npro1(i,j)*aprwate(i,j)*exp(-0.2466*Mpro(i));
  }

// *********************************     2nd PROJECTION     ************************************************

// this uses the recent recruitment years (Rdevuse1)

  Sapro2 = 0.0;
  Bpro2 = 0.0;
  Npro2 = 0.0;

  // first projection year

  Npro2(firstAge,nT+1) = mfexp(log_AR+Rdevuse1(nT+1));
  for (i=firstAge+1; i<=plusGroupAge-1; i++) Npro2(i,nT+1) = NnT(i-1)*mfexp(-(MnT(i-1)+selAge(i)*value(Ftg(nT,1))));
  Npro2(plusGroupAge,nT+1) =  NnT(plusGroupAge-1)*mfexp(-(MnT(plusGroupAge-1)+selAge(plusGroupAge-1)*value(Ftg(nT,1))));
  Npro2(plusGroupAge,nT+1) +=  NnT(plusGroupAge)*mfexp(-(MnT(plusGroupAge)+selAge(plusGroupAge)*value(Ftg(nT,1))));
  for (i=4; i<=plusGroupAge; i++) Bpro2(nT+1) += Npro2(i,nT+1)*usewate(i,nT+1);
  for (i=4; i<=plusGroupAge; i++) Sapro2(nT+1) += Npro2(i,nT+1)*aprwate(i,nT+1)*exp(-0.2466*Mpro(i));
  Fpro(2,nT+1)=getFforward( nT+1, 2, baranovIter, column(Npro2,nT+1));


  // second projection year: assume fishery selectivity in years nT+1 - pT is the same as in nT

  Npro2(firstAge,nT+2) = mfexp(log_AR+Rdevuse1(nT+2));
  for (i=firstAge+1; i<=plusGroupAge-1; i++) Npro2(i,nT+2) = Npro2(i-1,nT+1)*mfexp(-(Mpro(i-1) + selAge(i-1)*Fpro(2,nT+1)));
  Npro2(plusGroupAge,nT+2) =  Npro2(plusGroupAge-1,nT+1)*mfexp(-(Mpro(plusGroupAge-1)+ selAge(plusGroupAge-1)*Fpro(2,nT+1)));
  Npro2(plusGroupAge,nT+2) += Npro2(plusGroupAge,nT+1)*mfexp(-(Mpro(plusGroupAge) + selAge(plusGroupAge)*Fpro(2,nT+1)));
  for (i=4; i<=plusGroupAge; i++) Bpro2(nT+2) += Npro2(i,nT+2)*usewate(i,nT+2);
  for (i=4; i<=plusGroupAge; i++) Sapro2(nT+2) += Npro2(i,nT+2)*aprwate(i,nT+2)*exp(-0.2466*Mpro(i));
  Fpro(2,nT+2)=getFforward( nT+2, 2, baranovIter, column(Npro2,nT+2));

  // remaining projection years
  for( j=nT+3; j<=pT; j++)
  {

    Npro2(firstAge,j) = mfexp(log_AR+Rdevuse1(j));
    //cout<< Npro2(firstAge,j)<<endl;
    for (i=firstAge+1; i<=plusGroupAge-1; i++) Npro2(i,j) = Npro2(i-1,j-1)*mfexp(-(Mpro(i-1) + selAge(i-1)*Fpro(2,j-1)));
    Npro2(plusGroupAge,j) =  Npro2(plusGroupAge-1,j-1)*mfexp(-(Mpro(plusGroupAge-1) + selAge(plusGroupAge-1)*Fpro(2,j-1)));
    Npro2(plusGroupAge,j) += Npro2(plusGroupAge,j-1)*mfexp(-(Mpro(plusGroupAge) + selAge(plusGroupAge)*Fpro(2,j-1)));
    for (i=4; i<=plusGroupAge; i++) Bpro2(j) += Npro2(i,j)*usewate(i,j);
    for (i=4; i<=plusGroupAge; i++) Sapro2(j) += Npro2(i,j)*aprwate(i,j)*exp(-0.2466*Mpro(i));
    Fpro(2,j)=getFforward( j, 2, baranovIter, column(Npro2,j));
  }



// *********************************     3rd PROJECTION     ************************************************
// PROJECTION 3

  Sapro3 = 0.0;
  Bpro3 = 0.0;
  Npro3 = 0.0;

  // first projection year

  Npro3(firstAge,nT+1) = mfexp(log_AR+Rdevuse1(nT+1));
  for (i=firstAge+1; i<=plusGroupAge-1; i++) Npro3(i,nT+1) = NnT(i-1)*mfexp(-(MnT(i-1)+selAge(i)*value(Ftg(nT,1))));
  Npro3(plusGroupAge,nT+1) =  NnT(plusGroupAge-1)*mfexp(-(MnT(plusGroupAge-1)+selAge(plusGroupAge-1)*value(Ftg(nT,1))));
  Npro3(plusGroupAge,nT+1) +=  NnT(plusGroupAge)*mfexp(-(MnT(plusGroupAge)+selAge(plusGroupAge)*value(Ftg(nT,1))));
  for (i=4; i<=plusGroupAge; i++) Bpro3(nT+1) += Npro3(i,nT+1)*usewate(i,nT+1);
  for (i=4; i<=plusGroupAge; i++) Sapro3(nT+1) += Npro3(i,nT+1)*aprwate(i,nT+1)*exp(-0.2466*Mpro(i));
  Fpro(3,nT+1)=getFforward( nT+1, 3, baranovIter, column(Npro3,nT+1));


  // second projection year: assume fishery selectivity in years nT+1 - pT is the same as in nT

  Npro3(firstAge,nT+2) = mfexp(log_AR+Rdevuse1(nT+2));
  for (i=firstAge+1; i<=plusGroupAge-1; i++) Npro3(i,nT+2) = Npro3(i-1,nT+1)*mfexp(-(Mpro(i-1) + selAge(i-1)*Fpro(3,nT+1)));
  Npro3(plusGroupAge,nT+2) =  Npro3(plusGroupAge-1,nT+1)*mfexp(-(Mpro(plusGroupAge-1)+ selAge(plusGroupAge-1)*Fpro(3,nT+1)));
  Npro3(plusGroupAge,nT+2) += Npro3(plusGroupAge,nT+1)*mfexp(-(Mpro(plusGroupAge) + selAge(plusGroupAge)*Fpro(3,nT+1)));
  for (i=4; i<=plusGroupAge; i++) Bpro3(nT+2) += Npro3(i,nT+2)*usewate(i,nT+2);
  for (i=4; i<=plusGroupAge; i++) Sapro3(nT+2) += Npro3(i,nT+2)*aprwate(i,nT+2)*exp(-0.2466*Mpro(i));
  Fpro(3,nT+2)=getFforward( nT+2, 3, baranovIter, column(Npro3,nT+2));

  // remaining projection years
  for( j=nT+3; j<=pT; j++)
  {

    Npro3(firstAge,j) = mfexp(log_AR+Rdevuse1(j));
    for (i=firstAge+1; i<=plusGroupAge-1; i++) Npro3(i,j) = Npro3(i-1,j-1)*mfexp(-(Mpro(i-1) + selAge(i-1)*Fpro(3,j-1)));
    Npro3(plusGroupAge,j) =  Npro3(plusGroupAge-1,j-1)*mfexp(-(Mpro(plusGroupAge-1) + selAge(plusGroupAge-1)*Fpro(3,j-1)));
    Npro3(plusGroupAge,j) += Npro3(plusGroupAge,j-1)*mfexp(-(Mpro(plusGroupAge) + selAge(plusGroupAge)*Fpro(3,j-1)));
    for (i=4; i<=plusGroupAge; i++) Bpro3(j) += Npro3(i,j)*usewate(i,j);
    for (i=4; i<=plusGroupAge; i++) Sapro3(j) += Npro3(i,j)*aprwate(i,j)*exp(-0.2466*Mpro(i));
    Fpro(3,j)=getFforward( j, 3, baranovIter, column(Npro3,j));
  }


  // PROJECTION 4

    Sapro4 = 0.0;
    Bpro4 = 0.0;
    Npro4 = 0.0;

    // first projection year

    Npro4(firstAge,nT+1) = mfexp(log_AR+Rdevuse1(nT+1));
    for (i=firstAge+1; i<=plusGroupAge-1; i++) Npro4(i,nT+1) = NnT(i-1)*mfexp(-(MnT(i-1)+selAge(i)*value(Ftg(nT,1))));
    Npro4(plusGroupAge,nT+1) =  NnT(plusGroupAge-1)*mfexp(-(MnT(plusGroupAge-1)+selAge(plusGroupAge-1)*value(Ftg(nT,1))));
    Npro4(plusGroupAge,nT+1) +=  NnT(plusGroupAge)*mfexp(-(MnT(plusGroupAge)+selAge(plusGroupAge)*value(Ftg(nT,1))));
    for (i=4; i<=plusGroupAge; i++) Bpro4(nT+1) += Npro4(i,nT+1)*usewate(i,nT+1);
    for (i=4; i<=plusGroupAge; i++) Sapro4(nT+1) += Npro4(i,nT+1)*aprwate(i,nT+1)*exp(-0.2466*Mpro(i));
    Fpro(4,nT+1)=getFforward( nT+1, 4, baranovIter, column(Npro4,nT+1));


    // second projection year: assume fishery selectivity in years nT+1 - pT is the same as in nT

    Npro4(firstAge,nT+2) = mfexp(log_AR+Rdevuse1(nT+2));
    for (i=firstAge+1; i<=plusGroupAge-1; i++) Npro4(i,nT+2) = Npro4(i-1,nT+1)*mfexp(-(Mpro(i-1) + selAge(i-1)*Fpro(4,nT+1)));
    Npro4(plusGroupAge,nT+2) =  Npro4(plusGroupAge-1,nT+1)*mfexp(-(Mpro(plusGroupAge-1)+ selAge(plusGroupAge-1)*Fpro(4,nT+1)));
    Npro4(plusGroupAge,nT+2) += Npro4(plusGroupAge,nT+1)*mfexp(-(Mpro(plusGroupAge) + selAge(plusGroupAge)*Fpro(4,nT+1)));
    for (i=4; i<=plusGroupAge; i++) Bpro4(nT+2) += Npro4(i,nT+2)*usewate(i,nT+2);
    for (i=4; i<=plusGroupAge; i++) Sapro4(nT+2) += Npro4(i,nT+2)*aprwate(i,nT+2)*exp(-0.2466*Mpro(i));
    Fpro(4,nT+2)=getFforward( nT+2, 4, baranovIter, column(Npro4,nT+2));

    // remaining projection years
    for( j=nT+3; j<=pT; j++)
    {

      Npro4(firstAge,j) = mfexp(log_AR+Rdevuse1(j));
      for (i=firstAge+1; i<=plusGroupAge-1; i++) Npro4(i,j) = Npro4(i-1,j-1)*mfexp(-(Mpro(i-1) + selAge(i-1)*Fpro(4,j-1)));
      Npro4(plusGroupAge,j) =  Npro4(plusGroupAge-1,j-1)*mfexp(-(Mpro(plusGroupAge-1) + selAge(plusGroupAge-1)*Fpro(4,j-1)));
      Npro4(plusGroupAge,j) += Npro4(plusGroupAge,j-1)*mfexp(-(Mpro(plusGroupAge) + selAge(plusGroupAge)*Fpro(4,j-1)));
      for (i=4; i<=plusGroupAge; i++) Bpro4(j) += Npro4(i,j)*usewate(i,j);
      for (i=4; i<=plusGroupAge; i++) Sapro4(j) += Npro4(i,j)*aprwate(i,j)*exp(-0.2466*Mpro(i));
      Fpro(4,j)=getFforward( j, 4, baranovIter, column(Npro4,j));
    }


 }


REPORT_SECTION
  {
  // Parameter estimates
  report << "## Parameter estimates" << endl;
  report << "# avgR "<< endl;        report << avgR << endl;
  report << "# recRate "<< endl;     report << recRate << endl;
  report << "# SSB0" << endl;        report << SSBt(1) << endl;
  report << "# M" << endl;           report << mfexp(log_M) << endl;
  report << "# mDevs" << endl;       report << mDevs << endl;
  report << "# initRecDevs" << endl; report << init_recDevs << endl;
  report << "# recDevs" << endl;     report << recDevs << endl;

  // Selectivity parameters in blocks
  report << "## Selectivity parameters: S50_gi" << endl;
  for( int g=1;g<=nFisheries;g++ )
  {
    for( int j=1; j<=nSelBlocks(g); j++ )
    {
      report << "# S50_"<< g << j << endl;
      report << mfexp(log_S50_gi(g,j)) << endl;
    }
  }
  report << "## Selectivity parameters: S95_gi" << endl;
  for( int g=1;g<=nFisheries;g++ )
  {
    for( int j=1; j<=nSelBlocks(g); j++ )
    {
      report << "# S95_"<< g << j << endl;
      report << mfexp(log_S50_gi(g,j))+mfexp(log_S95_step_gi(g,j)) << endl;
    }
  }

  report << "# Finit" << endl;       report << Finit << endl;
  report << "# gamma_R" << endl;     report << gamma_R << endl;
  report << "# Qdevs" << endl;        report << Qdev << endl;
  report << "# logQ" << endl;        report << log_q << endl;
  report << "# qg" << endl;        report << q << endl;
  report << "# logSig" << endl;        report << log_sig << endl;
  report << "# sigIndex" << endl;        report << sigIndex << endl;

  // Derived likelihood components
  report << "## Standard error estimates" << endl;
  report << "# logSig" << endl;        report << log_sig << endl;
  report << "# sigIndex" << endl;        report << sigIndex << endl;
  report << "# tauAge" << endl;     report << sqrt(tauSquareAges) << endl;

  // Minimization performance
  report << endl;
  report << "## Minimization performance" << endl;
  report << "# indexLikelihood" << endl;  report << indexLikelihood << endl;
  report << "# ageLikelihood" << endl;    report << ageLikelihood << endl;
  report << "# recPrior" << endl;    report << recPrior << endl;
  report << "# qPrior" << endl;    report << qPrior << endl;
  report << "# objFun" << endl;       report << *objective_function_value::pobjfun << endl;
  report << "# maxGrad" << endl;      report << objective_function_value::gmax << endl;
  report << "# exitCode" << endl;     report << iexit << endl;
  report << "# funEvals" << endl;     report << neval << endl;

  // Derived variables
  report << endl;
  report << "## Derived variables" << endl;
  report << "# SSBt" << endl;       report << SSBt(1,nT) << endl;
  report << "# N4p" << endl;       report << N4p << endl;
  report << "# Mt" << endl;         report << exp(Mt) << endl;
  report <<"## Fishing mortality rates"<< endl;
  report << "# avgF"<< endl;        report << avgF << endl;
  report << "# Ftg" << endl;        report << column(Ftg,1)<< endl;

  report << "# Zta" << endl;        report << Zta << endl;
  report << "# Nta" << endl;        report << Nta << endl;
  report << "# Bta" << endl;        report << Bta << endl;
  report << "# B510" << endl;        report << B510 << endl;
  report <<"## Exploitable biomass by fishery (vals < 0 are missing...)"<< endl;
  report << "# expBtg"<< endl;
  for( int i=1;i<=nIndexSeries;i++ )
  {
    report << row(expBtg,i) << endl;
  }

  report <<"## Predicted Indices (vals < 0 are missing...)"<< endl;
  report << "# Itg"<< endl;
  for( int i=1;i<=nIndexSeries;i++ )
  {
    report << row(predNdx,i) << endl;
  }

  report <<"## Predicted age proportions by fishery (if obs series was used)"<< endl;
  for( int i=1;i<=nAgeSeries;i++ )
  {
    report << "# uCgat"<<ageIndex(i)<< endl;
    report << uCgat(i) << endl;
  }

  //--------------------------------------------------------------------------//
  // Echo other inputs
  //--------------------------------------------------------------------------//
  // scaControl.ctl
  report << "# nFisheries" << endl;      report << nFisheries << endl;
  report << "# firstAge" << endl;        report << firstAge << endl;
  report << "# plusGroupAge" << endl;    report << plusGroupAge << endl;
  report << "# nSelBlocks" << endl;      report << nSelBlocks << endl;
  // Time block boundaries for selectivity
  for( int g=1;g<=nFisheries;g++ )
  {
    for( int j=1; j<=nSelBlocks(g); j++ )
    {
      report << "# tBlock_"<< g << j << endl;
      report << tBlock(g,j) << endl;
    }
  }
  report << "# recType" << endl;         report << recType << endl;
  report << "# ph_log_recRate" << endl;  report << ph_log_recRate << endl;
  report << "# baranovIter" << endl;     report << baranovIter << endl;
  report << "# ph_log_AvgR" << endl;  report << ph_log_AvgR << endl;
  report << "# ph_initRecDevs" << endl;  report << ph_initRecDevs << endl;
  report << "# ph_recDevs" << endl;  report << ph_recDevs << endl;
  report << "# ph_M" << endl;  report << ph_M << endl;
  report << "# ph_mDevs" << endl;  report << ph_mDevs << endl;
  report << "# a1F" << endl;  report << a1F << endl;
  report << "# a2F" << endl;  report << a2F << endl;
  report << "# ph_S50" << endl;  report << ph_S50 << endl;
  report << "# ph_S95" << endl;  report << ph_S95 << endl;
  report << "# ph_gammaR" << endl;  report << ph_gammaR << endl;

  // Model parameter priors.
  report << "# priorSD_R" << endl;       report << priorSD_R << endl;
  report << "# prior_initM" << endl;     report << prior_initM << endl;
  report << "# priorSD_initM" << endl;   report << priorSD_initM << endl;
  report << "# priorSD_rwM" << endl;     report << priorSD_rwM << endl;

  // Inputs from scaCatch.dat.
  report << "# nT" << endl;              report << nT << endl;
  report << "# rwT" << endl;              report << rwT << endl;
  report << "# nFisheries" << endl;      report << nFisheries << endl;
  report << "# landCatchMatrix" << endl; report << landCatchMatrix << endl;

  // Inputs from scaIndex.dat.
  report << "# nIndexSeries" << endl;    report << nIndexSeries << endl;
  report << "# idxIndex" << endl;        report << idxIndex << endl;
  report << "# idxRelative" << endl;     report << idxRelative << endl;
  report << "# idxLikeWeight" << endl;   report << idxLikeWeight << endl;
  report << "# idxFirstYear" << endl;    report << idxFirstYear << endl;
  report << "# idxLastYear" << endl;     report << idxLastYear << endl;
  report << "# fracYearSurvey" << endl;  report << fracYearSurvey << endl;
  report << "## idxSeries" << endl;
  for ( int i=1; i<=nIndexSeries;i++ )
  {
    report << "# idxSeries" << idxIndex(i) << endl;
    report << idxSeries(i) << endl;
  }

  // Inputs from scaAges.dat.
  report << "# firstAge" << endl;        report << firstAge << endl;
  report << "# plusGroupAge" << endl;    report << plusGroupAge << endl;
  report << "# nAgeSeries" << endl;      report << nAgeSeries << endl;
  report << "# minAge" << endl;          report << minAge << endl;
  report << "# maxAge" << endl;          report << maxAge << endl;
  report << "# ageIndex" << endl;        report << ageIndex << endl;
  report << "# ageFirstYear" << endl;    report << ageFirstYear << endl;
  report << "# ageLastYear" << endl;     report << ageLastYear << endl;
  report << "## ageObsProp (written as nAgeSeries slices)" << endl;
  for ( int i=1; i<=nAgeSeries;i++ )
  {
    report << "# paaSeries" << ageIndex(i) << endl;
    report << ageObsProp(i) << endl;
  }
  report << "# init_tauAge " << endl;    report << init_tauAge << endl;

  report << "# ages " << endl;          report << age << endl;
  report << "# BwtAge" << endl;            report << BwtAge << endl;
  report << "# CwtAge" << endl;            report << CwtAge << endl;
  report << "# cpueWAA" << endl;            report << IwtAge(1) << endl;
  report << "# acouWAA" << endl;            report << IwtAge(2) << endl;
  report << "# matAge" << endl;             report << matAge << endl;
  report <<"## Selectivity function by fishery"<< endl;
  for( int g=1;g<=nFisheries;g++ )
  {
    report << "# sel_"<<g<< endl; report << sel_gta(g) << endl;
  }

  // projections parameters
  report << "# Mpro" << endl;             report << Mpro << endl;

  //--------------------------------------------------------------------------//
  // End echo of inputs                                                       //
  //--------------------------------------------------------------------------//
  }
