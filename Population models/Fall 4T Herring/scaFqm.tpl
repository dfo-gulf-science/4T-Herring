// Statistical catch-at-age model for 4T Fall spawning Atlantic Herring
// Acoustic index is in 1000t, expect q < 0.001
// model scale 1000 fish. input WAA in kg

DATA_SECTION
  init_int dum;
  // Control parameters -- switch to estimation control file
  !! ad_comm::change_datafile_name("scaControl.ctl");
  // Total number of years for the model.
  init_int nT;
  // Number of populations for the model.
  init_int nP;
  // Number of fisheries
  init_int nFisheries;
  // year rw M starts
  init_int rwT;
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
  // ph_M
  init_int ph_M;
  // ph_Mdevs
  init_int ph_mDevs;
  // prior_initM
  init_vector prior_initM(1,nSplitAge);
  // prior_SDinitM
  init_vector priorSD_initM(1,nSplitAge);  // std devs
  // priorSD_reM
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
  // assumed to be the same for all popns
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
  !! ad_comm::change_datafile_name("scaCatch2.dat");
  // Matrix of landed catch biomass (tonnes) values.
  init_matrix landCatchMatrix(1,nT,1,nP);
  vector propCatchFish(1,nP);
  // Catch by POPN
  matrix  Ctg(1,nT,1,nP);
  init_int eof2;

 LOC_CALCS
   if(eof2!=2222)
   {
     cout<<" data entry error catch file " <<endl;
     exit(1);
   }
 END_CALCS

  // Abundance index data -- switch to index data file
  //
  !! ad_comm::change_datafile_name("scaIndex.dat");
  // Number of 3P indices
  init_int nIndex3P;
  // Number of T indices
  init_int nIndexT;
  // Total number of indices
  init_int nIndex;
  //  Series index numbers.
  init_ivector idxIndex(1,nIndex);
  // Is index relative? (0=false, 1=true).
  init_vector idxRelative(1,nIndex);
  // Index weights in overall likelihood
  init_vector idxLikeWeight(1,nIndex);
  // Vector indicating first year for each index.
  init_ivector idxFirstYear(1,nIndex);
  // Vector indicating last year for each index.
  init_ivector idxLastYear(1,nIndex);
  // Fraction of year when survey occurs
  init_vector fracYearSurvey(1,nIndex);
  // Fraction of year before fishery Aug 1
  init_number fracAug1;
  // 3P Index series data
  init_3darray idxSeries3P(1,nIndex3P,1,nP,1,nT);
  // T Index series data
  init_matrix idxSeriesT(1,nIndexT,1,nT);
  // Create vector of max 3P index series lengths.
  ivector nobs3P(1,nIndex3P);
  !!nobs3P(1)=idxLastYear(1)-idxFirstYear(1)+1;
  !!nobs3P(2)=idxLastYear(2)-idxFirstYear(2)+1;
  matrix validObs3P(1,nIndex3P,1,nP);
  // Create vector of max T index series lengths.
  ivector nobsT(3,4);
  !!nobsT=idxLastYear(3,4)-idxFirstYear(3,4)+1;
  vector validObsT(1,nIndexT);

  init_int eof3;

 LOC_CALCS
   if(eof3!=3333)
   {
     cout<<" data entry error index file " <<endl;
     cout<<" data entry error eof3"<< eof3 <<endl;
     cout<<" data entry error fracYearSurvey "<< fracYearSurvey <<endl;
     cout <<"cpue(p=1,t=9) " << idxSeries3P(1,1,9) <<endl;
     cout <<"cpue(p=3,t=45) " << idxSeries3P(1,3,45) <<endl;
     cout <<"xnet(p=1,t=25) " << idxSeries3P(2,1,25) <<endl;
     cout <<"xnet(p=3,t=45) " << idxSeries3P(2,3,45) <<endl;
     cout <<"RV(t=17) " << idxSeriesT(1,17) <<endl;
     cout <<"RV(t=45) " << idxSeriesT(1,45) <<endl;
     cout <<"ACR(t=17) " << idxSeriesT(2,17) <<endl;
     cout <<"ACR(t=44) " << idxSeriesT(1,44) <<endl;
     exit(1);
   }
 END_CALCS


  // Age-related data -- prop-age, weight-age, mat-age
  !! ad_comm::change_datafile_name("scaAges.dat");
  // Number of 3P age series
  init_int n3PAgeSeries;
  // Number of T age series
  init_int nTAgeSeries;
  // Number of age series
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
  //Observed proportions-at-age Fishery
  init_3darray agePropFish(1,nP,firstAge,plusGroupAge,1,nT);
  //Observed proportions-at-age cpue
  init_3darray agePropCpue(1,nP,firstAge,plusGroupAge,1,nT);
  //Observed proportions-at-age xnet
  init_3darray agePropXnet(1,nP,firstAge,plusGroupAge,1,nT);
  //Observed proportions-at-age rv
  init_matrix agePropRV(4,6,17,45);
  //Observed proportions-at-age acorec
  init_matrix agePropACR(2,3,17,44);
  // Beginning-of-year Weight-at-age by year
  init_3darray BwtAge(1,nP,firstAge,plusGroupAge,1,nT);
  // Commercial fishery Weight-at-age by year
  init_3darray CwtAge(1,nP,firstAge,plusGroupAge,1,nT);
  // cpue Weight-at-age by year
  init_3darray FwtAge(1,nP,firstAge,plusGroupAge,1,nT);
  // xnet Weight-at-age by year
  init_3darray XwtAge(1,nP,firstAge,plusGroupAge,1,nT);
  // rv Weight-at-age by year
  init_matrix RwtAge(4,6,17,45);
  // ACR Weight-at-age by year
  init_matrix AwtAge(2,3,17,44);
  // Maturity-at-age
  init_vector matAge(firstAge,plusGroupAge);
  // Initial values for 3P age observation std errors
  init_matrix init_tauAge3P(1,n3PAgeSeries,1,nP);
  // Initial values for T age observation std errors
  init_vector init_tauAgeT(1,nTAgeSeries);
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
     cout<<" data entry error init_tauAge(1) "<< init_tauAge3P <<endl;
     cout<<" data entry error BwtAge 1 2 1_5  "<< BwtAge <<endl;
     cout<<" data entry error CwtAge 3 11 37  "<< CwtAge <<endl;
     cout<<" data entry error XwtAge 2 8 25 "<< XwtAge <<endl;
     cout<<" data entry error XwtAge 3 9 37 "<< AwtAge <<endl;
     cout<<" data entry error e0f4 "<< eof4 <<endl;
     exit(1);
   }
 END_CALCS

     // selectivity at age data -- switch to selectivity data file
    !! ad_comm::change_datafile_name("scaSel.dat");
    // Matrix of proportion 258 gillnets by popn
    init_matrix p258(1,nP,9,nT);
    // cpue calibration block
    init_int y1cal; init_int yTcal;
    init_int a1cal; init_int aTcal;
    // relative selectivity 258 mesh (1986, 2017,4,10)
    init_matrix rsel1(y1cal,yTcal,a1cal,aTcal);
    // relative selectivity 234 mesh (1986, 2017,4,10)
    init_matrix rsel2(y1cal,yTcal,a1cal,aTcal);
    // xnets calibration block
    init_int y1calx; init_int yTcalx;
    init_int a1calx; init_int aTcalx;
    // relative selectivity xnets (2002, 2017,3,9)
    init_matrix rsel3(y1calx,yTcalx,a1calx,aTcalx);

    // eof5
    init_int eof5;

 LOC_CALCS
    if(eof5!=5432)
    {
      cout<<" data entry error selectivity file " <<endl;
      cout<<" data entry error eof5 " << eof5  <<endl;
      cout<<" data entry error y1cal " << y1cal  <<endl;
      cout<<" data entry error y1calx " << y1calx  <<endl;
      cout<<" data entry error rsel1(y1cal,aTcal) " << rsel1  <<endl;
      cout<<" data entry error rsel2(y1cal,aTcal) " << rsel2  <<endl;
      cout<<" data entry error rsel3(y1cal,aTcal) " << rsel3  <<endl;
      exit(1);
    }
    //cout << " nP= " << nP <<  endl;
    cout<<" data read " <<endl;
 END_CALCS

//*******************************************************************/
PARAMETER_SECTION
  objective_function_value objFunc;
  // Average recruitment: initialized at 16 based on gbcod1fx2.par
  init_vector log_AvgR(1,nP,ph_log_AvgR);
  // Recruitment deviations: initial abundance 1978
  init_bounded_matrix init_recDevs(1,nP,firstAge,plusGroupAge,-5.,5.,ph_initRecDevs);
  // Recruitment deviations: 1978-2011
  init_bounded_matrix recDevs(1,nP,2,nT,-5.,5.,ph_recDevs);
  // fishery Selectivity parameters
  init_bounded_matrix log_S50_F(1,nP,1,3,0.6,2.0,ph_S50);
  init_bounded_matrix log_S95_step_F(1,nP,1,3,0.0,2.0,ph_S95);
  // CPUE Selectivity parameters
  init_bounded_vector log_S50_C(1,nP,0.6,2.0,ph_S50);
  init_bounded_vector log_S95_step_C(1,nP,0.0,2.0,ph_S95);
  // xnet Selectivity parameters
  init_bounded_vector log_S50_X(1,nP,0.6,2.0,ph_S50);
  init_bounded_vector log_S95_step_X(1,nP,0.0,2.0,ph_S95);
  // rv Selectivity parameters
  init_bounded_number log_S50_R(0.6,2.0,ph_S50);
  init_bounded_number log_S95_step_R(0.0,2.0,ph_S95);
  // ac Selectivity parameters
  init_bounded_number log_S50_A(0.6,2.0,ph_S50);
  init_bounded_number log_S95_step_A(0.0,2.0,ph_S95);
  // Initial fishing mortality prior to 1978
  init_vector Finit(1,nP);
  // Recruitment autocorrelation
  init_number logit_gamma_R(ph_gammaR);
  // q
  init_bounded_matrix Qdev(1,nP,QdevFirst,QdevLast,-5,5,phzQdev);
  init_bounded_matrix log_q3P(1,nIndex3P,1,nP,-12,1.5,1);
  init_bounded_matrix log_sig3P(1,nIndex3P,1,nP,-3.0,1.5,4);
  init_bounded_vector log_qT(1,nIndexT,-12,1.5,1);
  init_bounded_vector log_sigT(1,nIndexT,-3.0,1.5,4);
  // Initial natural mortality rates
  init_matrix log_M(1,nP,1,nSplitAge,ph_M);
  // Natural mortality ranWalk deviations
  init_bounded_matrix mDevN(1,nSplitAge,rwT,nT,-2.,2.,ph_mDevs);
  init_bounded_matrix mDevM(1,nSplitAge,rwT,nT,-2.,2.,ph_mDevs);
  init_bounded_matrix mDevS(1,nSplitAge,rwT,nT,-2.,2.,ph_mDevs);


  // Parameters: natural scale
  vector avgR(1,nP);
  number gamma_R;
  // Selectivity parameters: age-at-50%
  matrix S50_ptF(1,nP,1,nT);
  matrix S50_ptC(1,nP,1,nT);
  matrix S50_ptX(1,nP,1,nT);
  vector S50_tR(1,nT);
  vector S50_tA(1,nT);
  // Selectivity parameters: age-at-95%
  matrix S95_ptF(1,nP,1,nT);
  matrix S95_ptC(1,nP,1,nT);
  matrix S95_ptX(1,nP,1,nT);
  vector S95_tR(1,nT);
  vector S95_tA(1,nT);
  // natural mortality
  matrix M_init(1,nP,1,nSplitAge);
  3darray Mt(1,nP,1,nSplitAge,1,nT);
  3darray Mta(1,nP,1,nT,firstAge,plusGroupAge);
  3darray Mdv(1,nP,1,nSplitAge,rwT,nT);

  // Posterior quantities
  number recPrior;
  number qPrior;
  number mPrior;
  vector indexLikelihood(1,nIndex);
  matrix subNdxLike(1,nIndex3P,1,nP);
  vector ageLikelihood(1,nAgeSeries);
  matrix subAgeLike(1,n3PAgeSeries,1,nP);
  matrix tauSquareAges3P(1,n3PAgeSeries,1,nP);
  vector tauSquareAgesT(1,nTAgeSeries);
  matrix validObs3P(1,nIndex3P,1,nP); // in both data and parameter sections in spr sca also
  vector validObsT(1,nIndexT); // in both data and parameter sections in spr sca also
  3darray q3P(1,nIndex3P,1,nP,1,nT);
  vector qT(1,nIndexT);

  matrix ss3P(1,nIndex3P,1,nP);
  vector ssT(1,nIndexT);
  matrix sigIndex3P(1,nIndex3P,1,nP);
  vector sigIndexT(1,nIndexT);

  // Derived variables
  // Unfished SSB-per-recruit
  vector phiSSB(1,nP);
  matrix rw_recDevs(1,nP,1,nT);
  // change recRate to a number parameter if recType is 1
  matrix recRate2(1,nP,3,nT);
  matrix recRate4(1,nP,5,nT);
  matrix SSBt(1,nP,1,nT);
  matrix SSBa(1,nP,1,nT);
  3darray predNdx3P(1,nIndex3P,1,nP,1,nT);
  matrix predNdxT(1,nIndexT,1,nT);
  matrix N4p(1,nP,1,nT);
  matrix N4(1,nP,1,nT);
  matrix N2(1,nP,1,nT);
  matrix B510(1,nP,1,nT);
  matrix N510(1,nP,1,nT);
  3darray Nta(1,nP,1,nT,firstAge,plusGroupAge);
  3darray Bta(1,nP,1,nT,firstAge,plusGroupAge);
  3darray Zta(1,nP,1,nT,firstAge,plusGroupAge);
  matrix Ftg(1,nT,1,nP);
  3darray expBtg3P(1,nIndex3P,1,nP,1,nT);
  matrix expBtgT(1,nIndexT,1,nT);

  sdreport_matrix avgF(1,nP,1,nT);

  // Bprime deleted - not used outside baranov function
  3darray uCgatF(1,nP,firstAge,plusGroupAge,1,nT);
  3darray uCgatC(1,nP,firstAge,plusGroupAge,1,nT);
  3darray uCgatX(1,nP,firstAge,plusGroupAge,1,nT);
  matrix uCgatR(4,6,17,45);
  matrix uCgatA(2,3,17,44);

  // Selectivity by age
  //matrix sel(1,nFisheries,1,plusGroupAge);
  3darray selF_ta(1,nP,1,nT,firstAge,plusGroupAge);
  3darray selC_ta(1,nP,1,nT,firstAge,plusGroupAge);
  3darray selX_ta(1,nP,1,nT,firstAge,plusGroupAge);
  matrix selR_ta(1,nT,firstAge,plusGroupAge);
  matrix selA_ta(1,nT,firstAge,plusGroupAge);

  // total population
  vector TSSB(1,nT);
  vector TaSSB(1,nT);
  vector TN4p(1,nT);
  vector TN4(1,nT);
  vector TN2(1,nT);


 //*******************************************************************/
GLOBALS_SECTION
  void solveBaranov(const int& p, const int& t, const int& iter);

  // Flag to control whether header written to mcmc output file.
  int mcmcHeader = 0;

  // Flag to control which parameters are done in MCMC phase.
  int mcmcFlag = 1;

  // Uncomment this line for 32-bit compiler.
  // #include <fstream.h>
  #include <admodel.h>
  ofstream mcoutSSB1("mcoutSSB1.dat");
  ofstream mcoutSSB2("mcoutSSB2.dat");
  ofstream mcoutSSB3("mcoutSSB3.dat");
  ofstream mcoutSSBT("mcoutSSBT.dat");
  ofstream mcoutSSBa1("mcoutSSBa1.dat");
  ofstream mcoutSSBa2("mcoutSSBa2.dat");
  ofstream mcoutSSBa3("mcoutSSBa3.dat");
  ofstream mcoutSSBaT("mcoutSSBaT.dat");
  ofstream mcoutF1("mcoutF1.dat");
  ofstream mcoutF2("mcoutF2.dat");
  ofstream mcoutF3("mcoutF3.dat");
  ofstream mcoutRec2N("mcoutRec2N.dat");
  ofstream mcoutRec2M("mcoutRec2M.dat");
  ofstream mcoutRec2S("mcoutRec2S.dat");
  ofstream mcoutRec4N("mcoutRec4N.dat");
  ofstream mcoutRec4M("mcoutRec4M.dat");
  ofstream mcoutRec4S("mcoutRec4S.dat");
  ofstream mcoutRecR2N("mcoutRecR2N.dat");
  ofstream mcoutRecR2M("mcoutRecR2M.dat");
  ofstream mcoutRecR2S("mcoutRecR2S.dat");
  ofstream mcoutRecR4N("mcoutRecR4N.dat");
  ofstream mcoutRecR4M("mcoutRecR4M.dat");
  ofstream mcoutRecR4S("mcoutRecR4S.dat");
  ofstream mcoutN4p1("mcoutN4p1.dat");
  ofstream mcoutN4p2("mcoutN4p2.dat");
  ofstream mcoutN4p3("mcoutN4p3.dat");
  ofstream mcoutN4pT("mcoutN4pT.dat");
  ofstream mcoutPredC1("mcoutPredC1.dat");
  ofstream mcoutPredC2("mcoutPredC2.dat");
  ofstream mcoutPredC3("mcoutPredC3.dat");
  ofstream mcoutPredX1("mcoutPredX1.dat");
  ofstream mcoutPredX2("mcoutPredX2.dat");
  ofstream mcoutPredX3("mcoutPredX3.dat");
  ofstream mcoutPredR("mcoutPredR.dat");
  ofstream mcoutPredA("mcoutPredA.dat");
  ofstream mcoutqC1("mcoutq11.dat");
  ofstream mcoutqC2("mcoutq12.dat");
  ofstream mcoutqC3("mcoutq13.dat");
  ofstream mcoutqX1("mcoutq21.dat");
  ofstream mcoutqX2("mcoutq22.dat");
  ofstream mcoutqX3("mcoutq23.dat");
  ofstream mcoutqR("mcoutqR.dat");
  ofstream mcoutqA("mcoutqA.dat");
  ofstream mcoutMn1("mcoutMn1.dat");
  ofstream mcoutMn2("mcoutMn2.dat");
  ofstream mcoutMm1("mcoutMm1.dat");
  ofstream mcoutMm2("mcoutMm2.dat");
  ofstream mcoutMs1("mcoutMs1.dat");
  ofstream mcoutMs2("mcoutMs2.dat");

TOP_OF_MAIN_SECTION
  arrmblsize=20000000;
  gradient_structure::set_CMPDIF_BUFFER_SIZE(25000000);
  gradient_structure::set_GRADSTACK_BUFFER_SIZE(1000000);

PRELIMINARY_CALCS_SECTION
    for ( int t=1; t<=nT; t++ )
    {
      for ( int p=1; p<=nP; p++ )
      {
        //cout << "prelim calca t= p= "<< t << p  <<endl;
        Ctg(t,p) = landCatchMatrix(t,p);
      }
    }
  neval = 0;
  for (int i=1; i<=3; i++)
  {
    for (int p=1; p<=nP; p++)
    {
      tauSquareAges3P(i,p)  = square(init_tauAge3P(i,p));
    }
  }
  for (int i=1; i<=2; i++)
  {
      tauSquareAgesT(i)  = square(init_tauAgeT(i));
  }

  //cout << "prelim calc done" <<endl;


PROCEDURE_SECTION
  initModelParameters();
  //cout << "initModelParameters done..." << endl;
  popInit();
  //cout << "popInit done..." << endl;
  popDynamics();
  //cout << "popDynamics done..." << endl;
  getTotals();
  //cout << "getTotals done..." << endl;
  calc_index_likelihood();
  //cout << "index_likelihood done..." << endl;
  calc_age_likelihood();
  //cout << "age_likelihood done..." << endl;
  calc_rec_prior();
  //cout << "rec_prior done..." << endl;
  calc_q_prior();
  calc_M_prior();
  objFunc  = 0.;
  objFunc  = sum( elem_prod(idxLikeWeight, indexLikelihood) ) + sum( ageLikelihood );
  objFunc += recPrior;
  objFunc += qPrior;
  objFunc += mPrior;

  if ( mceval_phase() )
  {
    mcoutSSB1 << SSBt(1)    << endl;
    mcoutSSB2 << SSBt(2)    << endl;
    mcoutSSB3 << SSBt(3)    << endl;
    mcoutSSBT << TSSB    << endl;
	  mcoutSSBa1 << SSBa(1)    << endl;
    mcoutSSBa2 << SSBa(2)    << endl;
    mcoutSSBa3 << SSBa(3)    << endl;
    mcoutSSBaT << TaSSB    << endl;
    mcoutF1   << avgF(1)   << endl;
    mcoutF2   << avgF(2)   << endl;
    mcoutF3   << avgF(3)   << endl;
    mcoutRec2N << column(Nta(1),2)    << endl;
    mcoutRec2M << column(Nta(2),2)    << endl;
    mcoutRec2S << column(Nta(3),2)    << endl;
    mcoutRec4N << column(Nta(1),4)    << endl;
    mcoutRec4M << column(Nta(2),4)    << endl;
    mcoutRec4S << column(Nta(3),4)    << endl;
    mcoutRecR2N << recRate2(1)   << endl;
    mcoutRecR2M << recRate2(2)   << endl;
    mcoutRecR2S << recRate2(3)   << endl;
    mcoutRecR4N << recRate4(1)   << endl;
    mcoutRecR4M << recRate4(2)   << endl;
    mcoutRecR4S << recRate4(3)   << endl;
    mcoutN4p1 << N4p(1)    << endl;
    mcoutN4p2 << N4p(2)    << endl;
    mcoutN4p3 << N4p(3)    << endl;
    mcoutN4pT << TN4p    << endl;
    mcoutPredC1 << predNdx3P(1)(1)(9,45)    << endl;
    mcoutPredC2 << predNdx3P(1)(2)(9,45)    << endl;
    mcoutPredC3 << predNdx3P(1)(3)(9,45)    << endl;
    mcoutPredX1 << predNdx3P(2)(1)(5,45)    << endl;
    mcoutPredX2 << predNdx3P(2)(2)(5,45)    << endl;
    mcoutPredX3 << predNdx3P(2)(3)(5,45)    << endl;
    mcoutPredR << predNdxT(1)(17,45)    << endl;
    mcoutPredA << predNdxT(2)(17,44)    << endl;
    mcoutqC1 << q3P(1)(1)(9,45) << endl;
    mcoutqC2 << q3P(1)(2)(9,45) << endl;
    mcoutqC3 << q3P(1)(3)(9,45) << endl;
    mcoutqX1 << q3P(2)(1)(45) << endl;
    mcoutqX2 << q3P(2)(2)(45) << endl;
    mcoutqX3 << q3P(2)(3)(45) << endl;
    mcoutqR << qT(1) << endl;
    mcoutqA << qT(2) << endl;
    mcoutMn1 << exp(Mt(1)(1)) << endl;
    mcoutMn2 << exp(Mt(1)(2)) << endl;
    mcoutMm1 << exp(Mt(2)(1)) << endl;
    mcoutMm2 << exp(Mt(2)(2)) << endl;
    mcoutMs1 << exp(Mt(3)(1)) << endl;
    mcoutMs2 << exp(Mt(3)(2)) << endl;
    mcmcHeader = 1;
  }

FUNCTION initModelParameters
  {
  avgR    = mfexp( log_AvgR );
  qT  =  mfexp(log_qT);
  M_init = mfexp(log_M);
  sigIndex3P = mfexp(log_sig3P);
  sigIndexT = mfexp(log_sigT);
  gamma_R = exp( logit_gamma_R )/(1.+exp(logit_gamma_R));
  //cout << "avgR = " << avgR << endl;
  Mta.initialize();
  Mt.initialize();
  Mdv(1) = mDevN;
  Mdv(2) = mDevM;
  Mdv(3) = mDevS;
  calc_Mta_1Split();
  calc_sel_gta();
  //cout << "calc_sel_gta done"  <<endl;
  get_q();
  calcPhi();
  //cout << "calcPhi"  <<endl;
  }

FUNCTION get_q
    int i, t, p;
    for (p=1; p<=3; p++)
    {
      q3P(1)(p) = mfexp(log_q3P(1,p));
      q3P(2)(p) = mfexp(log_q3P(2,p));
      if(active(Qdev))
      {
        for (t=QdevFirst; t<=QdevLast; t++) q3P(1,p,t)=q3P(1,p,t-1)*mfexp(Qdev(p,t));
      }
    }

FUNCTION calcPhi
  {
  // SPC: calculating unfished SSB-per-recruit. Borrowing
  // Nta and Bta for now. Initializing with 1 recruit
  // to get through calcs.
  Nta.initialize(); SSBt.initialize();
  for (int p=1; p<=3; p++)
  {
    Nta(p)(1,firstAge) = 1.;
    for( int a=firstAge+1; a<=(plusGroupAge-1); a++ )
    {
      Nta(p)(1,a) = Nta(p)(1,a-1)*exp( -M);
    }
    Nta(p)(1,plusGroupAge) = Nta(p)(1,plusGroupAge-1)*exp(-M)/(1.-exp(-M));
    Bta(p)(1) = elem_prod( Nta(p)(1), column(BwtAge(p),1) ); // units here are units of Wat
    phiSSB(p) = sum( elem_prod( Bta(p)(1),matAge ) );
  // Unfished spawning biomass.
  SSBt(p,1) = avgR(p)*phiSSB(p);
  //cout << "phiSSB  = " << phiSSB << endl;
  //cout << "SSBt(1) = " << SSBt(1) << endl;
  //cout << "calcPhi done..." << endl;
  }
  }

FUNCTION calc_Mta_1Split
  {
  // Random walk for 2 age blocks for M
  // Build an Mta matrix that inclues age and time effects.
  // Set initial Mt's in years 1:9. All on log scale until the end.
   for( int p=1; p<=3; p++ )
    {
    for( int t=1; t<rwT; t++ )
     {
      Mt(p)(1)(t)      = log_M(p,1);
      Mt(p)(2)(t)      = log_M(p,2);
      Mta(p)(t)(firstAge,splitAge(2)-1)   = exp( Mt(p,1,t) );
      Mta(p)(t)(splitAge(2),plusGroupAge) = exp( Mt(p,2,t) );
     }
    // Fill random walk.
    for( int t=rwT; t<=nT; t++ )
     {
      Mt(p,1,t) = Mt(p,1,t-1) + Mdv(p,1,t);
      Mt(p,2,t) = Mt(p,2,t-1) + Mdv(p,2,t);

      Mta(p)(t)(firstAge,splitAge(2)-1)   = exp( Mt(p,1,t) );
      Mta(p)(t)(splitAge(2),plusGroupAge) = exp( Mt(p,2,t) );
     }
    }
  }

FUNCTION calc_sel_gta
  {
  //
  int lt, ut, g;
  dvariable tmp;
  selF_ta.initialize();
  selC_ta.initialize();
  selX_ta.initialize();
  selR_ta.initialize();
  selA_ta.initialize();
  S50_ptF.initialize();
  S50_ptC.initialize();
  S50_ptX.initialize();
  S50_tR.initialize();
  S50_tA.initialize();
  S95_ptF.initialize();
  S95_ptC.initialize();
  S95_ptX.initialize();
  S95_tR.initialize();
  S95_tA.initialize();

  for (int p=1; p<=3;  p++)
   {
    // commercial fishery selectivity
    // time-varying selectivity in blocks.
      g=1;
      for( int i=1; i<=nSelBlocks(g)-1; i++ )
      {
        lt = tBlock(g,i); ut = tBlock(g,i+1)-1;
        S50_ptF(p)(lt,ut) = mfexp( log_S50_F(p,i) );
        S95_ptF(p)(lt,ut) = S50_ptF(p)(lt,ut) + mfexp( log_S95_step_F(p,i) );
      }
      S50_ptF(p)(ut+1,nT) = mfexp( log_S50_F(p,nSelBlocks(g)) );
      S95_ptF(p)(ut+1,nT) = S50_ptF(p)(ut+1,nT) + mfexp( log_S95_step_F(p,nSelBlocks(g)) );
      //cout << "log_S50_F " << " p " << log_S50_F(p,1) <<endl;
      //cout << "log_S95_step_F " << " p "  << log_S95_step_F(p,1) <<endl;

      // time loop to fill
      for( int t=1; t<=nT; t++ )
      {
        tmp = log(19.)/( S95_ptF(p)(t) - S50_ptF(p)(t) );
        selF_ta(p)(t) = 1./( 1. + exp(-tmp*(age - S50_ptF(p)(t) ) ) );
      }
      //cout << "done selF_ta" <<endl;

    // cpue selectivity based on input rsel and average logistic selectivity
    // first calculate unadjusted selC_ta, 1 selection block
          S50_ptC(p)(1,nT) = mfexp( log_S50_C(p) );
          S95_ptC(p)(1,nT) = S50_ptC(p) + mfexp( log_S95_step_C(p) );
         //cout << "log_S50_C " << " p " << " i " << log_S50_C(p) <<endl;
         //cout << "log_S95_step_C " << " p " << " i " << log_S95_step_C(p) <<endl;

          // time loop to fill
          for( int t=1; t<=nT; t++ )
          {
            tmp = log(19.)/( S95_ptC(p)(t) - S50_ptC(p)(t) );
            selC_ta(p)(t) = 1./( 1. + exp(-tmp*(age - S50_ptC(p)(t) ) ) );
          }
          //cout << "done selC_ta" <<endl;

          // now adjust for changes in length at age
          for( int t=9; t<=nT; t++)
          {
            for ( int a=4; a<=10; a++)
            {
             selC_ta(p,t,a) *= (p258(p,t)*rsel1(t,a) + (1-p258(p,t))*rsel2(t,a));
            }
          }
          //cout << "done adjusting selC_ta" <<endl;
    // xnet selectivity based on input rsel3 and average logistic selectivity
    // first calculate unadjusted selX_ta, 1 selection block
      S50_ptX(p)(1,nT) = mfexp( log_S50_X(p) );
      S95_ptX(p)(1,nT) = S50_ptX(p) + mfexp( log_S95_step_X(p) );
      //cout << "log_S50_X " << " p " << log_S50_X(p) <<endl;
      //cout << "log_S95_step_X " << " p "  << log_S95_step_X(p) <<endl;

      // time loop to fill
      for( int t=1; t<=nT; t++ )
      {
        dvariable tmp = log(19.)/( S95_ptX(p)(t) - S50_ptX(p)(t) );
        selX_ta(p)(t) = 1./( 1. + exp(-tmp*(age - S50_ptX(p)(t) ) ) );
      }
      //cout << "done selX_ta" <<endl;

      // now adjust for changes in length at age
      for( int t=38; t<=nT; t++)
      {
        for ( int a=3; a<=9; a++)
        {
          selX_ta(p,t,a) *= rsel3(t,a);
        }
      }
      //cout << "done adjusting selX_ta" <<endl;
    } // end P loop
      // calculate RV selectivity
      S50_tR = mfexp( log_S50_R );
      S95_tR = S50_tR + mfexp( log_S95_step_R );
      // time loop to fill
      for( int t=1; t<=nT; t++ )
      {
        dvariable tmp = log(19.)/( S95_tR(t) - S50_tR(t) );
        selR_ta(t) = 1./( 1. + exp(-tmp*(age - S50_tR(t) ) ) );
      }

      // calculate AC selectivity
      S50_tA = mfexp( log_S50_A );
      S95_tA = S50_tA + mfexp( log_S95_step_A );
      // time loop to fill
      for( int t=1; t<=nT; t++ )
      {
        dvariable tmp = log(19.)/( S95_tA(t) - S50_tA(t) );
        selA_ta(t) = 1./( 1. + exp(-tmp*(age - S50_tA(t) ) ) );
      }
  }


FUNCTION popInit
  {
  int g;
  Bta.initialize(); SSBt.initialize(); N4p.initialize();
  uCgatF.initialize(); uCgatC.initialize(); uCgatX.initialize();
  expBtg3P.initialize(); B510.initialize(); N510.initialize();
  avgF.initialize(); expBtgT.initialize(); uCgatR.initialize();
  uCgatA.initialize(); Zta.initialize(); Ftg.initialize();
  SSBa.initialize();
  // Initialize abundance in first year assuming each cohort
  // (except plusGroup) had independent recruitment deviations
  // and age-dependent M.
  // age-1 does not involve M

  for (int p=1; p<=nP; p++)
  {
    Nta(p)(1)(firstAge) = avgR(p)*mfexp( init_recDevs(p,firstAge) );
    for( int a=firstAge+1; a<plusGroupAge; a++ )
    {
      Nta(p)(1)(a) = avgR(p)*mfexp( init_recDevs(p,a) - sum( selF_ta(p)(1)(firstAge,a-1)*Finit(p) + Mta(p)(1)(firstAge,a-1) ) );
    }
    Nta(p)(1)(plusGroupAge) = avgR(p)*mfexp( init_recDevs(p,plusGroupAge) - sum( selF_ta(p)(1)(firstAge,plusGroupAge-1)*Finit(p) + Mta(p)(1)(firstAge,plusGroupAge-1) ) );
    Nta(p)(1)(plusGroupAge) *= 1./( 1.-mfexp(-( selF_ta(p)(1)(plusGroupAge)*Finit(p) + Mta(p)(1)(plusGroupAge ) ) ));

    // Biomass-at-age and SSB:
    Bta(p)(1)   = elem_prod( Nta(p)(1),column(BwtAge(p),1) );
    SSBt(p)(1)  = sum( elem_prod( Bta(p)(1),matAge ) );
    N4p(p)(1) = sum(Nta(p)(1)(4,plusGroupAge));
    B510(p)(1) = sum(Bta(p)(1)(5,10));
    N510(p)(1) = sum(Nta(p)(1)(5,10));

    // Solve catch equations, t=1: returns Zta(nP,1) and Ftg(1)(1,nP)
    int t=1;
    Zta(p)(1) = selF_ta(p)(t)*Finit(p) + Mta(p)(1);
    solveBaranov(p,t,baranovIter);

	// calculate Aug 1 SSB
    for( int k=4; k<=11; k++)
        SSBa(p,1) += Nta(p,1,k)*CwtAge(p,k,1)*exp(-fracAug1*M);


    // omit expBtg calculation for cpue and xnet for t<9 or t<25

    // predicted age-proportions uCgat in this year's catch
    // based on numbers-at-age
    // omit cpue and xnet for t<9 or t<25
    // fishery only in year 1

    g = 1;
    dvar_vector N   = Nta(p)(1)(minAge(g),maxAge(g));
    dvar_vector S   = selF_ta(p)(t)(minAge(g),maxAge(g));
    dvar_vector tmp = elem_prod( N, S );
    for( int k=minAge(g); k<=maxAge(g); k++ )
      uCgatF(p)(k)(1) = tmp(k)/sum(tmp);
   }

  }

FUNCTION popDynamics
  {
  int t, g;
  for (int p=1; p<=nP; p++)
  {
    for( t=2; t<=nT; t++ )
    {
      // Age-2 recruitment
      if( recType==0 )
       Nta(p,t,firstAge) = mfexp( gamma_R*log(Nta(p,t-1,firstAge)) + (1.-gamma_R)*log_AvgR(p) + recDevs(p,t) );

          //if( recType==1 )
          //{
          //  rw_recDevs(t) = gamma_R*rw_recDevs(t-1) + recDevs(t);
          //  Nta(t,1)      = mfexp(recRate+rw_recDevs(t))*SSBt(t-1);
          //}

      // age-3 to age-(A-1)
      for( int a=firstAge+1; a<=plusGroupAge-1; a++ )
      {
        Nta(p,t,a) = Nta(p,t-1,a-1)*exp( -Zta(p,t-1,a-1) );
      }
      Nta(p,t,plusGroupAge) = Nta(p,t-1,plusGroupAge-1)*exp( -Zta(p,t-1,plusGroupAge-1) ) +
                            Nta(p,t-1,plusGroupAge)*exp( -Zta(p,t-1,plusGroupAge) );

      // Biomass-at-age and SSB
      Bta(p)(t)  = elem_prod( Nta(p)(t),column(BwtAge(p),t));
      SSBt(p,t) = sum( elem_prod( Bta(p)(t), matAge ) );
      N4p(p,t) = sum(Nta(p)(t)(4,plusGroupAge));
      B510(p,t) = sum(Bta(p)(t)(5,10));
      N510(p,t) = sum(Nta(p)(t)(5,10));
       //cout << "start baranov" <<endl;
      // Solve catch equations, t: returns Zta(t) and Ftg(t)(p)
      solveBaranov(p,t,baranovIter);
       //cout << "end baranov" <<endl;

	 // predicted exploitable biomass by index
     // cpue
     if(t >= 9)
     {
       g = idxIndex(1);
       for( int k=minAge(g); k<=maxAge(g); k++ )
         expBtg3P(1,p,t) += Nta(p,t,k)*FwtAge(p,k,t)*selC_ta(p,t,k)*exp(-fracYearSurvey(1)*Zta(p,t,k));
     }
      //cout << "done expBtg cpue" <<endl;
     // xnet
     if(t >= 38)
      {
       g = idxIndex(2);
       for( int k=minAge(g); k<=maxAge(g); k++ )
         expBtg3P(2,p,t) += Nta(p,t,k)*XwtAge(p,k,t)*selX_ta(p,t,k)*exp(-fracYearSurvey(2)*Zta(p,t,k));
      }
      //cout << "done expBtg xnet" <<endl;

     if(t >= 17)
     {
       g = 4;
       for( int k=minAge(g); k<=maxAge(g); k++ )
         expBtgT(1,t) += Nta(p,t,k)*RwtAge(k,t)*selR_ta(t,k)*exp(-fracYearSurvey(1)*Zta(p,t,k));
     }
      //cout << "done expBtg rv" <<endl;

     if(t >= 17 & t < 45)
     {
       g = 5;
       for( int k=minAge(g); k<=maxAge(g); k++ )
         expBtgT(2,t) += Nta(p,t,k)*AwtAge(k,t)*selA_ta(t,k)*exp(-fracYearSurvey(1)*Zta(p,t,k));
     }
      //cout << "done expBtg AC" <<endl;

	  // calculate Aug 1 SSB
      for( int k=4; k<=11; k++)
         SSBa(p,t) += Nta(p,t,k)*CwtAge(p,k,t)*exp(-fracAug1*M);

      // predicted age-proportions in this year's catch
      // Fishery
      g=1;
      dvar_vector N   = Nta(p)(t)(minAge(g),maxAge(g));
      dvar_vector S   = selF_ta(p)(t)(minAge(g),maxAge(g));
      dvar_vector tmp = elem_prod( N, S );
      for( int k=minAge(g); k<=maxAge(g); k++ )
        uCgatF(p,k,t) = tmp(k)/sum(tmp);
        //cout << "done uCgatF" <<endl;

      // cpue
      if(t >= 9)
      {
        g=2;
        dvar_vector Nc   = Nta(p)(t)(minAge(g),maxAge(g));
        dvar_vector Sc   = selC_ta(p)(t)(minAge(g),maxAge(g));
        dvar_vector tmpc = elem_prod( Nc, Sc );
        for( int k=minAge(g); k<=maxAge(g); k++ )
          uCgatC(p,k,t) = tmpc(k)/sum(tmpc);
      }
        //cout << "done uCgatC" <<endl;

      // xnet
      if(t >= 38)
      {
        g=3;
        dvar_vector Nx   = Nta(p)(t)(minAge(g),maxAge(g));
        dvar_vector Sx   = selX_ta(p)(t)(minAge(g),maxAge(g));
        dvar_vector tmpx = elem_prod( Nx, Sx );
        for( int k=minAge(g); k<=maxAge(g); k++ )
          uCgatX(p,k,t) = tmpx(k)/sum(tmpx);
          //cout << "done uCgatX" <<endl;
      }
    } // end t-loop
  } // end p-loop

      // rv
     for( t=17; t<=nT; t++ )
     {
        g=4;
        dvar_vector Nr = Nta(1)(t)(minAge(g),maxAge(g));
        Nr += Nta(2)(t)(minAge(g),maxAge(g));
        Nr += Nta(3)(t)(minAge(g),maxAge(g));
        dvar_vector Sr   = selR_ta(t)(minAge(g),maxAge(g));
        dvar_vector tmpr = elem_prod( Nr, Sr );
        for( int k=minAge(g); k<=maxAge(g); k++ )
          uCgatR(k,t) = tmpr(k)/sum(tmpr);
     }
        //cout << "done uCgatR" <<endl;

      // ac
     for( t=17; t<=44; t++ )
     {
        g=5;
        dvar_vector Nac = Nta(1)(t)(minAge(g),maxAge(g));
        Nac += Nta(2)(t)(minAge(g),maxAge(g));
        Nac += Nta(3)(t)(minAge(g),maxAge(g));
        dvar_vector Sac   = selA_ta(t)(minAge(g),maxAge(g));
        dvar_vector tmpa = elem_prod( Nac, Sac );
        for( int k=minAge(g); k<=maxAge(g); k++ )
          uCgatA(k,t) = tmpa(k)/sum(tmpa);
      }
        //cout << "done uCgatA" <<endl;

  if ( last_phase() ) neval+=1;
  else neval=-2.;  // start neval counter at -2 at start of last phase so it equals admb screen output
  }

FUNCTION getTotals
  {
    for (int t=1; t<=nT; t++)
    {
      TSSB(t) = SSBt(1,t) + SSBt(2,t) + SSBt(3,t);
	    TaSSB(t) = SSBa(1,t) + SSBa(2,t) + SSBa(3,t);
      TN4p(t) = N4p(1,t) +  N4p(2,t) + N4p(3,t);
      TN4(t) = N4(1,t) +  N4(2,t) + N4(3,t);
      TN2(t) = N2(1,t) +  N2(2,t) + N2(3,t);
    }
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
     {
       for( int p=1; p<=3; p++)
         qPrior += 0.5*norm2(Qdev(p))/square(QdevSD);
     }
    }

FUNCTION calc_M_prior
  {
  mPrior = 0.;
  if( active(mDevN) )
  {
    // Normal prior for M random walk residuals
    for( int j=1; j<=nSplitAge; j++ )
     {
      mPrior +=  0.5*norm2(mDevN(j))/square(priorSD_rwM(j));
      mPrior +=  0.5*norm2(mDevM(j))/square(priorSD_rwM(j));
      mPrior +=  0.5*norm2(mDevS(j))/square(priorSD_rwM(j));
     }
  }
  if( active(log_M) )
  {
  // Normal prior for initial M
  for( int p=1; p<=3; p++ )
   {
    for( int j=1; j<=nSplitAge; j++ )
     {
       mPrior += 0.5*square(mfexp(M_init(p,j)) - prior_initM(j))/square(priorSD_initM(j));
     }
   }
  }
  }


FUNCTION calc_index_likelihood
  {
  // Initialize likelihood function terms.
  int i,t,p;
  indexLikelihood.initialize();
  subNdxLike.initialize();
  predNdx3P.initialize();
  predNdxT.initialize();
  validObs3P.initialize();
  validObsT.initialize();
  ss3P.initialize();
  ssT.initialize();
  dvar_vector z(1,nT);  // residuals
  dvariable tmp;
  //cout<< " starting calc likelihood "<< i  <<endl;

  for ( i=1; i<=nIndex3P; i++ )
  {
    for ( p=1; p<=nP; p++)
    {
      z.initialize();
      for ( t=1; t<=nT; t++ )
      {
        predNdx3P(i,p,t) = q3P(i,p,t)*expBtg3P(i,p,t);
        if ( idxSeries3P(i,p,t) > 0. )
        {
          z[t] = log(idxSeries3P(i,p,t)) - log( predNdx3P(i,p,t) );
          ss3P[i][p] += pow(z[t], 2);
          validObs3P[i][p] += 1;
        }
      }
      subNdxLike(i,p) = validObs3P(i,p)*log_sig3P(i,p) + (ss3P(i,p)/(2.*square(sigIndex3P(i,p))));
    } // end p loop
  }   // end index series loop
  indexLikelihood(1) = sum(subNdxLike(1)(1,nP));
  indexLikelihood(2) = sum(subNdxLike(2)(1,nP));
    for ( i=1; i<=nIndexT; i++ )
    {
      z.initialize();
      for ( t=1; t<=nT; t++ )
      {
        predNdxT(i,t) = qT(i)*expBtgT(i,t);
        if ( idxSeriesT(i,t) > 0. )
        {
          z[t] = log(idxSeriesT(i,t)) - log( predNdxT(i,t) );
          ssT[i] += pow(z[t], 2);
          validObsT[i] += 1;
        }
      }
      tmp = validObsT(i)*log_sigT(i) + (ssT(i)/(2.*square(sigIndexT(i))));
      if(i == 1)
      {
          indexLikelihood(3) = tmp;
      }
      else
      {
          indexLikelihood(4) = tmp;
      }
    }
  }


FUNCTION calc_age_likelihood
  {
  dvariable meanDiff;
  dvariable nYearsAges;
  dvariable etaSumSq;
  dvariable mnLL;
  int g,p;
  tauSquareAges3P.initialize();
  tauSquareAgesT.initialize();
  ageLikelihood.initialize();
  subAgeLike.initialize();

  // 3P paa
  // Calculate predicted age-proportions for each index gear.
  // fishery
  for ( p=1; p<=nP; p++)
  {
    nYearsAges = 0.;
    etaSumSq   = 0.;
    mnLL       = 0.;
    g          = ageIndex(1);
    for ( int t=ageFirstYear(g); t<=ageLastYear(g); t++ )
    {
      dvar_vector res(1,maxAge(g));
      res.initialize();
      nYearsAges += 1;
      dvar_vector obsPat  = column(agePropFish(p),t)(minAge(g),maxAge(g));
      dvar_vector predPat = ( column(uCgatF(p),t)(minAge(g),maxAge(g)) )/sum( column(uCgatF(p),t)(minAge(g),maxAge(g)) );

      // Column means for p_at residuals.
      int iRes = 0; dvariable sumRes = 0.;
      for( int a=minAge(g); a<=maxAge(g); a++ )
      {
        if( obsPat(a) > 0 )
        {
          iRes   += 1;
          res(iRes) = log( obsPat(a) ) - log( predPat(a) );
          sumRes += res(iRes);
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
    tauSquareAges3P(g,p) = etaSumSq/( (maxAge(g)-minAge(g))*nYearsAges );
    subAgeLike(g,p) = (maxAge(g)-minAge(g))*nYearsAges*log(tauSquareAges3P(g,p)) + mnLL;
  }   // END OF P LOOP  for fishery

  // cpue
  for ( p=1; p<=nP; p++)
  {
    nYearsAges = 0.;
    etaSumSq   = 0.;
    mnLL       = 0.;
    g          = ageIndex(2);
    for ( int t=ageFirstYear(g); t<=ageLastYear(g); t++ )
    {
      dvar_vector res(1,maxAge(g));
      res.initialize();
      nYearsAges += 1;
      dvar_vector obsPat  = column(agePropCpue(p),t)(minAge(g),maxAge(g));
      dvar_vector predPat = ( column(uCgatC(p),t)(minAge(g),maxAge(g)) )/sum( column(uCgatC(p),t)(minAge(g),maxAge(g)) );

      // Column means for p_at residuals.
      int iRes = 0; dvariable sumRes = 0.;
      for( int a=minAge(g); a<=maxAge(g); a++ )
      {
        if( obsPat(a) > 0 )
        {
          iRes   += 1;
          res(iRes) = log( obsPat(a) ) - log( predPat(a) );
          sumRes += res(iRes);
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
    tauSquareAges3P(g,p) = etaSumSq/( (maxAge(g)-minAge(g))*nYearsAges );
    subAgeLike(g,p) = (maxAge(g)-minAge(g))*nYearsAges*log(tauSquareAges3P(g,p)) + mnLL;
  }   // END OF P LOOP  for Cpue

  // xnet
  for ( p=1; p<=nP; p++)
  {
    nYearsAges = 0.;
    etaSumSq   = 0.;
    mnLL       = 0.;
    g          = ageIndex(3);
    for ( int t=ageFirstYear(g); t<=ageLastYear(g); t++ )
    {
      dvar_vector res(1,maxAge(g));
      res.initialize();
      nYearsAges += 1;
      dvar_vector obsPat  = column(agePropXnet(p),t)(minAge(g),maxAge(g));
      dvar_vector predPat = ( column(uCgatX(p),t)(minAge(g),maxAge(g)) )/sum( column(uCgatX(p),t)(minAge(g),maxAge(g)) );

      // Column means for p_at residuals.
      int iRes = 0; dvariable sumRes = 0.;
      for( int a=minAge(g); a<=maxAge(g); a++ )
      {
        if( obsPat(a) > 0 )
        {
          iRes   += 1;
          res(iRes) = log( obsPat(a) ) - log( predPat(a) );
          sumRes += res(iRes);
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
    tauSquareAges3P(g,p) = etaSumSq/( (maxAge(g)-minAge(g))*nYearsAges );
    subAgeLike(g,p) = (maxAge(g)-minAge(g))*nYearsAges*log(tauSquareAges3P(g,p)) + mnLL;
  }   // END OF P LOOP for Xnets

    for( int g=1; g<=3; g++) ageLikelihood(g) = sum(subAgeLike(g)(1,nP));

  // T paa
  // RV
    nYearsAges = 0.;
    etaSumSq   = 0.;
    mnLL       = 0.;
    g          = ageIndex(4);
    for ( int t=ageFirstYear(g); t<=ageLastYear(g); t++ )
    {
      dvar_vector res(1,maxAge(g));
      res.initialize();
      nYearsAges += 1;
      dvar_vector obsPat  = column(agePropRV,t)(minAge(g),maxAge(g));
      dvar_vector predPat = ( column(uCgatR,t)(minAge(g),maxAge(g)) )/sum( column(uCgatR,t)(minAge(g),maxAge(g)) );

      // Column means for p_at residuals.
      int iRes = 0; dvariable sumRes = 0.;
      for( int a=minAge(g); a<=maxAge(g); a++ )
      {
        if( obsPat(a) > 0 )
        {
          iRes   += 1;
          res(iRes) = log( obsPat(a) ) - log( predPat(a) );
          sumRes += res(iRes);
        }
        else
        {
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
    tauSquareAgesT(1) = etaSumSq/( (maxAge(g)-minAge(g))*nYearsAges );
    ageLikelihood(4) = (maxAge(g)-minAge(g))*nYearsAges*log(tauSquareAgesT(1)) + mnLL;

  // AC
    nYearsAges = 0.;
    etaSumSq   = 0.;
    mnLL       = 0.;
    g          = ageIndex(5);
    for ( int t=ageFirstYear(g); t<=ageLastYear(g); t++ )
    {
      dvar_vector res(1,maxAge(g));
      res.initialize();
      nYearsAges += 1;
      dvar_vector obsPat  = column(agePropACR,t)(minAge(g),maxAge(g));
      dvar_vector predPat = ( column(uCgatA,t)(minAge(g),maxAge(g)) )/sum( column(uCgatA,t)(minAge(g),maxAge(g)) );

      // Column means for p_at residuals.
      int iRes = 0; dvariable sumRes = 0.;
      for( int a=minAge(g); a<=maxAge(g); a++ )
      {
        if( obsPat(a) > 0 )
        {
          iRes   += 1;
          res(iRes) = log( obsPat(a) ) - log( predPat(a) );
          sumRes += res(iRes);
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
    tauSquareAgesT(2) = etaSumSq/( (maxAge(g)-minAge(g))*nYearsAges );
    ageLikelihood(5) = (maxAge(g)-minAge(g))*nYearsAges*log(tauSquareAgesT(2)) + mnLL;

  }  // end of function


FUNCTION  void solveBaranov(int p, int t, int nIter )
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
  Bprime = elem_prod(elem_prod( Nta(p)(t),column(CwtAge(p),t)), selF_ta(p)(t) );
  Ftg(t,p) = Ctg(t,p)/sum( Bprime );
  ZaNew = Mta(p)(t) + selF_ta(p)(t)*Ftg(t,p);
  // refine F for fisheries
  for( int i=1; i<=nIter; i++ )
  {
    // Total mortality
    Za=ZaNew; ZaNew=Mta(p)(t);
    // Predicted catch given F
    tmp    = elem_div( elem_prod( Bprime*Ftg(t,p),1.-exp(-Za) ), Za );
    //cout << "predCatch = " << tmp << endl;
    // Function value: difference of pred - obs catch
    f =  sum(tmp) - Ctg(t,p); tmp=0;
    // Jacobian
    dvar_vector tmp1 = elem_div( Bprime, Za );
    dvar_vector tmp2 = elem_prod( selF_ta(p)(t), Za )*Ftg(t,p);
    dvar_vector tmp3 = 1. - mfexp( -Za );
    dvar_vector tmp4 = elem_prod( elem_div( selF_ta(p)(t)*Ftg(t,p), Za ), 1.-exp(-Za)  );

    tmp = elem_prod( tmp1, tmp2 + tmp3 - tmp4 );
    J = sum(tmp); tmp=0;
    Ftg(t,p) -= f/J;
    ZaNew += selF_ta(p)(t)*Ftg(t,p);
    //cout <<"t = " << t<< " iter = "<< i << " f = "<< f <<" J = "<< J << " f/J = " << f/J << endl;
    //cout <<"iter = "<< i << " Ftg = "<< Ftg(t, p) << endl;
  }
  Zta(p,t) = ZaNew;
  avgF(p,t) = sum( elem_prod( Nta(p)(t)(a1F,a2F),selF_ta(p)(t)(a1F,a2F)*Ftg(t,p) ) );
  avgF(p,t) /= sum(Nta(p)(t)(a1F,a2F));
  RETURN_ARRAYS_DECREMENT();
  }


REPORT_SECTION
  {
  // Parameter estimates
  report << "## Parameter estimates" << endl;
  report << "# avgR "<< endl;        report << avgR << endl;
  report << "# init_recDevs"<< endl; report << init_recDevs << endl;
  report << "# recDevs"<< endl; report << recDevs << endl;
  report << "# SSB0" << endl;        report << SSBt(1) << endl;
  report << "# SSBa" << endl;       report << SSBa << endl;
  report << "# Mn" << endl;           report << mfexp(Mt(1)) << endl;
  report << "# Mm" << endl;           report << mfexp(Mt(2)) << endl;
  report << "# Ms" << endl;           report << mfexp(Mt(3)) << endl;
  report << "# mDevN" << endl;           report << mDevN << endl;
  report << "# mDevM" << endl;           report << mDevM << endl;
  report << "# mDevS" << endl;           report << mDevS << endl;

  report << "## Selectivity parameters: S50_Fp" << endl;
  for( int p=1;p<=nP;p++ )
  {
    for( int j=1; j<=nSelBlocks(1); j++ )
    {
      report << "# S50_F"<< p << j << endl;
      report << mfexp(log_S50_F(p,j)) << endl;
    }
  }
  report << "## Selectivity parameters: S95_Fp" << endl;
  for( int p=1;p<=nP;p++ )
  {
    for( int j=1; j<=nSelBlocks(1); j++ )
    {
      report << "# S95_F"<< p << j << endl;
      report << mfexp(log_S50_F(p,j))+mfexp(log_S95_step_F(p,j)) << endl;
    }
  }
  report << "## Selectivity parameters: S50_Cp" << endl;
  for( int p=1;p<=nP;p++ )
  {
      report << "# S50_C"<< p << endl;
      report << mfexp(log_S50_C(p)) << endl;
  }
  report << "## Selectivity parameters: S95_Cp" << endl;
  for( int p=1;p<=nP;p++ )
  {
      report << "# S95_C"<< p << endl;
      report << mfexp(log_S50_C(p))+mfexp(log_S95_step_C(p)) << endl;
  }
  report << "## Selectivity parameters: S50_Xp" << endl;
  for( int p=1;p<=nP;p++ )
  {
      report << "# S50_X"<< p << endl;
      report << mfexp(log_S50_X(p)) << endl;
  }
  report << "## Selectivity parameters: S95_Xp" << endl;
  for( int p=1;p<=nP;p++ )
  {
      report << "# S95_X"<< p << endl;
      report << mfexp(log_S50_X(p))+mfexp(log_S95_step_X(p)) << endl;
  }
  report << "# S50_R" << endl; report <<mfexp(log_S50_R) << endl;
  report << "# S95_R"<< endl; report << mfexp(log_S50_R) + mfexp(log_S95_step_R) << endl;
  report << "# S50_A" << endl; report <<mfexp(log_S50_A) << endl;
  report << "# S95_A"<< endl; report << mfexp(log_S50_A) + mfexp(log_S95_step_A) << endl;
  report << "# Finit" << endl;       report << Finit << endl;
  report << "# gamma_R" << endl;     report << gamma_R << endl;
  report << "# logq3P" << endl;        report << log_q3P << endl;
  report << "# logqT" << endl;        report << log_qT << endl;
  report << "# qC" << endl;        report << q3P(1) << endl;
  report << "# qX" << endl;        report << q3P(2) << endl;
  report << "# qR" << endl;        report << qT(1) << endl;
  report << "# qA" << endl;        report << qT(2) << endl;

  // Derived likelihood components
  report << "## Standard error estimates" << endl;
  report << "# sigIndex3P" << endl;        report << log_sig3P << endl;
  report << "# sigIndexT" << endl;        report << log_sigT << endl;
  report << "# tauAge3P" << endl;     report << sqrt(tauSquareAges3P) << endl;
  report << "# tauAgeT" << endl;     report << sqrt(tauSquareAgesT) << endl;

  // Minimization performance
  report << endl;
  report << "## Minimization performance" << endl;
  report << "# indexLikelihood" << endl;  report << indexLikelihood << endl;
  report << "# subNdxLike" << endl;  report << subNdxLike << endl;
  report << "# ageLikelihood" << endl;    report << ageLikelihood << endl;
  report << "# subAgeLike" << endl;  report << subAgeLike << endl;
  report << "# recPrior" << endl;    report << recPrior << endl;
  report << "# qPrior" << endl;    report << qPrior << endl;
  report << "# mPrior" << endl;    report << mPrior << endl;
  report << "# objFun" << endl;       report << *objective_function_value::pobjfun << endl;
  report << "# maxGrad" << endl;      report << objective_function_value::gmax << endl;
  report << "# exitCode" << endl;     report << iexit << endl;
  report << "# funEvals" << endl;     report << neval << endl;

  // Derived variables
  report << endl;
  report << "## Derived variables" << endl;
  report << "# TSSB" << endl;       report << TSSB << endl;
  report << "# SSBt" << endl;       report << SSBt << endl;

  report << "# TN4p" << endl;       report << TN4p << endl;
  report << "# N4p" << endl;       report << N4p << endl;

  report <<"## Fishing mortality rates"<< endl;
  report << "# avgF"<< endl;        report << avgF << endl;
  report << "# FtgN" << endl;        report << column(Ftg,1)<< endl;
  report << "# FtgM" << endl;        report << column(Ftg,2)<< endl;
  report << "# FtgS" << endl;        report << column(Ftg,3)<< endl;
  report << "# Zta" << endl;         report << Zta << endl;

  report << "# NtaN" << endl;        report << Nta(1) << endl;
  report << "# NtaM" << endl;        report << Nta(2) << endl;
  report << "# NtaS" << endl;        report << Nta(3) << endl;

  report << "# BtaN" << endl;        report << Bta(1) << endl;
  report << "# BtaM" << endl;        report << Bta(2) << endl;
  report << "# BtaS" << endl;        report << Bta(3) << endl;
  report << "# B510" << endl;        report << B510 << endl;
  report <<"## Exploitable biomass by cpue (vals < 0 are missing...)"<< endl;
  report << "# expBtgC"<< endl;
  for( int p=1;p<=nP;p++ )
  {
    report << "# expBtgC" << p << endl;
    report << expBtg3P(1)(p) << endl;
  }
  report << "# expBtgX"<< endl;
  for( int p=1;p<=nP;p++ )
  {
    report << "# expBtgX" << p << endl;
    report << expBtg3P(2)(p) << endl;
  }
  report << "# expBtgR" << endl; report << expBtgT(1) << endl;
  report << "# expBtgA" << endl; report << expBtgT(2) << endl;
  report <<"## Predicted Indices (vals < 0 are missing...)"<< endl;
  for( int p=1;p<=nP;p++ )
  {
    report << "# predNdxC" << p << endl;
    report << predNdx3P(1)(p) << endl;
  }
  for( int p=1;p<=nP;p++ )
  {
    report << "# predNdxX" << p << endl;
    report << predNdx3P(2)(p) << endl;
  }
  report << "# predNdxR" << endl; report << predNdxT(1)<< endl;
  report << "# predNdxA" << endl; report << predNdxT(2)<< endl;

  report <<"## Predicted age proportions by fishery)"<< endl;
  for( int p=1;p<=nP;p++ )
  {
    report << "# uCgatF"<< p << endl;
    report << uCgatF(p) << endl;
  }
  for( int p=1;p<=nP;p++ )
  {
    report << "# uCgatC"<< p << endl;
    report << uCgatC(p) << endl;
  }
  for( int p=1;p<=nP;p++ )
  {
    report << "# uCgatX"<< p << endl;
    report << uCgatX(p) << endl;
  }
  report << "# uCgatR" << endl; report << uCgatR << endl;
  report << "# uCgatA" << endl; report << uCgatA << endl;

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

  report << "# a1F" << endl;  report << a1F << endl;
  report << "# a2F" << endl;  report << a2F << endl;

  report << "# ph_S50" << endl;  report << ph_S50 << endl;
  report << "# ph_S95" << endl;  report << ph_S95 << endl;
  report << "# ph_gammaR" << endl;  report << ph_gammaR << endl;

  // Model parameter priors.
  report << "# priorSD_R" << endl;       report << priorSD_R << endl;

  // Inputs from scaCatch.dat.
  report << "# nT" << endl;              report << nT << endl;
  report << "# nFisheries" << endl;      report << nFisheries << endl;
  report << "# landCatchMatrix" << endl; report << landCatchMatrix << endl;

  // Inputs from scaIndex.dat.
  report << "# nIndex3P" << endl;    report << nIndex3P << endl;
  report << "# nIndexT" << endl;    report << nIndexT << endl;
  report << "# nIndex" << endl;    report << nIndex << endl;
  report << "# idxIndex" << endl;        report << idxIndex << endl;
  report << "# idxRelative" << endl;     report << idxRelative << endl;
  report << "# idxLikeWeight" << endl;   report << idxLikeWeight << endl;
  report << "# idxFirstYear" << endl;    report << idxFirstYear << endl;
  report << "# idxLastYear" << endl;     report << idxLastYear << endl;
  report << "# fracYearSurvey" << endl;  report << fracYearSurvey << endl;
  report << "## idxSeriesC" << endl;
  for ( int p=1; p<=nP;p++ )
  {
    report << "# idxSeriesC" << p << endl;
    report << idxSeries3P(1)(p) << endl;
  }
  report << "## idxSeriesX" << endl;
  for ( int p=1; p<=nP;p++ )
  {
    report << "# idxSeriesX" << p << endl;
    report << idxSeries3P(2)(p) << endl;
  }
  report << "# idxSeriesR" << endl;  report << idxSeriesT(1) << endl;
  report << "# idxSeriesA" << endl;  report << idxSeriesT(2) << endl;

  // Inputs from scaAges.dat.
  report << "# firstAge" << endl;        report << firstAge << endl;
  report << "# plusGroupAge" << endl;    report << plusGroupAge << endl;
  report << "# nAgeSeries" << endl;      report << nAgeSeries << endl;
  report << "# minAge" << endl;          report << minAge << endl;
  report << "# maxAge" << endl;          report << maxAge << endl;
  report << "# ageIndex" << endl;        report << ageIndex << endl;
  report << "# ageFirstYear" << endl;    report << ageFirstYear << endl;
  report << "# ageLastYear" << endl;     report << ageLastYear << endl;
  report << "## ageObsProp Fishery" << endl;
  for ( int p=1; p<=nP;p++ )
  {
    report << "# agePropFish" << p << endl;
    report << agePropFish(p) << endl;
  }
  report << "## ageObsProp Cpue" << endl;
  for ( int p=1; p<=nP;p++ )
  {
    report << "# agePropCpue" << p << endl;
    report << agePropCpue(p) << endl;
  }
  report << "## ageObsPropXnet" << endl;
  for ( int p=1; p<=nP;p++ )
  {
    report << "# agePropXnet" << p << endl;
    report << agePropXnet(p) << endl;
  }
  report << "# agePropRV" << endl;  report << agePropRV << endl;
  report << "# agePropACR" << endl;  report << agePropACR << endl;

  report << "# init_tauAge3P" << endl;    report << init_tauAge3P << endl;
  report << "# init_tauAgeT" << endl;    report << init_tauAgeT << endl;

  report << "# ages" << endl;          report << age << endl;
  report << "# M" << endl;          report << M << endl;
  for ( int p=1; p<=nP;p++ )
  {
    report << "# BwtAge" << p << endl;
    report << BwtAge(p) << endl;
  }
  for ( int p=1; p<=nP;p++ )
  {
    report << "# CwtAge" << p << endl;
    report << CwtAge(p) << endl;
  }
  for ( int p=1; p<=nP;p++ )
  {
    report << "# FwtAge" << p << endl;
    report << FwtAge(p) << endl;
  }
  for ( int p=1; p<=nP;p++ )
  {
    report << "# XwtAge" << p << endl;
    report << XwtAge(p) << endl;
  }
  report << "# RwtAge" << endl;  report << RwtAge << endl;
  report << "# AwtAge" << endl;  report << AwtAge << endl;
  report << "# matAge" << endl;             report << matAge << endl;
  report <<"## Fishery Selectivity function"<< endl;
  for( int p=1;p<=nP;p++ )
  {
    report << "# selF_ta"<<p<< endl; report << selF_ta(p) << endl;
  }
  report <<"## CPUE Selectivity function"<< endl;
  for( int p=1;p<=nP;p++ )
  {
    report << "# selC_ta"<<p<< endl; report << selC_ta(p) << endl;
  }
  report <<"## Xnet Selectivity function"<< endl;
  for( int p=1;p<=nP;p++ )
  {
    report << "# selX_ta"<<p<< endl; report << selX_ta(p) << endl;
  }
  report << "# selR_ta" << endl;  report << selR_ta << endl;
  report << "# selA_ta" << endl;  report << selA_ta << endl;

  // Inputs from scaSelN.dat.
  report << "# p258" << endl;        report << p258 << endl;
  report << "# y1cal" << endl;    report << y1cal << endl;
  report << "# yTcal" << endl;      report << yTcal << endl;
  report << "# a1cal" << endl;          report << a1cal << endl;
  report << "# aTcal" << endl;          report << aTcal << endl;
  report << "# rsel1" << endl;        report << rsel1 << endl;
  report << "# rsel2" << endl;    report << rsel2 << endl;
  report << "# rsel3" << endl;    report << rsel3 << endl;
  //--------------------------------------------------------------------------//
  // End echo of inputs                                                       //
  //--------------------------------------------------------------------------//
  }
