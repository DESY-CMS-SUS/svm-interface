#include "libsvm-weights-3.18/svm.h"
#include "svm_interface.h"
#include <iostream>
#include <TRandom3.h>
#include <cmath>
#include <TH1D.h>
#include "fom.h"
#include "TMath.h"
#define fEpsilon 1.e-2
double TransformLikelihoodOutput( Double_t ps, Double_t pb )                                                                                                                                                   
{   
  if (ps < fEpsilon) ps = fEpsilon;
  if (pb < fEpsilon) pb = fEpsilon;
  Double_t r = ps/(ps + pb);                                                                                                                                                                                                                 
  if (r >= 1.0) r = 1. - 1.e-2;
  // inverse Fermi function
  // sanity check                                                                                                                                                                                                                         
    if      (r <= 0.0) r = fEpsilon;                                                                                                                                                                                                        
    else if (r >= 1.0) r = 1. - 1.e-2;
    Double_t tau = 2.;// 15.0;
    r = - TMath::Log10(1/r - 1.0)/tau; 
    //  }
  return r;                                                                                                                                                                                                                                  
}                    

void csvc_interface::obtain_probabilities(double c_p , double g_p)
{

  disc_S = new TH1D("svm_disc_signal", "SVM probability signal",25,0.,1.);
  disc_B = new TH1D("svm_disc_background", "SVM probability background",25,0.,1.);
  cuteff = new TH1D("svm_disc_cuteff", "SVM Cut Efficiency",25,0.,1.);

  disc_S->Sumw2();  
  disc_B->Sumw2();  

  double prob[2];
  set_para_gamma(g_p);
  set_para_c(c_p);
  svm_model * csvc_svm_model= (svm_model*)malloc(sizeof(svm_model));
  csvc_svm_model =  svm_train(&(sample_tra),&para);

  for(int comp = 0; comp < nsamp_eval; comp++)
    {
      svm_predict_probability(csvc_svm_model,sample_eval.x[comp],prob);
      disc_S->Fill(prob[1],sample_eval.W[comp]);
    }
  for(int comp = 0; comp < nbkg_test; comp++)
    {
      svm_predict_probability(csvc_svm_model,sample_test.x[comp],prob);
      disc_B->Fill(prob[1],sample_test.W[comp]);
    }

  for(int comp = 0; comp < nbkg_tra; comp++)
    {
      svm_predict_probability(csvc_svm_model,sample_tra.x[comp],prob);
      disc_B->Fill(prob[1],sample_tra.W[comp]);
    }
  fom FOM (0.25);
  std::cout<<" maximum significance with AsimovZ " <<  FOM.maxSignificance(disc_S, disc_B, true, cuteff) << std::endl;  
  /* for(int comp = 0; comp < nsamp_test; comp++)
    {
      svm_predict_probability(csvc_svm_model,sample_test.x[comp],prob);
      if(nbkg_test > comp) disc_B->Fill(prob[1],sample_test.W[comp]);
      else disc_S->Fill(prob[1],sample_test.W[comp]);
      }*/
}

double csvc_interface::train_test(double C, double gamma)
{
 
  svm_model * csvc_svm_model;
  svm_check_parameter(&(sample_tra),get_parameters(C,gamma));
  csvc_svm_model=  svm_train(&(sample_tra),get_parameters(C,gamma));
  double accur_S = 0;
  double accur_B = 0;
  
  for(int comp = 0; comp < nsamp_test; comp++)
    {

      if(fabs(svm_predict(csvc_svm_model, sample_test.x[comp])-2)<0.1 && comp> nbkg_test) accur_S;
      else if(fabs(svm_predict(csvc_svm_model,  sample_test.x[comp])-1)<0.1 && comp < nbkg_test) accur_B++;
    }
  svm_free_and_destroy_model(&csvc_svm_model);  

  return ((double)accur_S)/nsig_test + ((double)accur_B)/nbkg_test;
}

double csvc_interface::para_scan(double C, double pre_accur){ 

  double gamma[] = {0.01, 0.005, 0.001, 0.0005, 0.0001, 0.00005,0.00001,0.000005};
  double accur[8];
  int highest = 0;
  for(int indg = 0; indg < 8; indg++){
    accur[indg] = train_test(C,gamma[indg]);
  }
  for(int ind = 0; ind<7; ind++) {if(accur[ind] > accur[ind+1]) accur[ind+1] = accur[ind]; else highest=ind+1;}

  if(accur[highest] > pre_accur){ para_scan(2.*C,accur[highest]);}
  else {set_para_c(C/2);std::cout << " optimized C " << C/2 << " optimized gamma " << gamma[highest] <<  std::endl;  return pre_accur;}
}

//
void csvc_interface::set_gamma_array(std::vector<double> & gamma_arr, unsigned iter)
{
  static double n = 10; 
  if((iter%4)==0)
    {
      for(int ind = 0; ind < (gamma_arr.size()-1)/2; ind++)
	{
	  gamma_arr.at(ind) = highest_accur_gamma*pow(n,ind+1)/5;
	  gamma_arr.at(gamma_arr.size()-1-ind) = highest_accur_gamma*5/pow(n,ind+1);
	}
      gamma_arr.at((gamma_arr.size()-1)/2) = highest_accur_gamma;
      n = n-2;
    }
}

//omp compatible version of the one above :-O
double csvc_interface::para_scan_omp(double C, double pre_accur){ 
  //  double tmp_gamma[] = { 0.01, 0.005,0.003, 0.001,0.0008, 0.0005 , 0.0003, 0.0001, 0.00005};
  //optimize the gamma
  static int iteration = 0;
  static std::vector<double> gamma(9); 
  //dynamically setting gamma array
  set_gamma_array(gamma,iteration);
  double accur[gamma.size()];
  int highest = 0;
  std::vector<svm_parameter> scan_parameters;
  for(std::vector<double>::iterator g = gamma.begin(); g!= gamma.end();g++){
    scan_parameters.push_back(*get_parameters(C,*g));
  }
  
  svm_problem para_scan_problem_tra;
  deep_copy_svm_pro(sample_tra, nsamp_tra, svm_node_max, nsamp_tra, para_scan_problem_tra);
  svm_problem para_scan_problem_test;
  deep_copy_svm_pro(sample_test, nsamp_test,svm_node_max, 0 ,para_scan_problem_test);
  TH1D* hSig[gamma.size()];
  TH1D* hBkg[gamma.size()];
  
  for(int indg = 0; indg < gamma.size(); indg++){
    hSig[indg] = new TH1D("signal"+indg, "SVM probability signal"+indg,400,0.,1.);
    hBkg[indg] = new TH1D("background"+indg, "SVM probability background"+indg,400,0.,1.);
    hSig[indg]->SetDirectory(0);
    hBkg[indg]->SetDirectory(0);
 }

#pragma omp parallel for shared (para_scan_problem_test,para_scan_problem_tra, hSig, hBkg)
  for(int indg = 0; indg < gamma.size(); indg++){
    svm_model * csvc_svm_model;

    csvc_svm_model=  svm_train(&(para_scan_problem_tra),&(scan_parameters.at(indg)));
    double accur_S = 0.;
    double accur_B = 0.;
    double prob[] = {0,0} ;  
   
    for(int comp = 0; comp < nsamp_test; comp++)
      {
	
	svm_predict_probability(csvc_svm_model,para_scan_problem_test.x[comp],prob);
	if(nbkg_test > comp) hBkg[indg]->Fill(prob[1],para_scan_problem_test.W[comp]);
	else hSig[indg]->Fill(prob[1],para_scan_problem_test.W[comp]);
	/*    
	      if(fabs(svm_predict(csvc_svm_model, para_scan_problem_test.x[comp])-2)<0.1 && comp> nbkg_test) 
	      accur_S+=para_scan_problem_test.W[comp];
	      
	      else if(fabs(svm_predict(csvc_svm_model,  para_scan_problem_test.x[comp])-1)<0.1 && comp < nbkg_test) 
	      accur_B+=para_scan_problem_test.W[comp];
	*/
      }
    bool info = true;

    //    accur[indg] = ((double)accur_S)/nsig_test + ((double)accur_B)/nbkg_test;
    
    svm_free_and_destroy_model(&csvc_svm_model);    
  }
  //assume 25% unc
  fom FOM(0.25);
  //this can be included inside of the above loop but due to OMP thread exclusivity, it is better to run it outside 
  for(int ind = 0; ind < gamma.size(); ind++){
    accur[ind] = FOM.maxSignificance(hSig[ind], hBkg[ind], true);  //accur_S/nsig_test_w + accur_B/nbkg_test_w; /// sqrt(accur_B+pow(accur_B*0.2,2));
    hSig[ind]->Delete();
    hBkg[ind]->Delete();
  }
  for(int ind = 0; ind<gamma.size()-1; ind++) {
    if(accur[ind] > accur[ind+1]) accur[ind+1] = accur[ind]; 
    else highest=ind+1;
  }
  iteration++;
  if(((accur[highest] > pre_accur*0.99) || ( iteration < 5 && accur[highest] > pre_accur*0.9 && accur[highest] < 2.) )&&  iteration < max_iter){
    highest_accur_gamma = gamma[highest];
    std::cout <<" for given C value: " << C << "  highest accuracy: " << accur[highest] << "  obtained with gamma: " << gamma[highest] << "  in the iteration " << iteration << std::endl; 
    para_scan_omp(1.2*C,accur[highest]);
  }
  else {
    std::cout <<" for given C value: " << C << "  highest accuracy: " << accur[highest] << "  obtained with gamma: " << gamma[highest] << "  in the iteration " << iteration << std::endl; 
    highest_accur_C = C/1.2 ;
    std::cout << " optimized C " << highest_accur_C << " optimized gamma " << highest_accur_gamma <<  std::endl;  
    std::cout << "iteration " << iteration << std::endl;
    //  clean(para_scan_problem_tra2);
    return pre_accur;
  }
}

void csvc_interface::clean(svm_problem& tbc)
{
  free(tbc.W);
  free(tbc.x);
}

void csvc_interface::deep_copy_svm_pro(const svm_problem& tbc, int x1_max, int x2_max, int y_max, svm_problem& copy)
{
  copy.x = (svm_node **) calloc(x1_max,sizeof(svm_node*)); 
  for (unsigned indx = 0;indx< x1_max;indx++)
    { 
      copy.x[indx] = (svm_node*) calloc(x2_max,sizeof(svm_node)); 
      for(unsigned indy = 0; indy < x2_max; indy ++) {
	copy.x[indx][indy] = tbc.x[indx][indy]; 
      } 
    }

  copy.y = (double*)calloc(y_max,sizeof(double));
  for(unsigned ind =0; ind < y_max; ind++)
    {
      copy.y[ind] = tbc.y[ind];
    }
  copy.W = (double*)calloc(x1_max,sizeof(double));
  for(unsigned ind =0; ind < x1_max; ind++)
    {
      copy.W[ind] = tbc.W[ind];
    }

  copy.l = tbc.l;
};


bool csvc_interface::set_sample_random_split (svm_container & svm){
  TRandom3 sep;
  //c++ std vectors should be converted to the plain C arrays - allocating memory
  unsigned indtra = 0;
  unsigned indtest = 0;
  //looping over the container
  std::cout << " nsamp_tot " << nsamp_tot << " nbkg tot " << nbkg_tot << "  nsig_tot "<< nsig_tot << std::endl;

  nbkg_test_w = 0.;
  nsig_test_w = 0.;
  
  for(unsigned indx= 0 ; indx < nsamp_tot; indx++){
    //allocating memory for the node array
    if(indx < nbkg_tot)
      {
	if(sep.Binomial(1,0.5) == 0){
	  sample_tra.x[indtra] = (svm_node*) calloc(svm.svm_cont->at(0).size(),sizeof(svm_node));
	  sample_tra.W[indtra] = svm.weights->at(indx);
	  for(unsigned indy = 0; indy <svm.svm_cont->at(indx).size() ; indy ++) {
	    sample_tra.x[indtra][indy] = (svm.svm_cont->at(indx).at(indy)); 

	  } 
	  nbkg_tra++; indtra++;
	}
	else 
	  {
	    nbkg_test_w += svm.weights->at(indx);
	    sample_test.x[indtest] = (svm_node*) calloc(svm.svm_cont->at(0).size(),sizeof(svm_node));
	    sample_test.W[indtest] = svm.weights->at(indx);
	    for(unsigned indy = 0; indy <svm.svm_cont->at(indx).size() ; indy ++) {
	      sample_test.x[indtest][indy] = (svm.svm_cont->at(indx).at(indy));
	    }  
	    indtest++;
	  }
      }
    else
      {
	if(sep.Binomial(1,0.5) == 0){
	  sample_tra.x[indtra] = (svm_node*) calloc(svm.svm_cont->at(0).size(),sizeof(svm_node));
	  sample_tra.W[indtra] = svm.weights->at(indx);	  
	  for(unsigned indy = 0; indy <svm.svm_cont->at(indx).size() ; indy ++) {
	    sample_tra.x[indtra][indy] = (svm.svm_cont->at(indx).at(indy)); 
	  }
	  nsig_tra++; indtra++;
	}
	else 
	  {
	    nsig_test_w += svm.weights->at(indx);
	    sample_test.x[indtest] = (svm_node*) calloc(svm.svm_cont->at(0).size(),sizeof(svm_node));
	    sample_test.W[indtest] = svm.weights->at(indx);
	    for(unsigned indy = 0; indy <svm.svm_cont->at(indx).size() ; indy ++) {
	      sample_test.x[indtest][indy] = (svm.svm_cont->at(indx).at(indy));
	    }  
	    indtest++;
	  }
      }
  } 
  std::cout << nsig_test_w << "    " << nbkg_test_w << std::endl;
  nsamp_tra = nsig_tra+nbkg_tra; 
  svm_node_max = svm.svm_cont->at(0).size();

  nsamp_test = nsamp_tot-nsamp_tra;
  nsig_test = nsig_tot - nsig_tra;
  nbkg_test = nbkg_tot - nbkg_tra;
  svm_pro->set_test_samp(sample_test);
  svm_pro->set_training_samp(sample_tra);
  
  
  return true;
}

void csvc_interface::set_sample (const svm_container & svm, samp_type type)
{

  svm_problem sample;
  sample.x = (svm_node **) calloc(svm.svm_cont->size(),sizeof(svm_node *));
  sample.W = (double*) calloc (svm.svm_cont->size(),sizeof(double));

  unsigned indx= 0;
  for(indx = 0; indx < svm.svm_cont->size(); indx++){ 
    sample.x[indx] = (svm_node*) calloc(svm.svm_cont->at(0).size(),sizeof(svm_node));
    sample.W[indx] = svm.weights->at(indx);	    
    for(unsigned indy = 0; indy <svm.svm_cont->at(indx).size() ; indy ++) {
      sample.x[indx][indy] = (svm.svm_cont->at(indx).at(indy));
    }
  }
  if(EVALUATE == type){
    sample_eval = sample;
    nsamp_eval = indx;
  }
  else if(TRAIN == type){
    sample_tra = sample;
    nsamp_tra = indx;
  }
  else if(TEST == type){
    sample_test = sample;
    nsamp_test = indx;
  }
}
