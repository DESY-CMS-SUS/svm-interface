/*libsvm interface for the one and two class HEP problems */
#ifndef svm_interface_h
#define svm_interface_h
#include "libsvm-weights-3.18/svm.h"
#include "libsvm_container.h"
#include <cstdlib>
#include <string>
#include <iostream>
#include <TH1D.h>
#include <TFile.h>
#include "fom.h"
//svm problem container class
class svm_problem_container {
public:
  void set_training_samp (const svm_problem & set){sample_tra=set;};
  void set_test_samp (const svm_problem & set){sample_test=set;};
  void set_parameter_samp (const svm_parameter & set){para=set;};
private:
  svm_model * model;
  svm_problem sample_tra;
  svm_problem sample_eval;
  svm_problem sample_test;
  svm_parameter para;
};
//abstract class svm interface
class svm_interface {
protected:
  std::string svm_type;
  svm_problem_container *svm_pro;
  static const int hist_nbin = 400;
  static const int max_iter = 20;
  unsigned nsamp_tot, nsamp_tra, nsamp_test, nsamp_eval, nbkg_tot, nsig_tot , nbkg_tra, nbkg_test, nsig_tra, nsig_test;
  double nsig_test_w,nbkg_test_w;
  svm_interface(std::string type,unsigned n_tot, unsigned bkg_tot, unsigned sig_tot):svm_type(type),nsamp_tot(n_tot),nbkg_tot(bkg_tot),nsig_tot(sig_tot){
    svm_pro = new svm_problem_container;
    highest_accur_gamma = 0.0001;
    highest_accur_C = 0.1;
 };
  int svm_node_max;
  svm_problem sample_tra;
  svm_problem sample_test;
  svm_problem sample_eval;
  svm_parameter para; 
  bool parameters_set;
  
  //index vectors for the two class problem
  std::vector<double> * y1;
  std::vector<double> * y2;
  std::vector<double> * y;
  double* y_array;
public:

  enum samp_type {TRAIN,TEST,EVALUATE};

  double highest_accur_gamma;
  double highest_accur_C;

  TH1D* disc_S;
  TH1D* disc_B;
  TH1D* cuteff;
  //  virtual double maxSignificance(TH1D*, TH1D*, bool) = 0;
  virtual void set_sample(const svm_container&, samp_type)=0;
  virtual bool set_sample_random_split(svm_container &svm)=0;  
  virtual double para_scan(double,double)=0;
  virtual double para_scan_omp(double,double)=0;
  virtual double train_test(double,double)=0;
  virtual void set_parameters()=0;  
  virtual void gen_class_index()=0;
  virtual void set_indexes()=0;
  virtual ~svm_interface(){ delete svm_pro;};
  virtual void obtain_probabilities(double,double)=0;
};
//c-svc
class csvc_interface:public svm_interface {
private:
  void set_gamma_array(std::vector<double, std::allocator<double> >&, unsigned int);
  
public:
  virtual void set_sample(const svm_container&, samp_type);
  //  virtual double maxSignificance(TH1D*, TH1D*, bool info = false);
  //  virtual void set_sample(const svm_container&, samp_type); 
  virtual void obtain_probabilities(double,double);
  csvc_interface(unsigned n_tot, unsigned bkg_tot, unsigned sig_tot):svm_interface("c-svc",n_tot,bkg_tot,sig_tot){
    nbkg_tra = 0;
    nbkg_test = 0;
    nsig_tra = 0;
    nsig_test = 0;
    parameters_set = false;
    sample_tra.x = (svm_node **) calloc(nsamp_tot,sizeof(svm_node *));
    sample_test.x = (svm_node **) calloc(nsamp_tot,sizeof(svm_node *));
    sample_tra.W = (double*) calloc (nsamp_tot,sizeof(double));
    sample_test.W = (double*) calloc (nsamp_tot,sizeof(double));
  };
  virtual void gen_class_index () {
    if(nbkg_tra < 1) {std::cout<< "ERROR! while generating the class indexes: # of background events is 0 " << std:: endl; exit(2); }
    if(nsig_tra < 1) {std::cout<< "WARNING! while generating the class indexes: # of signal events is 0 " << std:: endl;}
    y1 = new std::vector<double>(nbkg_tra,1.); 
    y2 = new std::vector<double>(nsig_tra,2.);
    y = new std::vector<double>;
    //merging the class labels for background (-1) and signal (1) 
    y->insert(y->begin(),y2->begin(),y2->end());
    y->insert(y->begin(),y1->begin(),y1->end());
    //checking the size of the vectors
    if(y->size() != nsamp_tra)  
      {std::cout << "ERROR! while generating the class indexes: Check size of y (vector of labels) "  << y->size() << " total events reserved for the training: "  <<  nsamp_tra << std::endl; exit(2);}
  };
  virtual void set_indexes(){
    gen_class_index();
    sample_tra.l=(unsigned)nsamp_tra;
    y_array =(double*)calloc(nsamp_tra,sizeof(double));
    for(int indy = 0; indy < nsamp_tra; indy ++)
      y_array[indy] = y->at(indy);
    sample_tra.y= y_array;
  }
  virtual double para_scan(double,double);
  virtual bool set_sample_random_split(svm_container &svm);  
  virtual double para_scan_omp(double,double);
  void set_para_gamma(double g = 0.1){ if(!parameters_set) set_parameters() ; para.gamma=g;};
  void set_para_c(double C = 2.){ if(!parameters_set) set_parameters(); para.C=C;};
  void set_para_kernel(int type = RBF){ if(!parameters_set) set_parameters(); para.kernel_type=type;};     
  svm_parameter* get_parameters(double C, double g) { if(!parameters_set) set_parameters(); para.C=C;para.gamma=g; return &para;};
  void deep_copy_svm_pro(const svm_problem&,int,int,int,svm_problem&);
  virtual double train_test(double,double);

  virtual void set_parameters() {
    //setting svm model parameters
    parameters_set = true;
    set_para_gamma();
    set_para_c();
    set_para_kernel();
    para.svm_type = C_SVC;
    para.kernel_type = RBF;
    para.coef0 = 1.;
    para.eps = 0.001;
    para.nr_weight = 0;
    para.degree = 0;
    para.cache_size = 3000.;
    para.shrinking = 1;
    para.probability = 1;
    svm_pro->set_parameter_samp(para);
  };
  ~csvc_interface(){
    for(unsigned indx= 0 ; indx < nsamp_tra; indx++)
      free(sample_tra.x[indx]);
    free(sample_tra.x);
    for(unsigned indx= 0 ; indx < nsamp_test; indx++)
      free(sample_test.x[indx]);
    free(sample_test.x);
    free(y_array);
    delete y, y1, y2;
  }

  void clean(svm_problem& tbc);
};
  
class svm_analyze{
private:
  svm_interface * given_svm;  
  bool svm_inter_set;
public:
  svm_analyze():svm_inter_set(false){};
  void set_svm_interface(svm_interface *svm_int) {given_svm = svm_int;svm_inter_set=true;};
  bool setup_svm(svm_container& svm){
    if(!svm_inter_set) {std::cout<<"ERROR! svm_interface is not set! "<< std::endl;exit(2);}
    given_svm->set_parameters();
    given_svm->set_sample_random_split(svm);
    given_svm->set_indexes();
    return true;
  };
  void set_eval(const svm_container &eval){
    given_svm->set_sample(eval,svm_interface::EVALUATE);
  };
  double scan_parameters(){

    std::cout << " max accur " <<  given_svm->para_scan_omp(0.8,0) << std::endl; 
    obtain_probabilities(given_svm->highest_accur_C,given_svm->highest_accur_gamma);
  }
  bool obtain_probabilities(double C, double g){
    given_svm->obtain_probabilities(C,g);
    TFile *f = new TFile("svm.root","RECREATE");
    std::cout << "some text " << f << std::endl;
    given_svm->disc_S->Write();
    given_svm->disc_B->Write();
    given_svm->cuteff->Write();
    //    f->Close();
  }
};
#endif
