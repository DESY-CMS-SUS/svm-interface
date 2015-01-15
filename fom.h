//simple class to obrain Asimov significance
#ifndef fom_h
#define fom_h
#include <TH1D.h>
#include <cmath>
#include <iostream>
#include <cstdlib>
static const double epsilon = 1e-35;
class fom{
public:
  enum fom_type{
    asimov,
    stop
  };
  //constructor only with unc initially set
  fom(double uncertain){ unc = uncertain;}
  //destructor
  ~fom(){};
  //set bkg
  void setBackground(double bkg){background = bkg;}
  //set signal
  void setSignal(double sig){signal = sig;}

  double maxSignificance(TH1D* sig, TH1D* bkg, bool info, TH1D* cuteff = 0);

  double getSignificance(fom_type f_type) {if (f_type == asimov) asimovZ();
    else if(f_type == stop) Stop();
    else {std::cout << "wrong fom type, please enter one of the following significance types: asimov, stop" << std::endl; std::exit(0);}
    return significance;}
private:
  //asimovZ with uncertainty: http://www.pp.rhul.ac.uk/~cowan/atlas/cowan_statforum_8may12.pdf
  void asimovZ ();
  void Stop (); //Standard FOM used in sTop analysis
  double significance;
  double unc;
  double signal;
  double background;
  void check_vars();
  fom(){};
};

#endif
