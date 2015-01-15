#include "fom.h"
#include <cmath>
#include <iostream>
#include <cstdlib>
void fom::asimovZ(){
  check_vars();
  double varb = background*unc*background*unc; 
  double tot = signal + background;
  
  double asimovsig = sqrt(2*(tot*log((tot*(varb+background))/((background*background)+tot*varb))-(background*background/varb)*log(1+(varb*signal)/(background*(background+varb)))));
  significance = asimovsig;
}
void fom::Stop(){
  check_vars();
  double varb = background*unc*background*unc; 
  double stopsig = signal / sqrt(background+varb);
  significance = stopsig;
}
void fom::check_vars()
{
  if(signal < epsilon) std::cout << "WARNING! Number of signal events: " << signal << std::endl; 
  if(background < epsilon) std::cout << "WARNING! Number of background events: " << background << std::endl; 
  if(unc < epsilon) std::cout << "WARNING! Uncertainty is set to: " << unc << std::endl;
}
double fom::maxSignificance(TH1D* sig, TH1D* bkg, bool info, TH1D* cuteff )
{
  bool docuteff = false;
  if(cuteff) docuteff = true; 
  float sign = 0.;
  float max_sign = -1.;
  int max_bin = sig->GetSize() - 2; //substracting over-under flow bins
  if(max_bin != bkg->GetSize() - 2) {std::cout << " ERROR! Bin numbers are different. " << std::endl; exit(2);}
  int min_bin = 0;
  if(docuteff) min_bin = 0;
  else min_bin  =   (int)(bkg->GetMean()*(double)max_bin); // background -> 0 , signal -> 1 
  if(info) std::cout<<  " mean bin " << min_bin << "   ";
  float intSig = 0., intBkg = 0.;
  int maxBin = 0;
  for (int ind = min_bin; ind < max_bin ; ind++)
    {
      intSig = sig->Integral(ind,  max_bin);
      intBkg = bkg->Integral(ind,  max_bin);
      //well you have less than 3 signal events left...  
      if(intSig < 3.) 
	{
	  if(info) std::cout << " integrals sig " << intSig << " bkg " << intBkg << "    " ;
	  if(info) std::cout << " cut for the max significance is on bin " << maxBin << " with significance " << max_sign <<  std::endl;
	  return max_sign;
	}
      this->setSignal(intSig); this->setBackground(intBkg);
      sign = this->getSignificance(fom::asimov); //AsimovZ is better
      if(docuteff) cuteff->SetBinContent(ind,sign);
      if(max_sign < sign) {max_sign = sign; maxBin = ind;}
    }
  if(info) std::cout << " cut for the max significance is on bin " << maxBin << " with significance " << max_sign << std::endl;
  return max_sign;
}

