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
