#ifndef libsvm_container_h
#define libsvm_container_h
#include"libsvm-weights-3.18/svm.h"
#include <vector>
class svm_container{
public:
  svm_container(int nVar,int nEntries):nvar(nVar),nentries(nEntries) 
  {
    svm_cont = new std::vector<svm_node_array>[nvar*nentries*sizeof(svm_node)];
    //    weights = new std::vector<double>[nentries];
  };
  typedef std::vector<svm_node> svm_node_array;
  std::vector<svm_node_array>* svm_cont;
  std::vector<double>* weights;
  ~svm_container(){};//{delete svm_cont;};  
private:
  int nvar, nentries; 
  svm_container(){};

};

#endif
