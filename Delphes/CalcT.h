#include <TRandom3.h>

double wtFcn(int nevts, double d){
  double eps = 0.01/(double)nevts;
  return -log(d+eps);
}

double calcIntegral(Dataset &d1, Dataset &d2, bool same, int index=-1){
  double integral = 0;
  int n1 = d1.size(), n2 = d2.size();
  int ind_min = 0, ind_max = n1;
  if(index >= 0){
    ind_min = index; ind_max = index+1;
  }
  for(int ev1 = ind_min; ev1 < ind_max; ev1++){
    int min = 0;
    if(same) min = ev1+1;
    for(int ev2 = min; ev2 < n2; ev2++){
      double dist = d1.distance(ev1,d2,ev2);
      integral += wtFcn(n1+n2,dist);
    }
  }
  double num_terms = n1*n2;
  if(same) num_terms = n1*(n1-1);
  return integral/num_terms;
}

double calcT(Dataset &d1, Dataset &d2){
  double i11 = calcIntegral(d1,d1,true);
  double i22 = calcIntegral(d2,d2,true);
  double i12 = calcIntegral(d1,d2,false);
  return i11+i22-i12;
}

double permCalcT(Dataset &d1, Dataset &d2, TRandom3 *randGen){
  int ndim = d1.ndim();
  int n1 = d1.size(), n2 = d2.size(), ntot = n1+n2;
  Dataset dd1(ndim,n1),dd2(ndim,n2);
  // re-sample 
  vector<double> this_event(ndim);
  for(int i=0; i<ntot; i++){    
    if(i < n1){
      for(int d=0; d<ndim; d++) this_event[d] = d1.get(i,d);
    }
    else {
      for(int d=0; d<ndim; d++) this_event[d] = d2.get(i-n1,d);
    }
    double r = randGen->Rndm();
    if(dd1.size() >= n1) dd2.add(this_event);
    else if(dd2.size() >= n2) dd1.add(this_event);
    else if(r < n1/(double)ntot) dd1.add(this_event);
    else dd2.add(this_event);
  }

  // get T
  double tval = calcT(dd1,dd2);
  return tval;
}
