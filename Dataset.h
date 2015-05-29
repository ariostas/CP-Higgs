
#ifndef Dataset_h
#define Dataset_h

class Dataset {
  
protected:
  vector<vector<double > > _data;
  unsigned int _ndim;

public:
  
  Dataset():_ndim(0){}

  Dataset(int ndim, int maxevents): _ndim(ndim){
    _data.reserve(maxevents);
  }

  virtual ~Dataset(){}

  void add(vector<double> &event){
    if(event.size() != _ndim) {
      cout << "WTF are you doing adding an event with dimension " 
	   << event.size() << "?!?" << endl;
    }
    _data.push_back(event);
  }
  
  void setNDim(int n){ _ndim=n;}
  void reserve(int n){ _data.reserve(n);}
  int size(){return _data.size();}
  int ndim(){return _ndim;}
  double get(int ev, int d){return _data[ev][d];}

  double distance(int i, const Dataset &data, int j){
    double d2=0;
    for(int n=0; n < _ndim; n++){
      d2 += pow(_data[i][n]-data._data[j][n],2);
    }
    return sqrt(d2);
  }

};

#endif


