#include <Rcpp.h>
using namespace Rcpp;

#include "BCI2000FileReader.h"

SEXP paramListToSEXP(const ParamList &list);
SEXP paramToSEXP(const Param &list);

// [[Rcpp::export]]
Rcpp::List load_bcidat(std::string file, bool raw=false)
{
  BCI2000FileReader reader;
  reader.Open(file.c_str());
  if(!reader.IsOpen())
  {
    reader.Open((file+".dat").c_str());
    if(!reader.IsOpen())
    return Rcpp::List();
  }
  int samples = reader.NumSamples();
  int channels = reader.SignalProperties().Channels();
  
  //first - read all samples
  Rcpp::NumericMatrix signal(samples,channels);
  
  for(int s=0; s<samples; ++s)
  {
    for(int c=0; c<channels; ++c)
    {
      signal(s,c) = raw?reader.RawValue(c, s):reader.CalibratedValue(c, s);
    }
  }
  
  //read all states
  int numStates = reader.States()->Size();
  //  read all names
  
  Rcpp::NumericMatrix states(samples, numStates);
  
  //  read values
  for( int sample = 0; sample < samples; ++sample )
  {
    reader.ReadStateVector( sample );
    for( int i = 0; i < numStates; ++i )
    {
      State::ValueType value  = reader.StateVector()->StateValue( (*reader.States())[i].Location(), (*reader.States())[i].Length() );
      states(sample, i) = value;
    }
  }
  
  Rcpp::CharacterVector stateNames(numStates);
  for(int j=0; j< numStates; ++j)
    stateNames[j] = (*reader.States())[j].Name();
    
  Rcpp::List dimnms = Rcpp::List::create(R_NilValue, stateNames);

    
  states.attr("dimnames") = dimnms;
  
  //read parameters
  SEXP params = paramListToSEXP(*reader.Parameters());
  
  return Rcpp::List::create(Rcpp::Named("signal") = signal,
                            Rcpp::Named("states") = states,
                            Rcpp::Named("parameters") = params
                            );
}

SEXP paramListToSEXP(const ParamList &list)
{
  Rcpp::List params;
  for(int i=0; i<list.Size(); ++i)
  {
    const Param &param = list[i];
    params[param.Name()] = paramToSEXP(param);
  }
  
  return params;
}

SEXP paramToSEXP(const Param &param)
{
  if(param.NumRows() == 1 &&
    param.NumColumns() == 1)
  {
    return Rcpp::CharacterVector(param.Value().ToString());
  }
  else
  {
    bool isNested = false;
    for( int col = 0; !isNested && col < param.NumColumns(); ++col )
      for( int row = 0; !isNested && row < param.NumRows(); ++row )
        if( param.Value(row, col).Kind() != Param::ParamValue::Single )
          isNested = true;

    if(isNested)
    {
      // list of lists
      Rcpp::List out;
      for(int row=0; row < param.NumRows(); ++row)
      {
        Rcpp::List list;
        for(int col=0; col < param.NumColumns(); ++col)
        {
          list[col] = paramToSEXP(*param.Value(row, col).ToParam());
        }
        out[row] = list;
      }
      
      return out;
    }
    else
    {
      // simple matrix
      Rcpp::CharacterMatrix mat(param.NumRows(), param.NumColumns());
      for(int row=0; row < param.NumRows(); ++row)
      {
        for(int col=0; col < param.NumColumns(); ++col)
        {
          mat(row,col) = param.Value(row, col).ToString();
        }
      }
      return mat;
    }
  }
}
