//RCS: $Id: parameter_list.h,v 2.15 2005/04/18 10:22:52 nedelec Exp $

#ifndef PARAMETERS_LIST_H
#define PARAMETERS_LIST_H

// TODO: implement simple operation to calculate derived parameter
// that would accept:
// rule( "ratedt", &ratedt, 1, "rate * dt" );
// check( "ratedt > 0" );
// this would make assertions of parameter values more systematic
// and get rid of computeDerivedValues() in class ParamSim



//---------------------------------------------------------------------------
// ParameterList hold a static list of Parameter definitions, and can be used
// to read/write the corresponding parameter values from a file 
//---------------------------------------------------------------------------
//  F. Nedelec, nedelec@embl-heidelberg.de,  this version Mai 2003

/*
 Example:

 class ParamTest : public ParameterList
 {
public:
   
   int     an_int;
   char    string[256];
   int     array_of_int[6];
   double  a_double;
   real    a_real;
   
   ParamTest() {
     link( "name",      &an_int,                  "1");
     link( "string",    a_string,        256,     "Cytosim");
     link( "array",     array_of_int,      6,     "0 1 2 3");
     link( "areal",     &a_real,                  "3.14");
     
     //deprecated name are mapped to the given variable with a warning:
     link_deprecated( "oldname", &an_int);
   }
   
   virtual ~ParamTest() { }
 };
 
 
 ParamTest test;
 
 int main(int argc, char * argv[])
 {
   test.printDescriptions();
   test.parseFile("data.test");
   if ( argc > 1 ) test.parseLine(argv[1]);
   test.print(stdout);
   return EXIT_SUCCESS;
 }
 */

#include "types.h"
#include "parameter_typed.h"

/// Contains a lits of ParameterTyped
class ParameterList
{

public:
  //would be better if we could make datafile private
  ///name of the file from which the parameters were set:
  char        datafile[STRING_SIZE];
  
private:
  //we use a simple C array, with fixed limits for simplicity: ///TODO: use automatic array
  
  ///upper limit to the number of parameters
  static const int MAX = 228; 
  
  ///the current number of parameters
  int         nb_parameters;
  
  ///the array holding the definitions
  Parameter * defs[MAX];

  ///returns one of the defs[*] corresponding to 'name; or zero if not found
  Parameter * locate(const char * name) const;

  //add a new parameter to the list:
  int link(Parameter * p);
  
 public:

  ///Creator
  ParameterList();
  
  ///Destructor
  virtual ~ParameterList();

  /// Copy all the value from another list
  ParameterList & operator =(const ParameterList &);

  /// Add a array parameter to the list
  template<class T> 
    int link(const char * name, T * ptr, const int size, const char * def, const char * desc=0, int mod=PARAM_NORMAL ) { 
      return link( new ParameterTyped<T>( name, ptr, size, def, desc, mod ));
    }

  /// Add a parameter to the list, one value only
  template<class T> 
    int link(const char * name, T * ptr, const char * def, const char * desc=0, int mod=PARAM_NORMAL ) {
      return link( new ParameterTyped<T>( name, ptr, 1, def, desc, mod ));
    }

  /// Add a deprecated parameter to the list
  template<class T> 
    int link_deprecated(const char * name, T * ptr = 0) {
      return link( new ParameterTyped<T>( name, ptr, 0, 0, 0 ));
    }
  
  /// Add an invisible modifiable parameter to the list
  template<class T> 
    int link_invisible(const char * name, T * ptr = 0) {
      return link( new ParameterTyped<T>( name, ptr, 0, 0, 0, PARAM_NORMAL+PARAM_INVISIBLE ));
    }
  
  /// clear all values to their default
  int clearToDefaults(const char * parameterName) const;

  /// clear all values to their default
  void clearToDefaults() const;

  /// read a parameter and its value from a string, returns NO_ERROR or an error code
  int  parseLineNoCatch(const char *, const bool verbose = true) const;

  /// read a parameter and its value from a string, returns NO_ERROR or an error code
  int  parseLine(const char *, const bool verbose = true) const;
  
  /// real all parameters and their values form a file, this catches all exceptions
  void parseFile(FILE *, const bool verbose = true) const;
  
  /// read all parameters and their values from the file with name filename, resetting to default before
  /** return NO_ERROR if everythings OK */
  int  parseFile(const char * filename, const bool clearToDefaultBefore = true, const bool verbose = true);

  /// returns the name of the last file parsed
  const char * lastFileParsed() const { return datafile; }
  
  /// print all parameters and their values to the given string
  int  printLine(char *, const int lsize, ParameterModifier printLevel = PARAM_NORMAL) const;

  /// print all parameters and their values to the given file
  void printFile(FILE * = stdout, ParameterModifier printLevel = PARAM_NORMAL) const;
  
  /// print all parameters and their values to the file with name filename
  /** return NO_ERROR if everythings OK */
  int  printFile(const char * filename, ParameterModifier printLevel = PARAM_NORMAL) const;

  /// number of value the array can hold
  int  sizeOfParameter(const char * name) const; 
  
  /// this is non-zero if values have been set through parsing of file or line
  int  nbValuesRead(const char * name) const; 
  
  /// -1 sets to size, i.e. shows all
  void showThisManyValues(const char * name, int howmuch = -1) const;
  
  /// the parameter will not printed
  void hide(const char * name) const { showThisManyValues( name, 0); } 

  /// total number of values set by last parseFile()
  int  nbValuesRead() const; 

  /// print all the descriptions for a given parameters
  void printDescription(const char * name, FILE * = stdout, const ParameterModifier printLevel = PARAM_NORMAL) const;
  
  /// print all the descriptions for all parameters
  void printDescriptions(FILE * = stdout, const ParameterModifier printLevel = PARAM_NORMAL) const;

  /// print the descriptions of parameters which names start with the given string 'match'
  void printDescriptions(const char * match, FILE * = stdout, const ParameterModifier printLevel = PARAM_NORMAL) const;

  /// print a little help
  void showHelp();
};









#endif
