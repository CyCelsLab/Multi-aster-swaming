//RCS: $Id: parameter.h,v 2.10 2005/04/18 10:22:41 nedelec Exp $
//--------------------------------parameters.h-------------------------------
// Classes to read / write / change parameters by names from a file,
// or from the command line
//
//  F. Nedelec, nedelec@embl-heidelberg.de,  this version Mai 2003

#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <cstdio>

///The nature of the parameter and the output options is described by an enum
enum ParameterModifier { 
  PARAM_NORMAL      = 1,  ///<default parameter: modifiable, visible in all lists
  PARAM_ADVANCED    = 2,  ///<advanced parameter are hidden in basic lists
  PARAM_ALL         = PARAM_NORMAL + PARAM_ADVANCED,   ///<to print both normal and advanced parameters
  PARAM_INVISIBLE   = 4,  ///<to hide a parameter from the user
  PARAM_CONSTANT    = 8,  ///<a parameter can be made constant, to prevent user modification
  PARAM_ALL_VALUES  = 16, ///<to specify that all values should be written
  PARAM_NOT_DEFAULT = 32  ///<to print only values different from default
};

//---------------------------------------------------------------------------
//   class Parameter is the base class for all parameters
//---------------------------------------------------------------------------

///parent class for parameters, contains the descriptions but not the values, as they are type-dependent
class Parameter
{
  friend class ParameterList;
protected:
    
  /// a string with the name by which the parameter is refered
  char * paname;
    
  /// the default value, expressed as a string, from which value is scanned
  char * padeft;
  
  /// an optional short description of the parameters function
  char * padesc;
  
  /// array size; size=0 is used for a deprecated name
  int    pasize;
  
  /// the highest index of a value that was actually set for this parameter
  int    paread;
  
  /// if=1, the parameter cannot be changed, if=2, the parameter is hidden
  int    pamodifier;

public:

  /// Constructor 
  Parameter(const char * nam, const int siz, const char * def, const char * des, const int modifier );
  
  /// Destructor  
  virtual ~Parameter();
  
  /// a short string describing the type
  virtual char * typeCode()                            const  = 0;
  
  /// the addreses holding  value[ pos ]
  virtual char * addr(const int pos = 0)                      = 0;

  /// the addreses holding  value[ pos ]
  virtual const char * addrConst(const int pos = 0)     const = 0;

  /// clear the parameter to zero
  void clearToDefault();
  
  /// copy the values from another parameter
  virtual void copyValue( const int, const Parameter *, const int ) = 0;
  
  /// apply standard operations to the values
  virtual int operate( const char, const char *, const int ) = 0;
  
  /// print value[pos] to the file
  virtual int printWord( char * line, const int lsize, const int pos )             = 0;
  
  /// read value[pos] from the string
  virtual void parseWord( const char *, const int pos );
  
  /// call scanf( word, scanCode(), ptr+pos );
  virtual int callScanf( const char *, const int pos )        = 0;
  
  /// index of the last value with is different from the default
  virtual int  lastDifferentFromDefault()                     = 0;

  /// read the value[pos] and beyond from the string
  void parseLine( const char * );
  
  /// print the values to the file
  int printLine( char *, const int lsize, ParameterModifier printLevel );
  
  /// print a short description to the file, or an error-code
  int printDescription( char *, const int lsize );
  
  /// if=1, warning will be given if '.' found
  virtual bool   shouldBeInteger()  const { return 0; }
  
  /// if=1, issue a warning if a '-' is found
  virtual bool   shouldBePositive() const { return 0; }
  
  ///true if the character signals the end of a line, or beggining of comments
  bool terminatingCharacter(const char ch) const;
};

#endif
