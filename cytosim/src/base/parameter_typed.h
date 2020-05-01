//RCS: $Id: parameter_typed.h,v 2.10 2005/04/18 10:22:58 nedelec Exp $
//  F. Nedelec, nedelec@embl-heidelberg.de,  this version Mai 2003


#ifndef PARAMETER_TYPED_H
#define PARAMETER_TYPED_H

#include <cstring>
#include <stdlib.h>
#include "parameter.h"
#include "iomessages.h"

///templated class containing parameter values, and inheriting the descriptions from Parameter
template<class T>
class ParameterTyped : public Parameter
{

  /// a pointer to the array or variable holding the data
  T * ptr;

public:

  //-----------------------------------------------------------------
  /// Creator takes the adress of the data, together with its type
  ParameterTyped(const char * nam, T * p, const int size, const char * def = 0, const char * des = 0, int mod = 0 )
  : Parameter( nam, size, def, des, mod ) {
    ptr = p; 
    clearToDefault();
  }

  //-----------------------------------------------------------------
  /// returns the adress where the data is stored
  char * addr(const int pos = 0)  {
    return (char*)( ptr + pos ); 
  }

  //-----------------------------------------------------------------
  /// returns the adress where the data is stored
  const char * addrConst(const int pos = 0) const {
    return (const char*)( ptr + pos ); 
  }
  
  //-----------------------------------------------------------------
  /// copy the value from another parameter
  void copyValue(const int indx1, const Parameter * p, const int indx2) {
    //we check that the type is the same:
    if ( 0 == strcmp( typeCode(), p->typeCode()))
      ptr[indx1] = *((T*)(p->addrConst(indx2)));
  }

  //-----------------------------------------------------------------
  /// modify the values by applying operator 'op' with real stored in word
  int operate(const char op, const char * word, const int pos ) {
    if ( op == ':' )
      MSG.error("Parameter", "internal error: op == ':'");
    
    if ( pos >= pasize )
      MSG.error("Parameters","operate: wrong argument pos %i", pos);

    //we scan for a real value
    double rhs;
    if ( 0 == sscanf( word, "%lf", &rhs ))
      return 1;
    
    if ( paread < pos )
      paread = pos;
        
    switch( op ) {
      case '+': ptr[pos] = T( ptr[pos] + rhs ); break;
      case '-': ptr[pos] = T( ptr[pos] - rhs ); break;
      case '*': ptr[pos] = T( ptr[pos] * rhs ); break;
      case '/': ptr[pos] = T( ptr[pos] / rhs ); break;
    }
    return 0;
  }

 
  //-----------------------------------------------------------------
  int callScanf( const char * word, const int pos )
  {
    return  sscanf( word, scanCode(), ptr+pos );
  }  

  //-----------------------------------------------------------------
  /// print value[pos] to file
  int printWord( char * line, const int lsize, int pos ) {
    if (( pos >= 0 ) && ( pos < pasize ))
      return snprintf( line, lsize, printCode(), ptr[pos] );
    else {
      MSG.error("Parameter::printWord()", "index out of range");
      return 0;
    }
  }

  //-----------------------------------------------------------------
  /// return the index of the last value with is different from the default
  int lastDifferentFromDefault() {
    if ( pasize <= 0 ) return -1;
    int result = -1;
    T * bak_ptr = new T[ pasize ];
    if ( bak_ptr == 0 ) {
      fprintf(stderr, "ParameterTyped<>::lastDifferentFromDefault() memory allocation failed\n");
      exit(1);
    }
    ParameterTyped<T> bak(paname, bak_ptr, pasize, padeft);
    for( int ii = 0; ii < pasize; ++ii )
      if ( bak_ptr[ ii ] != ptr[ ii ] ) 
        result = ii;
    delete[] bak_ptr;
    return result;
  }

  //-----------------------------------------------------------------
  //these function should be specified for each template:
  
  /// a string describing the type e.g. "int"
  char * typeCode() const;
  
  /// the corresponding code for scanf()
  char * scanCode() const;
  
  /// the printing code for printf()
  char * printCode() const;
  
  /// if=1, warning will be given if '.' found
  bool   shouldBeInteger() const;
  
  /// if=1, issue a warning if a '-' is found
  bool   shouldBePositive() const;
};


#endif
