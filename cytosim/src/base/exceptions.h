//RCS $Id: exceptions.h,v 2.6 2005/03/23 21:02:28 nedelec Exp $
//Some error conditions are handled by throwing exceptions.
//here we define a very primite Exception class for cytosim

#ifndef EXCEPTIONS_H
#define EXCEPTIONS_H

#include <cstdio>
#include <cstring>

/// Exceptions are a mechanism to handle almost-error conditions. See C++ manual
class Exception {
 
  char * mMessage;
  
public:
    
    /// Creator with empty message
    Exception() {
      mMessage = strdup("\0");
    }
    
    /// Creator with a given message
    Exception(const char * m) {
      mMessage = strdup(m);
    }
    
    /// Prepend something to the message
    void addBeforeMessage(const char * m) {
      const int newSize = strlen(mMessage)+strlen(m);
      char * tmp = new char[newSize];
      snprintf( tmp, newSize, "%s%s", m, mMessage );
      delete [] mMessage;
      mMessage = tmp;
    }
    
    /// Append something to the message
    void addAfterMessage(const char * m) {
      const int newSize = strlen(mMessage)+strlen(m);
      char * tmp = new char[newSize];
      snprintf( tmp, newSize, "%s%s", mMessage, m );
      delete [] mMessage;
      mMessage = tmp;
    }
    
    /// return the message
    const char * getMessage() const { 
      return mMessage; 
    }
    
    /// printf() of the message in the exception
    void print() const {
      printf(" : %s\n", mMessage);
    }

    /// Destructor (exceptions should have an empty destructor)
    virtual ~Exception() {}
};



/// StuckException is thrown if an suitable initial configuration cannot be found
class StuckException : public Exception {
public:
  StuckException(const char * m) : Exception(m) {
    //printf("new StuckException [%s]\n", m);
  }
  virtual ~StuckException() {};
};

/// IOException is thrown during Input/Output of coordinate files
class IOException : public Exception {
public :
  IOException(const char * m) : Exception(m) {
    //printf("new IOException [%s]\n", m);
  }
  virtual ~IOException() {};
};

#endif
