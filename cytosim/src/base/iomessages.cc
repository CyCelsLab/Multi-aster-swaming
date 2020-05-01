//RCS: $Id: iomessages.cc,v 2.6 2005/03/16 14:53:39 nedelec Exp $

#include "iomessages.h"

//the global instantiation IOMessages used for input/output
IOMessages MSG;

#include <cstdlib>

//---------------------------------------------------------------------------
IOMessages::IOMessages()
{
    myout = stdout;
    
    //default values:
    ABORT_ON_ERROR = true;
    VERBOSE        = 4;
}

//---------------------------------------------------------------------------
IOMessages::~IOMessages()
{
    closeFile();
}


//---------------------------------------------------------------------------
void IOMessages::setFile( FILE * f )
{
    closeFile();
    if ( f && !ferror(f) ) 
        myout = f;
}

//---------------------------------------------------------------------------
FILE * IOMessages::openFile(const char * name)
{
    setFile( fopen(name, "w") );
    return myout;
}

//---------------------------------------------------------------------------
void IOMessages::closeFile()
{
    if ( myout && ( myout != stdout ))
        fclose( myout );
    myout = stdout;
}

//---------------------------------------------------------------------------
void IOMessages::vprint(const char * fmt, va_list & args)
{
  if ( VERBOSE >= 0 ) {
    vfprintf(myout, fmt, args);
    fflush(myout);
  }
}

//---------------------------------------------------------------------------
void IOMessages::print(const char * fmt, ...)
{
  if ( VERBOSE >= 0 ) {
    va_list args;
    va_start(args, fmt);
    vfprintf(myout, fmt, args);
    va_end(args);
    fflush(myout);
  }
}

//---------------------------------------------------------------------------
void IOMessages::print(const int level, const char * fmt, ...)
{
  if ( VERBOSE >= level ) {
    va_list args;
    va_start(args, fmt);
    vfprintf(myout, fmt, args);
    va_end(args);
    fflush(myout);
  }
}

//---------------------------------------------------------------------------
//the operators below are direct copies of print() above:
void IOMessages::operator()(const char * fmt, ...)
{
  if (VERBOSE >= 0) {
    va_list args;
    va_start(args, fmt);
    vfprintf(myout, fmt, args);  
    va_end(args);
    fflush(myout);
  }
}


//---------------------------------------------------------------------------
void IOMessages::operator()(const int level, const char * fmt, ...)
{
  if ( VERBOSE >= level ) {
    va_list args;
    va_start(args, fmt);
    vfprintf(myout, fmt, args);
    va_end(args);
    fflush(myout);
  }
}



//---------------------------------------------------------------------------
void IOMessages::error(const char * funct, const char * fmt, ...)
{
    static int errorCnt = 0;
    if (errorCnt > USER_LIMIT)
      return;
    
    va_list args;
    //print the error to stderr:
    fprintf(stderr, "ERROR in %s: ", funct);
    va_start(args, fmt);
    vfprintf(stderr, fmt, args);
    va_end(args);
    fprintf(stderr, "\n");
    
    if (( myout != stderr ) && ( myout != stdout )) {
      print("ERROR in %s: ", funct);
      va_start(args, fmt);
      vprint(fmt, args);
      va_end(args);
      print("\n");      
    }
    
    if (ABORT_ON_ERROR) 
      abort();
    else
      exit(EXIT_FAILURE);
    
    if (++errorCnt > USER_LIMIT) 
        print("error messages are now silent\n");
}


//---------------------------------------------------------------------------
void IOMessages::warning(const char * funct, const char * fmt, ...)
{
  static int WARNING_CNT = 0;
  
  if (( VERBOSE >= 0 ) && ( WARNING_CNT < USER_LIMIT )) {
    va_list args;
    print("warning : %s : ", funct);
    va_start(args, fmt);
    vfprintf(myout, fmt, args);
    va_end(args);
    fprintf(myout, "\n");
    fflush(myout);

    if (++WARNING_CNT >= USER_LIMIT) 
        fprintf(myout, "warning messages are now silent\n");
  }
}

//---------------------------------------------------------------------------
void IOMessages::warning(const int level, const char * funct, const char * fmt, ...)
{
  static int WARNING_CNT = 0;
  
  if (( VERBOSE >= level ) && ( WARNING_CNT < USER_LIMIT )) {
    va_list args;
    print("warning : %s : ", funct);
    va_start(args, fmt);
    vfprintf(myout, fmt, args);
    va_end(args);
    fprintf(myout, "\n");
    fflush(myout);
    
    if (++WARNING_CNT >= USER_LIMIT) 
      fprintf(myout, "warning messages are now silent\n");
  }
}


