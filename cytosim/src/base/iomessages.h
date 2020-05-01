//RCS: $Id: iomessages.h,v 2.4 2005/04/12 07:23:16 nedelec Exp $

#ifndef  IOMESSAGES_H
#define  IOMESSAGES_H

#include <cstdio>
#include <cstdarg>
#include <cctype>
#include <cstdlib>

/*
 \TODO: clean-up the verbose levels, introducing an enum for 
classifying different types of messages; I/O, initial state, 
math, Hand, Microtub, Solid, etc.
we could select only one level of message for example.
 */

/**three convenient functions to output to user:
print() or equivalently operator()
error()
warning()
*/

///A wrapper around FILE and printf() for Cytosim's output
class IOMessages
{
        
private:
    ///the output file
    FILE * myout;
    
    ///verbose level, PRINT has a option specifying a level which is compared to VERBOSE
    int VERBOSE;
    
    ///switch to call abort in ERROR
    bool ABORT_ON_ERROR;
    
    ///max. number of output-warnings
    static const int USER_LIMIT = 50;  
    
    ///writes messages to output file, variable number of arguments
    void    vprint(const char * fmt, va_list & args);

public:
    ///Constructor
    IOMessages();
	
    ///Destructor 
    virtual ~IOMessages();
    
    ///Opens a file
    FILE *  openFile(const char * name);
	
    ///Sets the output file to f
    void    setFile(FILE * f);
	
    ///Closes the output file
    void    closeFile();

    //functions to set the verbose level:
    ///suppresses output by setting Verbose to 0
    void    shutUp()                    { VERBOSE = 0; }
	
    ///Returns the verbose level
    int     getVerboseLevel()           { return VERBOSE; }
    
    ///Sets verbose to level m
    void    setVerboseLevel(int m)      { VERBOSE = m; }
    
    ///Sets verbose to level m, given in ascii
    void    setVerboseLevel(const char * m) {
      if ( isdigit( m[0] )) VERBOSE = (int)strtol(m, (char **)NULL, 10);
    }
	
    ///increases Verbose level by 1
    void    increaseVerboseLevel()      { ++VERBOSE; }
    
    ///decreases Verbose level by 1
    void    decreaseVerboseLevel()      { if ( VERBOSE > 0 ) --VERBOSE; }
    
    ///functions that print error warnings. All use the fprintf() format
    void    print(const char * fmt, ...);
    void    print(const int, const char * fmt, ...);
    
    ///convenient access to the same print() via by overriding the () operator
    void    operator()(const char * fmt, ...);
    void    operator()(const int, const char * fmt, ...);
    
    ///error normally stops the execution, calling abort() its first argument is the name of the calling function
    void    error(const char * funct, const char * fmt, ...);
    void    setAbortOnError(bool m)    { ABORT_ON_ERROR = m; }

    ///warning() issues a warning message, also with a verbose level
    void    warning(const char * funct, const char * fmt, ...);
    void    warning(const int, const char * funct, const char * fmt, ...);
};

///the global instantiation used for Cytosim's output
extern IOMessages MSG;

#endif
