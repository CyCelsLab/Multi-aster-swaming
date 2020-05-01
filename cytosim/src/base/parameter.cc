//RCS: $Id: parameter.cc,v 2.20 2005/04/18 10:22:35 nedelec Exp $

#include <cstdlib>
#include <cctype>
#include <cstring>
#include "types.h"
#include "parameter.h"
#include "iomessages.h"
#include "exceptions.h"

//----------------------------------------------------------------------------
bool Parameter::terminatingCharacter( const char c ) const
{
  //we cannot include C++ style comments '/', since that is used in file-path
  return ((c == '\n') || (c == '\0') || (c == '\r') || (c == '%') || (c == '#'));
}

//----------------------------------------------------------------------------
Parameter::Parameter(const char * nam, const int siz, const char * def, const char * des, const int mod )
{
  pasize     = siz;
  pamodifier = mod;
  paread     = -1;
  
  if ( nam ) {
    paname = strdup( nam );  //strdup() allocates, and a copy the string
    if ( paname == 0 ) {
      fprintf(stderr, "Parameter::Parameter() memory allocation failed\n");
      exit(1);
    }
  } else paname = 0;

  if ( def ) {
    padeft = strdup( def );
    if ( padeft == 0 ) {
      fprintf(stderr, "Parameter::Parameter() memory allocation failed\n");
      exit(1);
    }
  } else padeft = 0;

  if ( des ) {
    padesc = strdup( des );
    if ( padesc == 0 ) {
      fprintf(stderr, "Parameter::Parameter() memory allocation failed\n");
      exit(1);
    }
  } else padesc = 0;
}


//----------------------------------------------------------------------------
Parameter::~Parameter() {
  delete[] paname;
  if ( padesc ) delete[] padesc;
  if ( padeft ) delete[] padeft;
}

//-----------------------------------------------------------------
/// clear all values to the default
void Parameter::clearToDefault() {
  if ( pasize <= 0 ) return;
  
  try {
    parseLine( padeft );
  } catch( IOException e ) {
    MSG.error("Parameters","could not set [%s] from default value [%s]: %s", paname, padeft, e.getMessage());
  }
	
  if ( paread < 0 )
    MSG.error("Parameters","could not set [%s] from default value [%s]", paname, padeft);
  
  //pad the rest of the values with the last default
  for( int jj = paread+1; jj < pasize; ++jj )
    copyValue( jj, this, jj-1 );
  
  //reset the read count
  paread = -1;
}

//-----------------------------------------------------------------
/// read value[pos] from the string 'word'
void Parameter::parseWord( const char * word, const int pos ) {
  
  if ( pos < 0 )       throw IOException("array index out of range (<0)");
  if ( pos >= pasize ) throw IOException("array index out of range");
  
  if ( shouldBeInteger() && strchr( word, '.' ))
    MSG.warning("Parameter", "Ignored decimal part specified for integer [%s].", paname);
  
  if ( shouldBePositive() && strchr( word, '-' ))
    throw IOException("positive parameter, but [-] found");
  
  if ( 1 == callScanf( word, pos )) {
    //printf("Parameter: %15s[%i] : [%s] ", paname, pos, word ); printWord(stdout, pos); printf("\n");
    if ( paread < pos )
      paread = pos;
  } else {
    throw IOException("scanf() failed");
  }
}

//----------------------------------------------------------------------------
void Parameter::parseLine(const char * line)
{
  int indx = 0;
  //const int lsize = strlen(line);
  
  //printf("Parameter::%s: parseLine [%s]\n", paname, line);
  //skip initial space
  while (isspace(line[indx]) && !terminatingCharacter(line[indx]))
    ++indx;  
  
  char op = ':';
  int pos = 0;
  //--------------------------- optional array index
  if ( line[indx] == '[' ) {
    ++indx;
    if (1 != sscanf( line+indx, "%i", &pos ))
      throw IOException("missing/wrong array index");
    if (( pos < 0 ) || ( pos >= pasize ))
      throw IOException("array index out of range");
    
    //find and pass the closing bracket:
    const char * closing = strchr( line+indx, ']' );   // nk: on 15 july
    if ( closing == 0 )
      throw IOException("unbalanced bracket");
    indx = closing - line + 1;
  }
  
  //--------------------------- optional assignment operator
  if ((line[indx+1] == '=') &&  ((line[indx] == '+') || (line[indx] == '-') 
                                 || (line[indx] == '*') || (line[indx] == '/')))
    op = line[indx++];
  
  if (line[indx] == '=') ++indx;
  
  //--------------------------- skip separator (might be space)
  while(isspace(line[indx]) && !terminatingCharacter(line[indx]))
    ++indx;
  
  //string are a special case: we just copy the all line from there
  if ( 0 == strcmp( typeCode(), "char" )) {
    paread = 0;
    while ((paread < pasize-1) && !terminatingCharacter(line[indx]))
      *addr(paread++) = line[indx++];
    *addr(paread) = '\0';
    //printf("Parameter: set string %s = [%s] line = [%s]\n", paname, addr(), line );
    return;
  }
    
  if (line[indx] == '\0')
    throw IOException("unexpected end of line");
  
  //--------------------------- we should now be in the values    
  //we cut the line in words:
  char word[STRING_SIZE];

  while ( !terminatingCharacter(line[indx]) ) {
    //skip space
    while (!terminatingCharacter(line[indx]) && isspace(line[indx]))
      ++indx;

    //we stop if comments / can start with #, % or /: 
    if (terminatingCharacter(line[indx]))
      return;

    //normally, word is as big as line, so sscanf should not overflow
    if ( 1 != sscanf( line+indx, "%s", word )) 
      return;
    //printf("Parameter::%s: parseLine word [%s]\n", paname, word);
    
    indx += strlen( word );
      
    if ( op == ':' )
      parseWord( word, pos );
    else {
      if ( NO_ERROR != operate( op, word, pos ))
        throw IOException("unable to read the right-hand-side of the operator");
    }
    
    ++pos;
  }
}

//----------------------------------------------------------------------------
int Parameter::printLine( char * line, const int lsize, ParameterModifier printLevel )
{
  *line = '\0';
  if ( pamodifier & PARAM_INVISIBLE )
    return 0;
  
  //we treat strings as a special case, because of spaces:
  if ( 0 == strcmp( typeCode(), "char" )) {
    bool doprint = ( paread >= 0 );
    //level 1 prints all the values
    if ( printLevel & PARAM_ALL_VALUES )
      doprint = true;
    if (( printLevel & PARAM_NOT_DEFAULT ) && strcmp( addrConst(), padeft ))
      doprint = true;
    if ( doprint && ( *addrConst() != '\0' ))
      return snprintf( line, lsize, "%-16s %s", paname, addrConst());
    else
      return 0;
  }
    
  //by default, we print all the values set by reading
  int last = paread + 1;
  
  //level 1 prints all the values
  if ( printLevel & PARAM_ALL_VALUES ) 
    last = pasize;
  
  //level 2 prints only values different from the default setting
  if ( printLevel & PARAM_NOT_DEFAULT ) {
    int d = lastDifferentFromDefault() + 1;
    if ( last < d )
      last = d;
  }
    
  if ( last > pasize )
    last = pasize;

  if ( last > 0 ) {
    int pos = snprintf( line, lsize, "%-16s", paname );
    for( int s = 0; s < last; ++s )
      pos += printWord( line+pos, lsize-pos, s );
    return pos;
  }
  return 0;
}

//----------------------------------------------------------------------------
int Parameter::printDescription( char * line, const int lsize )
{
  *line = '\0';
  if ( pasize <= 0 )
    return 1;

  int pos = snprintf( line, lsize, "%-8s", typeCode() );

  if ( pasize > 1 ) {
    char namesize[128];
    snprintf( namesize, sizeof(namesize), "%s[%i]", paname, pasize);
    pos += snprintf( line+pos, lsize-pos, " %-18s", namesize );
  } else {
    pos += snprintf( line+pos, lsize-pos, " %-18s", paname );
  }
  
  if ( padesc ) 
    pos += snprintf( line+pos, lsize-pos, "  %-60s", padesc ); 
  else 
    pos += snprintf( line+pos, lsize-pos, "  %40s", "default" );

  if ( padeft )
    pos += snprintf( line+pos, lsize-pos, " = %s", padeft ); 
  else
    pos += snprintf( line+pos, lsize-pos, " %s", "(null)" );
  
  return NO_ERROR;
}

