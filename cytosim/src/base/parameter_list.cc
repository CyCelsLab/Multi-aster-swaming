//RCS: $Id: parameter_list.cc,v 2.22 2005/03/25 17:43:27 nedelec Exp $
//------------------------------parameters.cc---------------------------------

#include <cctype>
#include <cstring>

#include "types.h"
#include "parameter_list.h"
#include "iomessages.h"
#include "exceptions.h"

//---------------------------------------------------------------------------
ParameterList::ParameterList()
{
  nb_parameters = 0; 
  for(int ii = 0; ii < MAX; ++ii )
    defs[ ii ] = 0;
  snprintf( datafile, sizeof(datafile), "no-file" );
}

//---------------------------------------------------------------------------
int ParameterList::clearToDefaults(const char * name) const
{
  Parameter * p = locate( name );
  if ( p ) {
    p->clearToDefault();
    return NO_ERROR;
  }
  return 1;
}

//---------------------------------------------------------------------------
void ParameterList::clearToDefaults() const
{
  for(int ii = 0; ii < nb_parameters; ++ii )
    if ( defs[ ii ] ) 
      defs[ ii ] -> clearToDefault();
}


//---------------------------------------------------------------------------
ParameterList::~ParameterList()
{
  for(int ii = 0; ii < nb_parameters; ++ii )
    if ( defs[ ii ] ) delete defs[ ii ];
}

//---------------------------------------------------------------------------
int ParameterList::link(Parameter * p){ 
  if ( nb_parameters >= MAX ) {
    MSG.error("ParameterList", "Overflow");
    return 1; 
  }
  defs[ nb_parameters++ ] = p;
  return NO_ERROR;
}


//---------------------------------------------------------------------------
ParameterList & ParameterList::operator =(const ParameterList & other)
  //we copy whatever values are common to both.
{
  for(int ii = 0; ii < other.nb_parameters; ++ii ) {
    Parameter * c = other.defs[ii];
    if ( c && c->pasize ) {
      Parameter * p = locate( c->paname );
      if ( p ) {
        for(int jj = 0; ( jj < p->pasize ) && ( jj < c->pasize ); ++jj )
          p->copyValue(jj, c, jj);
        p->paread = c->paread;
      }
    }
  }
  return *this;
}

//---------------------------------------------------------------------------
int ParameterList::sizeOfParameter(const char * name) const
  //maximum number of values the parameter can hold
{
  Parameter * p = locate( name );
  if ( p ) return p->pasize;
  return 0;
}

//---------------------------------------------------------------------------
void ParameterList::showThisManyValues(const char * name, int howmuch) const
  //set the parameter number of read, which is also howmuch will be printed
{
  Parameter * p = locate( name );
  if ( p == 0 )
    throw IOException("parameter not found");

  if ( howmuch < 0 )
    p->paread = p->pasize;
  else
    p->paread = howmuch;
}

//---------------------------------------------------------------------------
int ParameterList::nbValuesRead(const char * name) const
  //how many times parameter 'name' has been read
{
  Parameter * p = locate( name );
  if ( p ) return p->paread;
  return 0;
}


//---------------------------------------------------------------------------
int ParameterList::nbValuesRead() const
{
  int result = 0;
  for (int ii = 0; ii < nb_parameters; ++ii )
    if (( defs[ii] ) && ( defs[ii]->paname )) {
      //printf("%lx %s %i\n", this, defs[ii]->paname, defs[ii]->paread);
      result += defs[ii]->paread + 1;
    }
  return result;
}

//---------------------------------------------------------------------------
Parameter * ParameterList::locate(const char *name) const
{
  //printf("%lx:locate(%s) %i\n", this, name, nb_parameters);
  //the name has to be a valid string
  if ( *name == '\0' )
    return 0;

  //the list has to contain at least one parameter:
  if ( nb_parameters <= 0 )
    throw IOException("empty ParameterList!");
  
  Parameter * result = 0;

  //scan through the list
  for (int ii = 0; ii < nb_parameters; ++ii) {
    if ((defs[ii]) && (defs[ii]->paname)
        && ( strlen(name) == strlen(defs[ii]->paname))
        && ( 0 == strncmp(name, defs[ii]->paname, strlen(name)))) {
      result = defs[ii];
      break;
    }
  }
  

  if ( result == 0 )
    return 0;    //the parameter name was not found!
  
  //if pasize==0, the parameter is deprecated, we search the up-to-date entry,
  //which is the parameter pointing to the same address (parameter values)
  if (( result -> pasize == 0 ) && ( result -> addrConst())) {
    for (int ii = 0; ii < nb_parameters; ++ii ) {
      if (( defs[ii] -> addrConst() == result -> addrConst() ) && ( defs[ii]->pasize )) {
        MSG(1, "Parameter [%s] is deprecated, use [%s] instead", result->paname, defs[ii]->paname );
        result = defs[ii];
        break;
      }
    }
  }
  
  //we did not find a up-to-date entry:
  if ( result->addrConst() == 0 )
    throw IOException("obsolete parameter (suppressed)");

  return result;
}



//---------------------------------------------------------------------------
int ParameterList::parseLineNoCatch(const char * line, const bool verbose) const
{
  int indx = 0;
  char name[256] = "\0";

  //--------------------------- skip initial white space
  while (isspace(line[indx]) && (line[indx] != '\0')) 
    ++indx;
  
  //--------------------------- skip entire line if comments (# or /)
  if ( ! isalpha( line[indx] )) {
    if ((line[indx] == '#' ) || (line[indx] == '/') || (line[indx] == '\0'))
      return NO_ERROR;
    throw IOException("non-alpha first character");
  }

  //--------------------------- get the name as an alphanum sequence:
  unsigned int npos = 0;
  while( isalnum(line[indx] ) && ( npos+1 < sizeof(name) ))
    name[npos++] = line[indx++];

  //terminate the string:
  name[npos] = '\0';
  
  if ( npos+1 >= sizeof(name) )
    throw IOException("the parameter name is too long");

  //find the corresponding parameter:
  Parameter * p = locate( name );
  
  if (( p == 0 ) || ( p->pamodifier & PARAM_INVISIBLE )){
    //printf("%p::parseLineNoCatch [%s] not found, among %i names\n", this, name, nb_parameters);
    return 1;  //throw IOException("unknown parameter!");
  }
  //printf("%p::parseLineNoCatch found name [%s] in line [%s]\n", this, p->paname, line);

  if ( p->pamodifier & PARAM_CONSTANT )
    throw IOException("parameter is declared constant and cannot be modified");

  if (( p->paread >= 0 ) && verbose )
    MSG.warning("Parameter", "[%s] might be multiply defined", name);
  
  //the rest is handled by the parameter:  
  p->parseLine( line+indx );
  return NO_ERROR;
}

//---------------------------------------------------------------------------
int ParameterList::parseLine(const char * line, const bool verbose) const
{
  try {
    return parseLineNoCatch( line, verbose );
  } catch( IOException e ) {
    if ( verbose )
      MSG(1, "Parameter::parseLine %s: ignored [%s]\n", e.getMessage(), line);
    return 1;
  }
}

//----------------------------------------------------------------------------
void ParameterList::parseFile(FILE * file, const bool verbose) const
{
  char line[STRING_SIZE];
  do {
    
    if ( ferror(file) ) {
      if ( verbose ) 
        MSG("ParameterList::ParseFile %s: file error\n");
      return;
    }
    
    if (0 == fgets( line, STRING_SIZE, file )) {
      if ( ferror(file) && verbose )
        MSG("ParameterList::ParseFile %s: file error\n");
      return;
    }
    
    // we check for line[] overflow
    if ( strlen(line) >= STRING_SIZE-1 )
      MSG.error("parseFile","line overflow: you should recompile with a bigger STRING_SIZE");

    //substitute new-line / carriage return and clean up some of the junk:
    for(unsigned int ii = 0; (ii < STRING_SIZE) && (line[ii] != '\0') ; ++ii ) {
      if (( line[ii] == '\n' ) || ( line[ii] == '\r' ) || (! isprint(line[ii])))
        line[ii]=' ';
    }
        
    //printf("%p:parseFile considering line=[%s] %i\n", this, line, nb_parameters);
    try {
      if ( parseLineNoCatch(line, verbose) && verbose )
        MSG(1, "Parameter::parseFile unknown parameter: ignored [%s]\n", line);
    } catch( IOException e ) {
      if ( verbose )
        MSG(1, "Parameter::parseFile %s: ignored [%s]\n", e.getMessage(), line);
    }
  } while( !feof(file) );
}


//---------------------------------------------------------------------------
int ParameterList::parseFile(const char * filename, const bool clearToDefaultBefore, const bool verbose)
{
  //copy the file name:
  char filename_copy[STRING_SIZE];
  snprintf( filename_copy, sizeof(filename_copy), "%s", filename );

  //open the file
  FILE * file = fopen( filename, "r" );
  if ( file == 0 ) {
    if ( verbose )
      MSG(6, "ParameterList: fopen(%s) == 0\n", filename );
    return 1;
  }
  if ( ferror(file) ) {
    fclose(file);
    if ( verbose ) 
      MSG(6, "ParameterList: ferror(%s)\n", filename );
    return 2; 
  }
  
  if ( verbose )
    MSG(6, "ParameterList: reading file [%s]\n", filename );
  
  //clear the parameter values if desired
  if ( clearToDefaultBefore ) clearToDefaults();
  
  //parse the file:
  parseFile(file, verbose);
  
  //store the name of the file in datafile again:
  snprintf( datafile, sizeof(datafile), "%s", filename_copy);
  if ( verbose )
    MSG(6, "ParameterList: reading file [%s] done\n", filename );
  
  //cleanup
  fclose(file);
  return NO_ERROR;
}

//---------------------------------------------------------------------------
int ParameterList::printLine( char * line, const int lsize, ParameterModifier printLevel ) const
{
  int pos = 0;
  for (int ii = 0; ii < nb_parameters; ++ii )
    if ( defs[ii] ) {
      int lpos = defs[ii]->printLine( line+pos, lsize-pos, printLevel );
      if ( lpos > 0 ) {
        pos += lpos;
        pos += snprintf( line+pos, lsize-pos, "\n" );
      }
    }
  return pos;
}

//---------------------------------------------------------------------------
void ParameterList::printFile( FILE * file, ParameterModifier printLevel ) const
{
  char line[STRING_SIZE];
  for (int ii = 0; ii < nb_parameters; ++ii )
    if ( defs[ii] ) {
      int lpos = defs[ii]->printLine( line, sizeof(line), printLevel );
      if ( lpos > 0 )
        fprintf( file, "%s\n", line );
    }
}


//---------------------------------------------------------------------------
int ParameterList::printFile( const char * filename, ParameterModifier printLevel ) const
{
  FILE * file = fopen( filename, "w" );
  if (( file == 0 ) || ( ferror(file)) ) 
    return 1;
  printFile( file, printLevel );
  fclose(file);
  return NO_ERROR;
}

//---------------------------------------------------------------------------
void ParameterList::printDescription(const char * parameterName, FILE * file, const ParameterModifier printLevel ) const
{
  if (( parameterName == 0 ) || ( * parameterName == '\0' )) {
    printDescriptions( file, printLevel );
  } else {
    char line[STRING_SIZE];
    Parameter * p = locate( parameterName );
    if ( p && ( NO_ERROR == p -> printDescription( line, sizeof(line) ))) {
      fprintf(file, "%s\n", line );
    } else {
      fprintf(file, "unknown parameter [%s]\n", parameterName);
    }     
      
  }
}

//---------------------------------------------------------------------------
void ParameterList::printDescriptions( FILE * file, const ParameterModifier printLevel ) const
{
  char line[STRING_SIZE];
  for (int ii = 0; ii < nb_parameters; ++ii ) {
    if ( defs[ii] && ( defs[ii]->pamodifier & printLevel ))
      //check that the printLevel matches the parameter type (bit-wise)
      if ( NO_ERROR == defs[ii]->printDescription( line, sizeof(line) ))
        fprintf(file, "%s\n", line );
  }
}

//---------------------------------------------------------------------------
void ParameterList::printDescriptions( const char * match, FILE * file, const ParameterModifier printLevel ) const
{
  char line[STRING_SIZE];
  for (int ii = 0; ii < nb_parameters; ++ii ) {
    if ( defs[ii] && ( defs[ii]->pamodifier & printLevel ))
      if ( defs[ii]->paname == strstr( defs[ii]->paname, match ))
        if ( NO_ERROR == defs[ii]->printDescription( line, sizeof(line) ))
          fprintf(file, "%s\n", line );
    }
}


//---------------------------------------------------------------------------
void ParameterList::showHelp()
{
  printf("------------------------------ parameters ----------------------------------\n");
  printf(" The parameters below can be set or modified on the command line:\n");
  printf(" to set:       name=value\n");
  printf(" to modify:    name*=value, idem with +=, -=, /=\n");
  printf(" examples:     <play zoom*=4> starts play with a zoom-in of 4\n");
  printf(" for an array like focus, you may use on the command line:\n");
  printf(" \"focus 1 2\"    focus=\"1 2\"   focus=1:2    focus[0]=1\n");
  printf("---------------------------- list of parameters: ---------------------------\n");
}


