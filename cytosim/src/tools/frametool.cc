//RCS: $Id: frametool.cc,v 2.18 2005/04/22 19:30:57 nedelec Exp $
//this is an utility to quickly read and manipulate frames in 'result.out',
//frametool only uses the START and END Tags of frames, and does not interpret
//or verify the logical content of the file.
//for these task, use tools/readwrite or test/testreader

#include <cstdio>
#include <cctype>
#include <cstring>
#include <cstdlib>
#include "types.h"
#include "iowrapper.h"

const int LINEWIDTH = 256*1024;
const int MAX_FRAME = 16384;


void ERROR(char * message)
{
  fprintf(stderr, "ERROR: %s\n", message);
  exit( EXIT_FAILURE );
}

//returns a code to specify which type of line was recognized
int getLine(FILE * infile, char * result)
{
  fgets(result, LINEWIDTH, infile);
  if ( feof(infile) ) return 3;
  if ( result == strstr(result, "#frm " ))        return 1;
  if ( result == strstr(result, "#end " ))        return 2;
  if ( result == strstr(result, " #end " ))       return 2;
  if ( result == strstr(result, FRAME_START_TAG)) return 1;
  if ( result == strstr(result, FRAME_END_TAG))   return 2;
  if ( result == strstr(result, "#F " ))          return 1;
  return 0;
}


//=============================================================================
void countFrame(FILE * infile)
{
  char line[LINEWIDTH];
  int frame  = 0, linecnt = 0;
  long start = 0, end, size=0, size_old;
  int code;

  do {
    ++linecnt;
    code = getLine(infile, line);
    
    if (code && ( linecnt > start+8 )) {
      
      end = linecnt;
      size_old = size;
      size = end - start;
      printf(" (lines %7li  %+7li)\n", size, size - size_old);
      start = linecnt + 1;
    
    }
    
    if ( code == 1 ) {
      
      start = linecnt;
      printf("F %3d :", frame++);
      line[ strlen(line)-1 ]='\0';
      printf("[%s]", line);
      
    }
    
  } while ( !feof(infile) );
}


//=============================================================================
void extractFile(FILE * infile, const unsigned long start, const unsigned long finish)
{
  clearerr(infile);
  fseek(infile, start, SEEK_SET);
  for( unsigned long T = start; T < finish; ++T ) 
    putchar( fgetc( infile ) );
  fseek(infile, finish, SEEK_SET);
  printf("\n\n\n");
}

//=============================================================================
void extractFrame(FILE * infile, const char list[])
{
  char line[LINEWIDTH];
  unsigned long lpos, start = 0;
  int  code, frame = -1;
  bool inside = false;
  
  do {
    lpos = ftell(infile);
    code = getLine(infile, line);

    switch( code ) {
      case 0: break;
      case 1: {
        if ( inside ) {
          if ( 'h' != list[ frame ] ) 
            extractFile(infile, start, ftell(infile));
          if ( 'e' == list[ frame ] )
            return;
        }
        inside = true;
        start  = lpos;
        ++frame;
        //fprintf(stderr, "frame %i : %i\n", frame, start);
        if ( frame >= MAX_FRAME ) 
          return;
      } break;

      case 2: {
        inside = false;
        if ( 'h' != list[ frame ] ) 
          extractFile(infile, start, ftell(infile));
        if ( 'e' == list[ frame ] ) 
          return;          
      } break;
    }
    
  } while ( !feof(infile) );
}


//=============================================================================
void extractLastFrame(FILE * infile)
{
  char line[LINEWIDTH];
  long pos, start;

  do {
    pos = ftell( infile );
    if ( getLine(infile, line) == 1 ) 
      start = pos;
  } while ( !feof(infile) );

  clearerr(infile);
  fseek(infile, start, SEEK_SET);

  do {
    putchar( fgetc( infile ) );
  } while (!feof(infile));
  printf("\n");
}


//=============================================================================
void showHelp()
{
  printf("Usage: frametool -option[nb nb2...] filename\n");
  printf("possible options:\n");
  printf(" -c   Counting\n");
  printf(" -e   Extracts\n");
  printf("        ex.: -e10  -e10,13  -e10-13  -e10-  -e-13\n");
  printf(" -l   Extract last frame\n");
  printf(" -p   Extract periodic\n");
  printf("        ex.: -p2  -p2,1\n");
  printf(" -r   Remove, same syntax as -e\n");
}



//=============================================================================
int main(int argc, char* argv[])
{
  const char * filename = RESULT_OUT;
  FILE * infile;
  
  //a list of frame to collect:
  char list[MAX_FRAME];
  for( int ii = 0; ii < MAX_FRAME; ++ii )
    list[ii] = 'h';
  
  char action='c';

  
  //get the list of frame from the command line arguments:
  for( int ar = 1; ar < argc; ++ar ) {
    
    if (argv[ar][0]=='-') {
      
      switch (argv[ar][1]) {
        case 0:
	    case 'h':
	      showHelp();
	      return 0;
          
	    case 'c':
	      action='c';
	      break;
	      
	    case 'e':
	    case 'r': {
	      action='e';
        
        if (strlen(argv[ar]) <= 2 ) {
          fprintf(stderr, "Do  not use space between the -e and the numbers\n");
          ERROR("you must specify a frame, or a range as low-high / low- / -high");
        }
        
        int first = 0, last = MAX_FRAME-1;
        char middle = ' ';
        
        if ( argv[ar][2] == '-' ) {
          if (1 != sscanf(argv[ar]+3,"%u", &last))
            ERROR("you must specify a range as low-high / low- / -high");
        } else {
          switch (sscanf(argv[ar]+2,"%d%c%d", &first, &middle, &last)) {

            case 3: {
              if ( first < 0 )           ERROR("frame number must be positive");
              if ( last  < first )       ERROR("range must be positive, and specified in ascending order");
              if ( last  > MAX_FRAME-1 ) ERROR("specified range is above limit specified at compile time");
              //fprintf(stderr, "3 : %d %c %d \n", first, middle, last );
            } break;

            case 2: {
              if ( first < 0 )           ERROR("frame number must be positive");
              if ( last  < first )       ERROR("range must be positive, and specified in ascending order");
              if ( last  > MAX_FRAME-1 ) ERROR("specified range is above limit specified at compile time");
              //fprintf(stderr, "2 : %d %c %d \n", first, middle, last );
           } break;

            case 1: {
              if ( first < 0 )           ERROR("frame number must be positive");
              if ( first > MAX_FRAME-1 ) ERROR("specified range is above limit specified at compile time");
              last = first;
              //fprintf(stderr, "1 : %u %c %u \n", first, middle, last );
            } break;

            case 0: {
              ERROR("you must specify a frame, or a range as low-high / low- / -high");
            } break;
          }
        }
        
        //store the selection in the list[]
        for (int ii = first; ii < last; ii++) 
          list[ii] = 'p';  //print
        list[last] = 'e';  //print & exit

        //invert the selection if the action is '-r':
        if ( argv[ar][1] == 'r' )  {
          for (int ii=0; ii < MAX_FRAME; ++ii) 
            list[ii] = (list[ii] == 'h') ? 'p' : 'h';
        }
                  
      } break;
	      
	    case 'p': {
	      action='e';
	      int period=2, offset = 0;
	      sscanf(argv[ar]+2, "%d,%d", &period, &offset);
	      if ( period < 2 ) ERROR("period < 2");
	      for (int ii=offset; ii < MAX_FRAME; ii += period)
          list[ii] = 'p';
        } break;
	      
	    case 'l':
	      action='l';
	      break;
	      
	    default:
	      showHelp();
	      break;
	    }
    }
    else
      filename = argv[ar];
  }

  //----------------------------------------------
  infile=fopen(filename, "r");
  if ((infile==0) || ferror(infile)) {
    printf("Error opening file [%s]\n", filename);
    return 0;
  }

  /*
  // output of selection list
  for ( int ii=0; ii < MAX_FRAME; ++ii) {
    if( list[ii] != 'h' ) printf("%i ", ii);
    if( list[ii] == 'e' ) break;
  }
  printf("\n");
  exit(0);
  */
  
  //----------------------------------------------
  switch (action) {
    case 'c':
      countFrame(infile);
      break;
      
    case 'e':
      extractFrame(infile, list);
      break;
      
    case 'l':
      extractLastFrame(infile);
      break;

    default:
      printf("error\n");
      break;
    }
  
  fclose(infile);
  return 0;
}
