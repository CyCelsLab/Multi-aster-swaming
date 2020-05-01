//RCS: $Id: readwrite.cc,v 2.12 2005/04/22 19:31:02 nedelec Exp $
//------------------just read and write coordinate file----------------------
//    mostly for testing, but can also make some modifications on state
//    F. Nedelec, EMBL 2002

#include "sim.h"
#include "exceptions.h"
#include "iowrapper.h"

void showHelp()
{
    printf("readwrite [options] file_in file_out\n");
    printf("line command options:\n");
    printf("        -h     display this message\n");
    printf("        -t     generate a text file <file_out>\n");
    printf("        -b     generate a binary file <file_out>\n");
    printf("        -v#    set the verbose level\n");
    printf("        -n     renumber the frame consecutively\n");
}

int main(int argc, char * argv[])
{
    int renumber = 0;
    int cnt = 0;
    const char * namein  = RESULT_OUT;
    char * nameout = "result2.out";
    
    for( int ii = 1 ; ii < argc ; ++ii ) {
      if (( *argv[ii] == '-' ) && ( strlen( argv[ii] ) > 1 )) {
        switch( argv[ii][1] ) {
          case 'n': renumber = 1;                     break;
          case 't': IO.produceBinaryOutput(false);    break;
          case 'b': IO.produceBinaryOutput(true);     break;
          case 'h': showHelp();                       return 0;
          case 'v': MSG.setVerboseLevel(argv[ii]+2);  break;
          default: printf("ignored [%s] on command line\n", argv[ii]);
        }
      } else {
        switch( cnt++ ) {
          case 0: namein = argv[ii]; break;
          case 1: nameout = argv[ii]; break;
        }
      }
    }
    
    try {
      IO.openInputFile(namein);
    } catch( IOException e ) {
      printf("could not open input file [%s]: %s\n", namein, e.getMessage());
      return EXIT_FAILURE;
    }
      
    try {
      IO.openOutputFile(nameout);
    } catch( IOException e ) {
      printf("could not open output file [%s]: %s\n", nameout, e.getMessage());
      return EXIT_FAILURE;
    }
      
    printf("reading %s -> writing %s\n", namein, nameout);
    cnt = 0;
    
    while ( ! IO.eof() )   {
        
      try {
        if ( NO_ERROR == sim.readState() ) {
          MSG("Frame %3i with %i mts\n", sim.frameInBuffer(), sim.nbMicrotub());
          if ( renumber ) sim.setFrameInBuffer( cnt );
          sim.writeState();
          ++cnt;
        }
      } catch( IOException e ) {
        MSG("ERROR reading frame %i: %s\n", sim.frameInBuffer(), e.getMessage());
      }
      
    }
    
    MSG("%i frames written to [%s]\n", cnt, nameout);
    return 0;
}
