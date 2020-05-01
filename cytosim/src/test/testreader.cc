//RCS: $Id: testreader.cc,v 2.11 2005/04/18 11:41:48 foethke Exp $
//----------------------------------------------------------------------------
//                             testreader.cc
//
//     this is mostly a test for the class defined in reader.h
//     but it can be used to navigate from frame to frame in a result.out file
//----------------------------------------------------------------------------
#include <cstring>
#include <cctype>
#include <cstdlib>

#include "iomessages.h"
#include "iowrapper.h"
#include "reader.h"
#include "sim.h"

void showHelp()
{
    printf("readwrite [options] file_in file_out\n");
    printf("line command options:\n");
    printf("        -h     display this message\n");
    printf("        -t     generate a text file <file_out>\n");
    printf("        -t     generate a binary file <file_out>\n");
    printf("        -v#    set the verbose level\n");
    printf("\n");
    printf("enter frame number to jump to a specified frame in file_in\n");
    printf("enter w to write it to file_out\n");
    printf("enter w to quit\n");
}

void showCommands()
{
  printf("would you mind some instructions ?\n");
  printf("  'q'    quit\n");
  printf("  'n'    read next frame\n");
  printf("  'w'    write frame\n");
  printf("  'c'    clear buffer of frame starting positions\n");
  printf("  'e'    erase state\n");
  printf(" number  try to read the specified frame\n");
  printf(" hint:  try starting testreader with -v11 for verbose indications\n");  
}



int main(int argc, char * argv[])
{
    int cnt = 0;
    const char * namein  = RESULT_OUT;
    char * nameout = "result2.out";
    
    //parse the command line:
    for( int i = 1 ; i < argc ; ++i ) {
        if ( ( *argv[i] == '-' ) && ( strlen( argv[i] ) > 1 ) )
            switch( argv[i][1] ) {
                case 't': IO.produceBinaryOutput(false); break;
                case 'b': IO.produceBinaryOutput(true);  break;
                case 'h': showHelp(); return 0;
                case 'v': MSG.setVerboseLevel( argv[i]+2 ); break;
            }
            else switch( cnt++ )  {
                case 0: namein = argv[i]; break;
                case 1: nameout = argv[i]; break;
            }
    }
        
    char userinput[STRING_SIZE];
    SimReader reader( namein );
    int frame;
    
    printf("TestReader: read/write frame for cytosim. Enter (h) for help\n");
    while( true ) {
        
      if ( sim.frameInBuffer() < 0 )
        printf("No frame in buffer\n");
      else
        printf("Frame %i in buffer, time = %7.3f sec. %i Microtubs\n", 
               sim.frameInBuffer(), sim.simTime(), sim.nbMicrotub());
      printf("%s", sim.frameInfo());
      //sim.printShortDescription();
        
      printf("?");
      fgets( userinput, sizeof(userinput), stdin );
        
      if ( isdigit( *userinput )) {
        if ( 1 == sscanf(userinput, "%i", &frame ) ) {
          if ( NO_ERROR != reader.readFrame( frame ) )
            printf("frame not found : ");
        }
      } else {
        switch( *userinput ) {
          case '\n':
          case 'n':
            if ( NO_ERROR == reader.readNextFrame() )
              printf("next frame : ");
            break;
            
          case 'w':
            if ( IO.getOutputFile() == 0 ) 
              IO.openOutputFile(nameout);
            sim.writeState();
            break;
            
          case 'e':
            sim.eraseState();
            break;
            
          case 'c':
            reader.clearBuffer();
            break;
            
          case 'q': case 'Q': case 27:
            return 0;
            
          default:
            showCommands();
            break;
        }
      }
    }
    return 0;
}
