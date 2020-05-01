//RCS: $Id: analyse_nucpos.cc,v 2.8 2005/04/22 19:30:53 nedelec Exp $

#include "sim.h"
#include "reader.h"
#include "iowrapper.h"
#include "iomessages.h"

const unsigned int STRLEN = 512;
char commands[STRLEN];

int first_frame = 0;
int last_frame  = SimReader::MAX_FRAME;

//===================================================================


void NucleusPosition()
{
  Nucleus* nu = sim.findNucleus(1);
  //if( DIM == 2 )
  //printf("%f %f\n", nu->coord(0), nu->coord(1));
  //else if( DIM == 3 )
  //printf("%f %f %f\n", nu->coord(0), nu->coord(1), nu->coord(2));
  printf("%f %f\n", sim.simTime(), nu->coord(0));
}

void MTTouchTips()
{
  Vecteur plusend;
  
  //printf("reached MTTouchTips!\n");
  
  for ( Microtub * mt=sim.firstMicrotub(); mt ; mt=mt->next() ) {
    plusend = mt->whereEnd(MT_PLUS_END);
    //if (plusend.XX > fabs(MP.boxsize[0] - MP.boxsize[1])) {
    if (plusend.XX > fabs(MP.boxsize[0])) {
      if(mt->space->isOutside(plusend)) {
        printf("Microtubule is hitting cortex at %fsec!\n", sim.simTime());
      }
    }
    }
  }

//=====================================================================

void showHelp()
{
  printf("Cytosim / analyse:\n");
  printf("Reads [%s] and prints some statistics on the objects\n", RESULT_OUT);
  printf("USAGE: analyse [ -f(index) ] [-f(index1-index2)] [ command_name1 command_name2 ... ]\n");
  printf("\n");
  printf("Commands supported:\n");
  printf(" nucleus          report position of nucleus\n");
  printf(" mttips           report microtubule touching cortex\n");
  printf("\n");
  printf("analyse_nucpos  DIM=%i, compiled %s, %s\n", DIM, __TIME__, __DATE__ );
}

//=====================================================================

void parseCommandLine(int argc, char * argv[])
{
  for(int ii = 1; ii < argc; ++ii ) {
    
    if ( '-' == argv[ii][0] ) {
      switch( argv[ii][1] ) {
        case 'h': showHelp(); exit(EXIT_SUCCESS);
        case 'v': MSG.setVerboseLevel( argv[ii]+2 ); continue;
        case 'f':   //frame index
          if ( argv[ii][2] )
            if ( 1 == sscanf(argv[ii]+2, "%u-%u", &first_frame, &last_frame ) )
              last_frame = first_frame; //if only one value set, make window square
          continue;
          
        default:
          printf("ignored option [%s] on command line\n", argv[ii] );
          continue;
      }
    }
    
    if ( strlen( commands ) + strlen( argv[ii] ) + 1 < STRLEN ) {
      strcat( commands, " " );
      strcat( commands, argv[ii] );
    }
  }
  strcat( commands, " " );
}

//=====================================================================
int main(int argc, char * argv[])
{
  if ( argc < 2 ) { 
    showHelp(); 
    return EXIT_SUCCESS;
  } 
  
  parseCommandLine(argc, argv);

  MSG.shutUp();
    
  if ( NO_ERROR != MP.parseFile(DATA_OUT) )
    if ( NO_ERROR != MP.parseFile(MP.datafile) )
      MSG.error("analyse_nucpos.cc","Cannot read %s or data.out", MP.datafile);

  MP.massageParameters();

  SimReader reader;
  if ( NO_ERROR != reader.openFile( RESULT_OUT ) ) {
    fprintf(stderr, "analyse_nucpos.cc: cannot open %s\n", RESULT_OUT);
    return EXIT_FAILURE; 
  }

  for( int frame = first_frame; frame <= last_frame ; ++frame )
    if ( NO_ERROR == reader.readFrame( frame ) ) {
      
      if ( strstr(commands, " nucleus " ))       NucleusPosition();
      if ( strstr(commands, " mttips " ))        MTTouchTips();
      //add a line here to call your analyse procedures
      
    }
  return EXIT_SUCCESS;
}
