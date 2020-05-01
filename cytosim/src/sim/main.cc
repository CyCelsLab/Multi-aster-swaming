//RCS: $Id: main.cc,v 2.27 2005/04/22 19:30:43 nedelec Exp $
//---------------------------main.cc-------------------------------//
// main for Cytosim - non interactive
// Francois Nedelec 1997-2004  Email: nedelec@embl.de
// EMBL, Meyerhofstrasse. 1, 69117 Heidelberg, Germany.
// http://www.cytosim.org



// Neha Khetan, 3 September 2016 - adding function to output Gh properties
// Neha Khetan, 15 September 2016 -adding to output Gh properties in specific interval

#include "sim.h"
#include "iowrapper.h"
#include <ctime>
#include <cctype>


void reportCompileOptions()
{
  MSG(3, "Compiled at %s, %s, DIM=%i\n", __TIME__, __DATE__, DIM);
  MSG(3, "Hydrodynamic correction (MT) = %.3f\n", HYDRO);

#ifdef NEAR_ATTACH_FASTER
  MSG(3, "NEAR_ATTACH_FASTER : YES\n");
#else
  MSG(3, "NEAR_ATTACH_FASTER : NO\n");
#endif

#ifdef SINGLE
  MSG(3, "Precision: SINGLE (real=float, %i bytes) Epsilon = %e\n", sizeof(real), EPSILON);
#else
  MSG(3, "Precision: DOUBLE (real=double, %i bytes) Epsilon = %e\n", sizeof(real), EPSILON);
#endif

#ifdef UNIFORM_FREE_COMPLEX
  MSG(3, "UNIFORM_FREE_COMPLEX approximation\n");
#endif

#ifdef ASSERT         //this keyword is defined in "assert_macro.h"
  MSG(3, "This is a SAFE version : assertions are enabled (in assert_macro.h)\n");
#endif
}


//===========================================================================


void reportCPUtime()
{
  static clock_t old_clock  = clock();
  static unsigned long cpu = 0;
  time_t now = time(0);

  MSG("F %3i  %7.2fs %3i  ", sim.frameInBuffer(), sim.simTime(),
      (100*sim.iterationCnt())/MP.nbiter );

  char date[32];
  strncpy( date, asctime( localtime( &now )), 19);
  date[19]='\0';
  MSG("%s ", date );

  unsigned long dcpu = ((unsigned long)clock() - (unsigned long)old_clock );
  dcpu  /= CLOCKS_PER_SEC / 10;
  old_clock = clock();
  cpu += dcpu;
  MSG("CPU %5d : %7d  %5.2fH\n", dcpu/10, cpu/10, (cpu/360)/100.0);

}


//===========================================================================


void showHelp()
{
  MSG("------------------------------------------------------------\n");
  MSG("CYTOSIM version 2.4-1 Neha Khetan, 2014-2019 \n");
  MSG("CYTOSIM version 2.4chai  June   2008     http://www.embl.de/~athale\n");
  MSG("Chaitanya Athale\n");
  MSG("based on code by Francois J. Nedelec (nedelec@embl.de) with\n");
  MSG("Dietrich Foethke, Cleo Kozlowski, Maria Mora, Thomas Clausen\n");
  MSG("------------------------------------------------------------\n");
  MSG("command line options:\n");
  MSG("     [int]    set the number of iterations to the specified integer\n");
  MSG("     -t       produce text coordinate file in [%s]\n", RESULT_OUT);
  MSG("     -b       produce binary coordinate file in [%s]\n", RESULT_OUT);
  MSG("     -o       send output to terminal rather than in file <sim.out>\n");
  MSG("     -vINT    set the verbose level\n");
  MSG("     -p       print common parameters descriptions and default values\n");
  MSG("     -P       print all parameter descriptions and default values\n");
  MSG("     -A       make all parameters visible in <data.out>\n");
  MSG("     -h       print this message\n");
  MSG("\n(compiled at %s, %s, DIM=%i)\n", __TIME__, __DATE__, DIM);
}


//===========================================================================
//================================  MAIN  ===================================
//===========================================================================

int main(int argc, char* argv[])
{
  ParameterModifier  print_all_parameters = PARAM_ALL;     //to file 'data.out'
  int                output_to_file       = 1;             //to file 'sim.out'

  //-----------------parse the command line for global options--------------

  for( int ii = 1 ; ii < argc ; ++ii ) {
    if ( argv[ii][0] == '-' ) {
      switch( argv[ii][1] ) {
        case 'A': print_all_parameters = PARAM_ALL_VALUES;  continue;
        case 't': IO.produceBinaryOutput(false);            continue;
        case 'b': IO.produceBinaryOutput(true);             continue;
        case 'p': MP.printDescriptions(argv[ii]+2);                     return EXIT_SUCCESS;
        case 'P': MP.printDescriptions(argv[ii]+2, stdout, PARAM_ALL);  return EXIT_SUCCESS;
        case 'h': showHelp();                                           return EXIT_SUCCESS;
        case 'o': output_to_file = 0;                       continue;
        case 'v': MSG.setVerboseLevel( argv[ii]+2 );        continue;
        default:
          MSG("ignored option [%s] on command line\n", argv[ii] );
          continue;
      }
    } else {
      //we check for parameter sets, only for MP.datafile,
      //which defines the file from which parameters will be read
      if ( isalpha(*argv[ii]) ) {
        if ( NO_ERROR != MP.parseLine( argv[ii] ))
          MSG("ignored option [%s] on command line\n", argv[ii] );
      }
    }
  }

  if ( output_to_file )
    MSG.openFile("sim.out");

  MSG("========================COMPILE OPTIONS============================\n");

  reportCompileOptions();

  MSG("=============================PARAMS================================\n");

  if ( NO_ERROR != MP.parseFile( MP.datafile )) {
    MSG("input parameter file [%s] not found or invalid\n", MP.datafile);
    return EXIT_FAILURE;
  }

  if ( MP.nbValuesRead() == 0 ) {
    MSG("parameter file [%s] empty or invalid\n", MP.lastFileParsed());
    return EXIT_FAILURE;
  }

  //-------------------------parse the command line for parameter values:
  //we scan the command line for parameter sets, to override the file
  for( int ii = 1 ; ii < argc ; ++ii ) {
    if ( isalpha(*argv[ii]) ) MP.parseLine( argv[ii] );
    if ( isdigit(*argv[ii]) ) sscanf(argv[ii], "%lu", &MP.nbiter);
  }

  MP.massageParameters();
  MP.printSome();

  MSG("========================== INITIALIZING ===========================\n");

  sim.getReady();
  sim.populateState();
 //sim.printShortDescription();
 
	

  if ( MP.record == 0 ) MP.record = 1;

  int frame     = 0;
  int statFrame = 0;
  int ghstatFrame = 0;

  real steps_per_frame       = real( MP.nbiter ) / MP.record;
  real steps_per_statFrame   = 1;
  real steps_per_GhstatFrame = 1;
 


  if( MP.statistics[1] > 0 )
    //record only MP.statistics[1] statistical frames
    steps_per_statFrame = real( MP.nbiter ) / MP.statistics[1];

  Step next_record     = 0; //when we should record the next frame
  Step next_recordStat = 0; //when we should record the next statistics frame


  // Neha Khetan, 15 September 2016 - adding the interval at which ghstats should be outputted
  if( MP.ghstatistics[1] > 0 )
    //record only MP.sghtatistics[1] statistical frames
     steps_per_GhstatFrame = real( MP.nbiter ) / MP.ghstatistics[1];   
  Step next_GhrecordStat = 0; //when we should record the next statistics frame






  MSG("total real time %.3f s, dt=%.2e s, %i intervals of %.3f s\n",
      MP.nbiter*MP.dt, MP.dt, MP.record, MP.dt * steps_per_frame );

  MSG("=============================RUNNING===============================\n");

  

  

  do {

    //recording of a frame to disc:
    if ( sim.iterationCnt() >= next_record ) {

      if ( sim.iterationCnt() >= (Step)MP.nbiter ) break;

      reportCPUtime();
      sim.writeStateAuto();

      next_record = Step( ++frame * steps_per_frame );
      if ( next_record > (Step)MP.nbiter )
        next_record = MP.nbiter;
    }

    //record statistics on the fly:
    if ( MP.statistics[0] > 0 ) {
      if ( sim.iterationCnt() >= next_recordStat ) {
        sim.recordStatistics();		// Neha Khetan, 3 September 2016: to output statistics of motors
        next_recordStat = Step( ++statFrame * steps_per_statFrame );
        if ( next_recordStat > (Step)MP.nbiter )
          next_recordStat = MP.nbiter;
      }
    }
    // Neha Khetan, 15 September 2016 - adding the interval at which ghstats should be outputted
    //record Ghstatistics on the fly:
    if ( MP.ghstatistics[0] > 0 ) {
      if ( sim.iterationCnt() >= next_GhrecordStat ) {        
	       sim.recordGhStatistics();   // Neha Khetan, 3 September 2016: to output statistics of motors
         next_GhrecordStat = Step( ++ghstatFrame * steps_per_GhstatFrame );
         if ( next_GhrecordStat > (Step)MP.nbiter )
          next_GhrecordStat = MP.nbiter;
      }
    }




    //perform the simulation step:
    sim.stepAll();

  } while ( finish_now == false );



  //clean-up:
  reportCPUtime();
  sim.writeStateAuto();
  if ( MP.statistics[0] > 0 )
    sim.recordStatistics(true);
  MP.printFile("data.out", print_all_parameters);
  MSG("end\n");

  sim.eraseState();
  sim.clearSpace();
  sim.deleteMTGrid();
  return EXIT_SUCCESS;
}
