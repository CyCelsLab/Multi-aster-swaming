//RCS: $Id: play.cc,v 2.34 2005/04/22 19:28:53 nedelec Exp $
//-----------------------------------------------------------------------
//                               play.cc
//              main file for the graphical part of Cytosim
//-----------------------------------------------------------------------
#include <cctype>

#include "play.h"
#include "play_param.h"
#include "exceptions.h"
#include "offscreen.h"
#include "glextensions.h"
using namespace glExtensions;

#include "play_draw1.cc"
#include "play_draw2.cc"
#include "play_draw3.cc"
#include "play_disp.cc"
#include "play_glut.cc"
#include "play_keys.cc"
#include "play_mouse.cc"
#include "play_menus.cc"

//TODO: use OpenGL accumulation buffer for 'displacement' effects

//-----------------------------------------------------------------------
void showHelp()
{
  printf("------------------------- CYTOSIM-PLAY -------------------------\n");
  printf(" Graphical interactive CYTOSIM, EMBL Heidelberg, nedelec@embl.de\n");
  printf(" Command line options:\n");
  printf("   -l             enter live simulation mode directly\n");
  printf("   -m             offscreen rendering of the movie (all frames)\n");
  printf("   -i             offscreen rendering of one image (combine with frame=[int])\n");
  printf("   -s[int]        set window size (square)\n");
  printf("   -s[int]x[int]  set window width and height\n");
  printf("   -p             display a list of parameters for play\n");
  printf("   -P             display a list of parameters for sim\n");
  printf("   -k             display the list of keyboard commands\n");
  printf("   -v[int]        set the verbose level\n");
  printf("   -h             display this help\n");
  printf(" There should be no space between an option and its argument\n");
  printf(" (compiled at %s on %s, DIM = %i)\n", __TIME__, __DATE__, DIM);
  printf("----------------------------------------------------------------\n");
}


//-----------------------------------------------------------------------
//the simulations goes in live mode only if set from the command line:
bool enterLiveMode = false;

//-----------------------------------------------------------------------
void parseCommandLine(int argc, char * argv[], const bool verbose = true)
{
  for( int ii = 1; ii < argc; ++ii ) {
    
    if ( isalpha( *argv[ii] )) {
      try {
        //first see if argv[ii] specifies a Play parameter:
        if ( PP.parseLineNoCatch( argv[ii], verbose ) )
          //try to interpret argv[ii] as a sim parameter, in a second option:
          if ( MP.parseLine( argv[ii], false ) && verbose )
            MSG("cytosim-play: ignored [%s] on command line\n", argv[ii]);        
      } catch( IOException e ) {
          MSG("cytosim-play: ignored [%s] on command line: %s\n", argv[ii], e.getMessage());
      }
            
    } else {
      
      //check for play major options
      if ( argv[ii][0] == '-' )
        switch( argv[ii][1] ) {
          case 'h': showHelp();                                 exit( EXIT_SUCCESS );
          case 'P': MP.printDescriptions(argv[ii]+2);           exit( EXIT_SUCCESS );
          case 'k': printf("%s", Player::descriptionOfKeys());  exit( EXIT_SUCCESS );
          case 'p': PP.printDescriptions(argv[ii]+2);           exit( EXIT_SUCCESS );
          case 'l': enterLiveMode = true;                       continue;    //live mode
          case 'i': PP.offscreen = 1; PP.showframe = 0;         continue;    //offscreen image
          case 'm': PP.offscreen = 2; PP.showframe = 0;         continue;    //offscreen movie
          case 'v': MSG.setVerboseLevel( argv[ii]+2 );          continue;    //verbose level
          case 's':                                                  //window size
            if (isdigit( argv[ii][2] ))
              if ( 1 == sscanf(argv[ii]+2, "%ux%u", &PP.width, &PP.height ))
                PP.height = PP.width; //if only one value set, make window square
            continue;
        }
      
      //we could do nothing with that option:
      if ( verbose )
        MSG("cytosim-play: WARNING: ignored [%s] on command line\n", argv[ii] );
    }
  }
}


//-----------------------------------------------------------------------
int fixFilePath()
{
  //this function finds the parameter file in the same directory as PP.file
  //and sets MP.datafile accordingly
  
  int code = NO_ERROR;

  //if the filepath is not set:
  if ( * PP.filepath == '\0' ) {
    //we get the path either from MP.datafile  or MP.file:
    if ( enterLiveMode )
      code = IO.resolvePath(PP.filepath, sizeof(PP.filepath), MP.datafile);
    else 
      code = IO.resolvePath(PP.filepath, sizeof(PP.filepath), PP.file);
  }
  if ( code ) return code;
  
  //if not live, we fix the path for the result.out file:
  if ( ! enterLiveMode )
    code = IO.fixPath(PP.file, sizeof(PP.file), PP.filepath);
  
  //we fix the path of MP.datafile:
  code = IO.fixPath(MP.datafile, sizeof(MP.datafile), PP.filepath);

  //we try to open data.in, and if that fails, we try data.out:
  FILE * test = fopen( MP.datafile, "r" );
  if (( test == 0 ) || ferror( test )) {
    MSG("play: MP.datafile=[%s] not found, opening [data.out] instead\n", MP.datafile);
    //we set MP.datafile to try a file "data.out":
    snprintf(MP.datafile, sizeof(MP.datafile), DATA_OUT);
    code = IO.fixPath(MP.datafile, sizeof(MP.datafile), PP.filepath);
  }
  if ( test ) fclose( test );
  
  MSG(6, "fixFilePath(): filepath=%s\n", PP.filepath);
  MSG(6, "fixFilePath(): datafile=%s\n", MP.datafile);
  MSG(6, "fixFilePath(): file=%s\n", PP.file);
  return code;
}
  
  
//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
//-----------------------------------------------------------------------

int main(int argc, char *argv[])
{
  //this suppresses printout of Warnings:
  MSG.shutUp();
    
  //parse the command line for options to play and parameters for sim/play
  //this first parsing is useful to set MP.datafile, the name of the parameter file
  parseCommandLine(argc, argv);

  //try to find the full file-name of the parameter file
  fixFilePath();  

  //we parse the same parameter file for additional play parameters setting:
  //the boolean arguments meaning is (no reset, no verbose)
  PP.parseFile(MP.datafile, false, false);
  
  //---------------start simulation if PP.live = 1, oterwise we open the file
  if ( enterLiveMode ) {

    //set the parameter PP.live, since that is the one driving play:
    if ( PP.live == 0 ) PP.live = 1;
    
    //try to open a valid data file for parameters for SIM:
    if ( NO_ERROR != MP.parseFile(MP.datafile, false) ) {
      MSG("play: Failed to read file [%s]\n", MP.datafile);
      return EXIT_FAILURE; 
    }
    
    //try to open a valid data file for parameters for SIM:
    if ( 0 == MP.nbValuesRead() ) {
      MSG("play: No parameter could be parsed from file [%s]\n", MP.datafile);
      return EXIT_FAILURE; 
    }
    
    //parse again the command line, to overide any value set in data.in
    parseCommandLine(argc, argv, false);

    //we need to correct the path, if datafile was specified on the command line
    fixFilePath();  

    //prepare for the simulation:
    if ( NO_ERROR != Player::prepareToSimulate() ) {
      MSG("play: prepareToSimulate() failed\n");
      return EXIT_FAILURE;
    }
    Player::resetSimulation();
  
  } else {

    //we reset PP.live, which might have been set by parsing datafile
    //live mode must be specified on the command line!
    PP.live = 0;
    
    //if not live, we need the coordinate file:
    if ( NO_ERROR != Player::reader.openFile( PP.file ) ) {
      printf("Cannot open file [%s]\n", PP.file);
      printf("TIP: try 'play -l' to start a live simulation\n");
      printf("      or 'play file=myfile', if the coordinate file is not [%s]\n", RESULT_OUT);
      return EXIT_FAILURE;
    }
    
    //parse sim parameter file, to prepare for a keyboard-triggered sim:
    if ( NO_ERROR != MP.parseFile(MP.datafile) ) {
      MSG.warning("play", "Cannot read sim-parameters from file [%s]", MP.datafile);
    }
    
    //parse again the command line, to overide any value set in data.in
    parseCommandLine(argc, argv, false);
    
    //set the frame buffer, to force loading by display()
    sim.setFrameInBuffer(-1);
  }

 //------------------ off-screen non interactive rendering ---------------
  
  if ( PP.offscreen ) {
    if ( OffScreen::open(PP.width, PP.height) ) 
      return EXIT_FAILURE;
    
    //initialize the window:
    Player::reshapeWindow(PP.width, PP.height);
    Player::initGL();
    
    switch ( PP.offscreen ) {
    case 1:
      Player::display();
      Player::saveImageAsPPM();
      break;
    case 2:
      Player::saveMovieAsPPM();
      break;
    }
    OffScreen::close();
    return EXIT_SUCCESS;
  }

  //--------------------- on-screen interactive rendering ------------------
  //------------- calls GLUT to open a screen & start displaying -----------
  
  glutInit(&argc, argv);

  //set the display mode to initialize GLUT rendering:
  unsigned int mode = GLUT_RGBA;
  if ( PP.buffered ) mode = mode | GLUT_DOUBLE;
  if ( DIM == 3 )    mode = mode | GLUT_DEPTH;
  glutInitDisplayMode( mode );

  //init the window:
  glutInitWindowSize(PP.width, PP.height); 
  glutInitWindowPosition(50, 50);
  glutCreateWindow(argv[0]);
  if ( PP.fullscreen ) glutFullScreen();

  //init the other stuff:
  Player::initGL();
  Player::buildMenus();

  //register all the GLUT callback functions:
  glutDisplayFunc(Player::display);
  glutReshapeFunc(Player::reshapeWindow); 
  glutMouseFunc(Player::processMouse);
  glutMotionFunc(Player::processMotion);
  glutSpecialFunc(Player::processInputKey);
  glutKeyboardFunc(Player::processNormalKey);
  glutTimerFunc(PP.delay, Player::timerFunction, 1);
  
  //start the GLUT window handling:
  glutMainLoop();

  return EXIT_SUCCESS;
}
