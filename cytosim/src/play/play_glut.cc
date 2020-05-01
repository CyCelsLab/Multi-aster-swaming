//RCS: $Id: play_glut.cc,v 2.31 2005/04/24 22:48:26 nedelec Exp $
//------------------------------play_glut.cc-----------------------------
#include <cstdarg>

//-----------------------------------------------------------------------
void Player::flashText( const char fmt[], ...)
{
  va_list args;
  
  if ( PP.live ) {
    //show the text for one second on the screen:
    if ( PP.delay > 10 )
      PP.showflashtext = 1 + 1000 / PP.delay;
    else
      PP.showflashtext = 30;
  } else {
    PP.showflashtext = 1;
  }
  
  va_start(args, fmt);
  vsnprintf(PP.flashtext, sizeof(PP.flashtext), fmt, args);
  va_end(args);
}


//-----------------------------------------------------------------------
void Player::setModelView(const bool update)
{
  static GLdouble mat[16];
  glMatrixMode(GL_MODELVIEW);
  PP.quat.normalizeit();
  PP.quat.setThisMatrix16( mat, cameraTranslation );
  glLoadMatrixd( mat );
  glScaled( PP.zoom, PP.zoom, PP.zoom );
  glTranslate( -PP.focus[0], -PP.focus[1], -PP.focus[2]);
  if ( update ) 
    glutPostRedisplay();
}

//-----------------------------------------------------------------------
void Player::resetView()
{
  PP.quat.set(1);
  PP.focus.set(0,0,0);
  PP.zoom = 1;
  setModelView();
}

//-----------------------------------------------------------------------
void Player::swapFullScreen()
{
  static int width_saved = 512, height_saved = 512;
  static int position_x = 10, position_y = 10;
  
  //swap between full-screen and small-window
  if (PP.fullscreen) {
    //invoke small-window with saved sizes
    glutReshapeWindow( width_saved, height_saved );
    glutPositionWindow( position_x, position_y );
    PP.fullscreen = false;
  } else {
    //save the current window size for the next call of swapFullScreen:
    width_saved  = PP.width;
    height_saved = PP.height;
    position_x   = glutGet(GLUT_WINDOW_X);
    position_y   = glutGet(GLUT_WINDOW_Y);
    //invoke full screen function from GLUT
    glutFullScreen();
    PP.fullscreen = true;
  }
}

//-----------------------------------------------------------------------
//--------------------------  General commands --------------------------
//-----------------------------------------------------------------------

void Player::prevFrame()
{
  if ( PP.frame > 0 )
    { --PP.frame; glutPostRedisplay(); }
  else {
    if ( PP.loop )
      { PP.frame = reader.lastGoodFrame(); glutPostRedisplay(); }
    else
      PP.dir = PLAY_STOP;
  }
}

//-----------------------------------------------------------------------
void Player::nextFrame()
{
  //if PP.frame is set to 0, and we have a different
  //frame in the buffer, the best guess is to rewind the file:
  if ( reader.eof() )
    reader.rewindFile();

  //we do not just do ++PP.frame, but read the next frame from the
  //current file position. This allows allow play to jump missing frames
  switch ( reader.readNextFrame() ) { 
      case NO_ERROR:
          //we set PP.frame to avoid the automatic read in display():
          PP.frame = sim.frameInBuffer(); 
          glutPostRedisplay();
          return;
          
      case 1: //EOF:
          if ( PP.loop )
            PP.frame = 0;
          else
            PP.dir = PLAY_STOP;
          glutPostRedisplay();
          return;
          
      default: //ERROR
          ++ PP.frame;
          glutPostRedisplay();
          return;
  }
}

//-----------------------------------------------------------------------
//save the image in PPM binary format. PPM can be converted later using
//for example the library http://www.acme.com/software/pbmplus/

int Player::saveImageAsPPM()
{
  static int lastFrame = -1;
  int indxFrame = sim.frameInBuffer();
  
  //create a name for the image file
  char filename[STRING_SIZE];
  //printf("filepath = %s frame_in_buffer %i\n", PP.filepath, indxFrame);
  
  //we make sure we do not overwrite a frame already produced previously: 
  if ( indxFrame <= lastFrame )
    indxFrame = ++lastFrame;
  else
    lastFrame = indxFrame;

  if ( * PP.filepath == '\0' )
    snprintf(filename, sizeof(filename), "image%04i.ppm", indxFrame);
  else
    snprintf(filename, sizeof(filename), "%s/image%04i.ppm", PP.filepath, indxFrame);

  //open the file for binary writing:
  FILE * file = fopen(filename, "wb");
  if (( file == 0 ) || ( ferror( file ))) {
    printf("Cytosim::play could not open file [%s]\n", filename);
    return 1;
  }
  
  int code = glExtensions::saveImageAsPPM(file, 0, 0, PP.width, PP.height);
  if ( NO_ERROR == code ) {
    //printf to the terminal window for user info:
    MSG("cytosim: saved image at time = %.3f sec. in %s\n", sim.simTime(), filename );
    return 0;
  } else return code;
}


//-----------------------------------------------------------------------

void Player::saveMovieAsPPM()
{
  if ( reader.getInputFile() == 0 ) return;
  flashText("");
  
  //make sure the current frame is loaded:
  reader.readFrame(PP.frame);

  try {
    do {
      //we can escape the automatic readFrame() in display(),
      //by setting the requested frame to what is in the buffer:
      PP.frame = sim.frameInBuffer();
      //we directly call display();
      display();
      saveImageAsPPM();
      //read next frame until EOF comes
    } while ( 1 != sim.readState() );
  } catch( IOException e ) {
    MSG("Exception in saveMovieAsPPM");
    e.print();
  }
}

//-----------------------------------------------------------------------
void Player::writeToFile()
{
  try {
    sim.writeState("start.out");
    MSG("State written to file <start.out>\n");
  } catch( IOException e ) {
    MSG("Exception in writeToFile()");
    e.print();
  }
}

//-----------------------------------------------------------------------
void Player::appendToFile()
{
  try {
    IO.openOutputFile("start.out", 1);
    sim.writeState();
    IO.closeOutputFile();
    MSG("State written to file <start.out>\n");
  } catch( IOException e ) {
    MSG("Exception in appendToFile()");
    e.print();
  }
}

//-----------------------------------------------------------------------
void Player::readFromFile()
{
  try {
    sim.readState(START_IN);
    MSG("Simulation state from file %s\n", START_IN);
  } catch( IOException e ) {
    MSG("Exception in readFromFile()");
    e.print();
  }
}

//-----------------------------------------------------------------------
int Player::prepareToSimulate(const bool forceIt)
{
  int code = NO_ERROR;
  //check that the parameter file was read:
  if ( 0 == MP.nbValuesRead() ) {
    flashText("Cannot start simulation: unset parameter values (file %s)\n", MP.datafile);
    return 1;
  }
  //initialize the simulation:
  if ( ! sim.isReady() || forceIt ) {
    code = sim.getReady();
    //we avoid reloading the sim-state from the file:
    sim.setFrameInBuffer( PP.frame );
  }
  return code;
}

//-----------------------------------------------------------------------
void Player::resetSimulation()
{
  if ( sim.isReady() ) {
    sim.eraseState();
    sim.populateState();
    sim.resetCounters();
    playGraftedList.clear();
    playGrafted=0;
  }
}

//-----------------------------------------------------------------------
//---------------------------------- initGL -----------------------------
//-----------------------------------------------------------------------

void Player::initGL()
{
  switch( PP.style ) {
  case 3: initGL3(); break;
  case 2: initGL2(); break;
  default:    
  case 1: initGL1(); break;
  }
  setModelView();
  glPrintErrors("initGL");
}


//-----------------------------------------------------------------------
void Player::reshapeWindow(const int w, const int h)
{
  PP.width  = w;
  PP.height = h;

  glViewport(0, 0, w, h);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  real ratio = w * visibleRegion[1] / real( visibleRegion[0] * h );

  real znear = -cameraTranslation[2] - visibleRegion[2];
  real zfar  = -cameraTranslation[2] + visibleRegion[2];

  if ( ratio > 1 )
    glOrtho(-visibleRegion[0], visibleRegion[0],
            -visibleRegion[1]/ratio, visibleRegion[1]/ratio, znear, zfar);
  else
    glOrtho(-visibleRegion[0]*ratio, visibleRegion[0]*ratio,
            -visibleRegion[1], visibleRegion[1], znear, zfar);

  if ( PP.fog ) {
    glFogf(GL_FOG_DENSITY, 0.20/(zfar - znear) );
    glFogf(GL_FOG_START, znear);
    glFogf(GL_FOG_END,   zfar);
  }

  //get the transformation matrices, to be used for mouse control
  glGetIntegerv(GL_VIEWPORT, viewport);
  glGetDoublev(GL_PROJECTION_MATRIX, projectionMat);
  glMatrixMode(GL_MODELVIEW);
}



//-----------------------------------------------------------------------
void Player::timerFunction(const int value)
{
 //MSG("timer live = %i  sim.isReady() = %i\n", PP.live, sim.isReady());
    if ( PP.live ) {
        if ( sim.isReady() ) {
            //perform 2^live simulation steps before displaying:
            for( int s = 1 << ( PP.live - 1 ); s > 0; --s ) 
                sim.stepAll();
            glutPostRedisplay();
        }
    } else {
        switch( PP.dir ) {
            case PLAY_REVERSE:
                prevFrame();
                break;
                
            case PLAY_STOP: 
                break;
                
            case PLAY_FORWARD:
                nextFrame();
                break;
              
            case PLAY_FORWARD_WRITE:
                saveImageAsPPM();
                nextFrame();
                break;
        }
    }

  //register another timer call back in PP.delay milli-sec:
  if ( ! finish_now )
    glutTimerFunc( PP.delay, Player::timerFunction, 1);
}
