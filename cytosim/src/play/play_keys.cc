//RCS: $Id: play_keys.cc,v 2.38 2005/04/12 20:42:07 nedelec Exp $

//--------------------------------------------------------------------------
//-------------------------- keyboard commands -----------------------------
//--------------------------------------------------------------------------

const char * Player::descriptionOfKeys()
{
  const unsigned int SZ = 2660;
  static char * hS = 0;

  if ( hS == 0 ) {
    hS = new char[SZ];
    if (hS == 0) {
      fprintf(stderr, "play_keys.cc::buildMenus() memory allocation failed\n");
      exit(1);
    }

    *hS = '\0';

    ///\todo: replace strcat by snprintf() for safety reasons
    strcat(hS,"--------------------------------Display-------------------------------------\n");
    strcat(hS," Use the mouse to rotate and translate the objects space:\n");
    strcat(hS,"      left-button                 + drag: translation\n");
    strcat(hS,"      middle-button ( alt-click)  + drag: zoom in and out\n");
    strcat(hS,"      right-button  (ctrl-click)        : pop-up menu\n\n");
    strcat(hS," + -         Zoom in and out (small increament if SHIFT is down)\n");
    strcat(hS,"             +CTRL : increase / decrease window size\n");
    strcat(hS," arrow keys  Translate in (X,Y) plane, press SHIFT for fine motion\n");
    strcat(hS," v space     Reset view, refresh display\n");
    strcat(hS," h           Hide/Show this message on display, show the parameters\n");
    strcat(hS,"Microtubules:\n");
    strcat(hS," 1           Rotate between different display modes\n");
    strcat(hS," !           Rotate +/- ends display: none, plus-ends, minus-ends\n");
    strcat(hS," 2 3         Decrease/increase line width\n");
    strcat(hS," @ #         Decrease/increase point size\n");
    strcat(hS,"Complexes:\n");
    strcat(hS," 8           Rotate between displaying all, free, bound or bridging complexes\n");
    strcat(hS," 9 0         Decrease/increase hand point size\n");
    strcat(hS," ( )         Decrease/increase complexes link width\n");
    strcat(hS,"Grafteds:\n");
    strcat(hS," 5           Rotate between displaying all, free or bound grafted\n");
    strcat(hS," 6 7         Decrease/increase hand size\n");
    strcat(hS," ( )         Decrease/increase hand link width\n");
    strcat(hS,"Solids:\n");
    strcat(hS," %           Switch between different solid display mode\n");
    strcat(hS," ^ &         Decrease/increase point size\n");
    strcat(hS,"-------------------------------Animation------------------------------------\n");
    strcat(hS," < >         Show previous / next frame ( , . also work)\n");
    strcat(hS," u i o p     reverse, stop, slower, play, faster\n");
    strcat(hS," z l         Rewind movie to first frame, toggle looping\n");
    strcat(hS," q           Quit \n");
    strcat(hS,"-------------------------------Simulation-----------------------------------\n");
    strcat(hS," a A         Start live simulation, increase steps/display ratio\n");
    strcat(hS," s           Perform one simulation step, also stops live mode\n");
    strcat(hS," x           Read file <data.in>, and update parameters\n");
    strcat(hS," z           Reset the simulation, creating a new initial state\n");
    strcat(hS,"-------------------------------Read/Write-----------------------------------\n");
    strcat(hS," y Y         Save image as a PPM file, save all movie at PPM\n");
    //strcat(hS," w W         Append/write current state to file <start.out>\n");
    //strcat(hS," r           Read state from file <start.in>\n");
    strcat(hS," B b         Increase verbose level of cytosim messages, supress it\n");
    strcat(hS,"------------------------------Interaction-----------------------------------\n");
    strcat(hS," in 2D, shift-click controls a strong undetachable grafted-motor\n");
    strcat(hS," g G         Create a new mouse-controlled grafted motor / delete them all\n");

    //We check a posteriori that the allocation size was sufficient
    if ( strlen(hS) >= SZ )
      MSG.error("play_keys.cc:descriptionOfKeys()", "internal error: set SZ greater than %i\n", strlen(hS));
  }
  return hS;
}


//--------------------------------------------------------------------------------------
void Player::processNormalKey(const unsigned char c, const int x, const int y)
{
  if ( PP.keyboard == 0 ) return;

  //in the switch below:
  //  use break;  glutPostRedisplay() will be called, to redraw the frame
  //  use return; if redrawing is not necessary.

  switch (c) {

  case 'q': case 'Q': case 27: // ascii 27 is the escape key
    exit( EXIT_SUCCESS );

  case 'm':  //this is the magic key
    sim.link( new Microtub(1) );
    flashText("one more Microtub!");
    break;

  case 'M':  //this is the cursed key
    delete( sim.firstMicrotub() );
    flashText("lost one Microtub...");
    break;

  case ' ':  //force reload of the frame from the file
    reader.clearBuffer();
    sim.setFrameInBuffer(-1);
    flashText("");
    break;

  case 'h':
    if ( PP.showhelp ) {
      PP.showhelp = false;
      PP.showparameters = true;
    } else {
      if ( PP.showparameters )
        PP.showparameters = false;
      else
        PP.showhelp = true;
    }
    break;

  case 'f':
    swapFullScreen();
    break;

  case 'y':
    saveImageAsPPM();
    break;

  case 'Y':
    PP.dir = PLAY_FORWARD_WRITE; //will save all frames (cf timer function).
    //saveMovieAsPPM();
    return;

  case 'v':           //reset view
    resetView();
    PP.autozoom = 1;
    flashText("Reset view");
    break;

  //------------------------- next / prev frame:


  case '<': case ',':
    prevFrame();
    return;

  case '>': case '.':
    nextFrame();
    return;

    //------------------------- live simulation mode:

  case 'x': {  //read and parse parameters from 'data.in' or equivalent
    if ( NO_ERROR == MP.parseFile(MP.datafile) ) {
      if ( NO_ERROR == prepareToSimulate( true )) //force the reloading
        flashText("Parameters set from %s", MP.datafile);
      else
        flashText("Read file %s, but failed to prepareToSimulate()", MP.datafile);
    } else {
      flashText("Failed to read %s", MP.datafile);
    }
  } break;

  case 'z':  //reset the state, or rewind to frame=0, and stop the movie
    PP.frame = 0;
    if ( PP.live && sim.isReady() ) {
      resetSimulation();
      flashText("Reset state");
    } else {
      PP.dir = PLAY_STOP;
    }
    break;

  case 'a':    //starts the live mode
    if ( NO_ERROR == prepareToSimulate() )
      PP.live = 1;
    else
      flashText("prepareToSimulate() failed");
    break;

  case 'A':    //increase the number of simulation between displays.
    if ( NO_ERROR == prepareToSimulate() )
      PP.live = 1 + PP.live % 6;
    else
      flashText("prepareToSimulate() failed");
    break;

  case 's':   //a step of the simulation engine. Also stops the animation.
    PP.dir  = PLAY_STOP;
    if ( PP.live ) {
      PP.live = 0;
    } else {
      if ( sim.isReady() )
        sim.stepAll();
    }
    break;

  case 'S':   //a step of the simulation engine, starting it if necessary
    PP.dir  = PLAY_STOP;
    PP.live = 0;
    if ( NO_ERROR == prepareToSimulate() )
      sim.stepAll();
    else
      flashText("prepareToSimulate() failed");
    break;

  case 'b':
    PP.scalebar = ! PP.scalebar;
    //MSG.shutUp();
    return;

  case 'R':
    PP.showGrad = ! PP.showGrad;
    //MSG.shutUp();
    return;


  //case 'B':
    //MSG.increaseVerboseLevel();
    //return;

    /* temporarily disabled
  case 'w':
    appendToFile();
    return;

  case 'W':
    writeToFile();
    return;

  case 'r':
    readFromFile();
    break;
    */

  case 'g':
    flashText("Released mouse-controled grafted");
    playGrafted = 0;
    break;

  case 'G':
    flashText("Deleted mouse-controled grafted");
    for(int ii = 0; ii < playGraftedList.size(); ++ii) {
      if (playGraftedList[ii])
        delete( playGraftedList[ii] );
    }
    playGraftedList.clear();
    playGrafted = 0;
    break;

  //------------------------- play / stop / reverse:

  case 'p':
      if ( PP.dir != PLAY_FORWARD ) {
        //rewind if at the end of the file:
        if ( PP.frame >= reader.lastFrameBeforeEof() ) {
          PP.frame = 0;
        }
        PP.dir = PLAY_FORWARD;
      } else {
        PP.delay /= 2;
        //limit to 120 Hz max:
        if ( PP.delay < 8 )
          PP.delay = 8;
      }
      return;

  case 'o':
    PP.delay *= 2;
    return;

  case 'i':
    PP.dir  = PLAY_STOP;
    PP.live = 0;
    return;

  case 'u':
    if ( PP.dir != PLAY_REVERSE ) {
      PP.dir = PLAY_REVERSE;
      flashText("Play reverse");
    }
    else if ( PP.delay > 16 )
      PP.delay /= 2;
    return;

  case 'l':
    PP.loop = ! PP.loop;
    if ( PP.loop )
      flashText("Loop movie");
    return;

  case 'k':
    PP.frame = 0;
    break;

    //------------------------- zooming & window resizing:

  case '-':
  case '_':
    switch( glutGetModifiers() ) {
    default:
      PP.zoom /= 1.4142;
      flashText("Zoom %.2f", PP.zoom);
      break;
    case GLUT_ACTIVE_SHIFT:
      PP.zoom /= 1.04142;
      flashText("Zoom %.2f", PP.zoom);
      break;
    case GLUT_ACTIVE_ALT:
      if (( PP.width > 128 ) && ( PP.height > 128 ))
        glutReshapeWindow(PP.width-64, PP.height-64);
      break;
    }
    setModelView();
    break;

  case '=':
  case '+':
    switch( glutGetModifiers() ) {
    default:
      PP.zoom *= 1.4142;
      flashText("Zoom %.2f", PP.zoom);
      break;
    case GLUT_ACTIVE_SHIFT:
      PP.zoom *= 1.04142;
      flashText("Zoom %.2f", PP.zoom);
      break;
    case GLUT_ACTIVE_ALT:
      glutReshapeWindow(PP.width+64, PP.height+64);
      break;
    }
    setModelView();
    break;

    //------------------------- microtubules:

  case '1':
    if ( glutGetModifiers() & GLUT_ACTIVE_ALT ) {
      user_control_keys[1] = !user_control_keys[1]; ///< press ALT-1
      break;
    }

    //Rotate between 5 different MT display modes:
    if ( PP.mtlines ) {
      PP.mtlines = 0; PP.mtratchets = 1; PP.mtpoints = 0; PP.mtspeckles = 0; PP.mtrainbow = 0;
      flashText("microtubules:  ratchets");
    } else if ( PP.mtratchets ) {
      PP.mtlines = 0; PP.mtratchets = 0; PP.mtpoints = 1; PP.mtspeckles = 0; PP.mtrainbow = 0;
      flashText("microtubules:  model points");
    } else if ( PP.mtpoints ) {
      PP.mtlines = 0; PP.mtratchets = 0; PP.mtpoints = 0; PP.mtspeckles = 1; PP.mtrainbow = 0;
      flashText("microtubules:  speckles");
    } else if ( PP.mtspeckles ) {
      PP.mtlines = 0; PP.mtratchets = 0; PP.mtpoints = 0; PP.mtspeckles = 0; PP.mtrainbow = 1;
      flashText("microtubules:  rainbow pressure");
    } else {
      PP.mtlines = 1; PP.mtratchets = 0; PP.mtpoints = 0; PP.mtspeckles = 0; PP.mtrainbow = 0;
      flashText("microtubules:  lines");
    }

    break;

    // pressing the '!' key cycles between different display modes:
    // showing the plus ends -> the minus ends -> none
  case '!':
    if ( PP.mtends[MT_PLUS_END] ) {
      PP.mtends[MT_MINUS_END] = 1;
      PP.mtends[MT_PLUS_END]  = 0;
      flashText("MT minus-ends");
    } else if ( PP.mtends[MT_MINUS_END] ) {
      PP.mtends[MT_MINUS_END] = 0;
      PP.mtends[MT_PLUS_END]  = 0;
      flashText("MT no ends");
    } else {
      PP.mtends[MT_PLUS_END]  = 1;
      flashText("MT plus-ends");
    }
    break;

  case '2':
    if ( glutGetModifiers() & GLUT_ACTIVE_ALT ) {
      user_control_keys[2] = !user_control_keys[2]; ///< press ALT-2
      break;
    }

    if ( PP.mtwidth > SIZE_INC )  PP.mtwidth -= SIZE_INC;
    flashText("MT line width %.1f", PP.mtwidth);
    break;

  case '3':
    if ( glutGetModifiers() & GLUT_ACTIVE_ALT ) {
      user_control_keys[3] = !user_control_keys[3]; ///< press ALT-3
      break;
    }

    PP.mtwidth += SIZE_INC;
    flashText("MT line width %.1f", PP.mtwidth);
    break;

  case '@':
    if ( PP.mtsize > SIZE_INC )  PP.mtsize -= SIZE_INC;
    flashText("MT point size %.1f", PP.mtsize);
    break;

  case '#':
    PP.mtsize += SIZE_INC;
    flashText("MT point size %.1f", PP.mtsize);
    if ( !PP.mtends[MT_MINUS_END] && !PP.mtends[MT_PLUS_END] && !PP.mtspeckles )
      PP.mtpoints = 1;
    break;

  //-------------------------
  case '4':
    if ( glutGetModifiers() & GLUT_ACTIVE_ALT ) {
      user_control_keys[4] = !user_control_keys[4]; ///< press ALT-4
      break;
    }
    break;

  //------------------------ points size for solids:

  case '%':
    PP.solinks = ! PP.solinks;
    break;

  case '^':
    if ( PP.sosize > SIZE_INC )  PP.sosize -= SIZE_INC;
    flashText("Solids: size %.1f", PP.sosize);
    break;

  case '&':
    PP.sosize += SIZE_INC;
    flashText("Solids: point size %.1f", PP.sosize);
    PP.sopoints = 1;
    break;

    //------------------------ points size for grafted:

  case '5':
    switch( PP.ghselect ) {
      case 255: PP.ghselect = 2;   flashText("Grafteds: only bound"); break;
      case 2:   PP.ghselect = 1;   flashText("Grafteds: only free");  break;
      default:  PP.ghselect = 255; flashText("Grafteds: all");        break;
    }
    break;

  case '6':
    if ( PP.ghsize > SIZE_INC )  PP.ghsize -= SIZE_INC;
    flashText("Grafted: size %.1f", PP.ghsize);
    break;

  case '7':
    PP.ghsize += SIZE_INC;
    flashText("Grafted: size %.1f", PP.ghsize);
    PP.ghhands = 1;
    break;


    //------------------------ points size for complexes:

  case '8':
    switch( PP.cxselect ) {
      case 255: PP.cxselect = 4;   flashText("Complexes: only bridge"); break;
      case 4:   PP.cxselect = 2;   flashText("Complexes: only bound");  break;
      case 2:   PP.cxselect = 1;   flashText("Complexes: only free");   break;
      case 1:   PP.cxselect = 255; flashText("Complexes: all");         break;
    }
    break;

  case '9':
    if ( PP.cxsize > SIZE_INC )  PP.cxsize -= SIZE_INC;
    flashText("Complexes: size %.1f", PP.cxsize);
    break;

  case '0':
    if ( glutGetModifiers() & GLUT_ACTIVE_ALT ) {
      user_control_keys[0] = !user_control_keys[0]; ///< press ALT-0
      break;
    }

    PP.cxsize += SIZE_INC;
    flashText("Complexes: size %.1f", PP.cxsize);
    PP.cxhands = 1;
    break;

    //----------------- link width for both Complex & Grafted links

  case '(':
    if ( PP.cxwidth > SIZE_INC )  PP.cxwidth -= SIZE_INC;
    if ( PP.ghwidth > SIZE_INC )  PP.ghwidth -= SIZE_INC;
    flashText("Complexes: link width %.1f", PP.cxwidth);
    break;

  case ')':
    PP.cxwidth += SIZE_INC;
    PP.ghwidth += SIZE_INC;
    flashText("Complexes: link width %.1f", PP.cxwidth);
    PP.cxlinks = 1;
    break;

  default:
    flashText("ignored key [%c]", c);
    //MSG("unknown command key <%c> %i %i\n", c, x, y);
    return;
  }

  //if break was called, we call for a redisplay:
  glutPostRedisplay();
  buildMenus();
}


//--------------------------------------------------------------------------------------
void Player::processInputKey(const int key, const int x, const int y)
{
  //the amount of motion by one key-press, this is relative to
  //what is currently being displayed on the window (i.e. depends on Zoom)
  real percent = 0.2;

  // if shift is down, we make smaller moves:
  if ( glutGetModifiers() == GLUT_ACTIVE_SHIFT )
    percent /= 20;

  // we move the focus point according to the arrow key pressed:
  switch ( key ) {
    case GLUT_KEY_LEFT:  PP.focus[0] += percent * visibleRegion[0] / PP.zoom;  break;
    case GLUT_KEY_RIGHT: PP.focus[0] -= percent * visibleRegion[0] / PP.zoom;  break;
    case GLUT_KEY_UP:    PP.focus[1] -= percent * visibleRegion[1] / PP.zoom;  break;
    case GLUT_KEY_DOWN:  PP.focus[1] += percent * visibleRegion[1] / PP.zoom;  break;
    default:
      //MSG("ignored input key %i mouse=( %i %i )\n", key, x, y);
      return;
  }

  flashText("Focus %.2f %.2f %.2f", PP.focus[0], PP.focus[1], PP.focus[2]);
  setModelView(1);
}

