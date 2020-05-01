//RCS: $Id: play.h,v 1.5 2005/03/30 21:01:40 nedelec Exp $

#ifndef PLAY_H
#define PLAY_H

#include "types.h"
#include "sim.h"
#include "glut.h"
#include "quaternion.h"
#include "vecteur3.h"
#include "reader.h"
#include "array.h"


///Player class
namespace Player {

  ///reader use to access result.out file
  SimReader reader;

  ///amount added/subtracted to size/width when a key is pressed
  static const real SIZE_INC = 0.5;

  //==========================GLOBAL VARIABLES=========================

  ///size of the OpenGL visible region
  real visibleRegion[3] = { 100, 100, 100 };

  ///translation of camera with respect to the center of the region
  real cameraTranslation[3] = { 0, 0, -100 };

  ///backup copy of the OpenGL transformation used in unproject()
  GLint    viewport[4];
  ///backup copy of the OpenGL transformation used in unproject()
  GLdouble modelviewMat[16];
  ///backup copy of the OpenGL transformation used in unproject()
  GLdouble projectionMat[16];

  /// list of Grafted used to grab microtubules with the mouse
  Array<Grafted *> playGraftedList;
  /// the current Grafted being controlled with the mouse
  Grafted * playGrafted = 0;
  /// the distance to which a grafted can be controled
  static const real playGraftedDistanceThreshold = 0.5;

  ///amplification factor for mouse controlled rotation
  static const real mouse_amplification = 3.0;

  ///the current action being performed by the mouse
  int mouse_action;

  //----------------variables for the mouse driven rotation/transformation
  real       zoom_save;           ///<saved value for the scaling factor
  Quaternion quat_save;           ///<saved value for the rotation quaternion
  Vecteur3   focus_save;          ///<saved value for focus point

  Vecteur3 mouse_origin;          ///<origin for mouse drag action
  Vecteur3 mouse_normal;          ///<vector normal to desired rotation
  Vecteur3 mouse_axis;            ///<axis of rotation for MOUSE_SPIN
  real     mouse_zoom_scalar;     ///<accessory scalar for zooming


  //--------------------------------- MENUS ----------------------------------

  ///list of GLUT menus handles
  Array<int> playMenus;

  ///function to store a menu handle
  int storeMenu(const int menu) { playMenus.pushFront(menu); return menu; }

  ///clear all menu handles
  void eraseAllMenus();

  ///build all the menus from scratch
  void buildMenus();

  ///menu callback function
  void processMenuEvents(const int menu);

  //-----------------------------------KEYS------------------------------------

  ///returns a string with some help on what pressing key do
  const char * descriptionOfKeys();

  ///GLUT callback function for most keys
  void processNormalKey(const unsigned char, const int = 0, const int = 0 );

  ///GLUT callback function for arrow keys
  void processInputKey(const int key, const int x, const int y);

  //-----------------------------------MOUSE-----------------------------------

  ///list of different actions controlled by the mouse
  enum MOUSE_ACTION { MOUSE_NOTHING, MOUSE_TRANSLATE, MOUSE_TRANSLATE_AXIS,
                      MOUSE_ZOOM, MOUSE_ROTATE, MOUSE_SPIN, MOUSE_GRAB };

  ///the buttons under which menu pops up
  static const int MENU_BUTTON = GLUT_RIGHT_BUTTON;

  ///this function maps the button pressed to the MOUSE_ACTION
  MOUSE_ACTION mouseAction(const int button, const int modifier);

  ///GLUT callback function for mouse button pressed
  void processMouse(const int button, const int state, const int x, const int y);

  ///GLUT callback function for mouse motion
  void processMotion(const int x, const int y);

  //----------------------------------GRAPHICS---------------------------------

  ///init OpenGL for PP.style=1
  void initGL1();
  ///display function for PP.style=1
  void drawFrame1();

  ///init OpenGL for PP.style=2
  void initGL2();
  ///display function for PP.style=2
  void drawFrame2();

  ///init OpenGL for PP.style=3
  void initGL3();
  ///display function for PP.style=3
  void drawFrame3();

  ///init function for any style, calls one of the initGL?()
  void initGL();
  ///display function for any style, calls one of the drawFrame?()
  void drawFrame();

  //---------------------------------------------------------------------------

  ///function to show text on the display, for a short period of time
  void flashText(const char fmt[], ...);

  ///set the model-view transformation matrix from zoom and quat
  void setModelView(const bool update = false);

  ///reset the view (rotation+zoom)
  void resetView();

  ///GLUT callback function for window resize event
  void reshapeWindow(const int, const int);

  ///GLUT callback function for timed events
  void timerFunction(const int);

  ///swap in/out of full-screen mode
  void swapFullScreen();


  ///repeat the display, for periodic boundary conditions
  void tileFrame();
  ///display a portion of a scale bar
  void displaySubScaleBar(const int, const real, const real);
  ///display a scale bar
  void displayScaleBar(const int);
  ///display gradient
  void displayGradient(const int);
  ///display the message on the window: time, frame nb, etc.
  void displayMessage();
  ///display the parameter values used in the live mode
  void displayParameters(const bool displayValues = true);
  ///GLUT display function
  void display();

  //------------------------------------I/O------------------------------------

  ///go to the next frame in result.out
  void prevFrame();
  ///go to the previous frame in result.out
  void nextFrame();
  ///function called after a frame is read, set the Microtub colors
  void processState();


  ///save the displayed imaged in a graphic file, PPM format
  int  saveImageAsPPM();
  ///save all frame stored in result.out as graphic files
  void saveMovieAsPPM();

  //---------------------------------------------------------------------------

  ///write current state of sim to a file
  void writeToFile();
  ///append current state of sim to a file
  void appendToFile();
  ///read sim-state from a file
  void readFromFile();

  //---------------------------------------------------------------------------

  ///to be called when switching to live simulation mode
  int  prepareToSimulate(const bool forceIt = false);
  ///reset the sim-state and timer
  void resetSimulation();

  ///the different values for the play-mode in PP.dir
  enum PlayMode { PLAY_STOP=0, PLAY_FORWARD=1, PLAY_REVERSE=-1, PLAY_FORWARD_WRITE=2 };
};


#endif
