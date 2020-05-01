//RCS: $Id: testglut2d.cc,v 2.5 2005/03/30 12:31:10 nedelec Exp $
//--------------------------------------------------------------------------
//               simplest GLUT interface with 2D graphics
//
//      Francois Nedelec nedelec@embl.de  December 2003
// compilation on mac:
// g++ testglut2d.cc -framework GLUT -framework openGL -framework Foundation
// Linux:
// g++ testglut2d.cc -L/usr/X11R6/lib -lglut -lGL -lGLU -lXt -lX11 -lXi -lXmu

// This program displays a triangle with a blinking dot. You can move the dot
// by clicking the left button, and you can zoom in/out by draging with the
// mouse with the middle button down. The blinking uses a timer callback.

#include <cstdlib>
#include <cstdio>
#include <cmath>

#ifdef __APPLE__
   #include <GLUT/glut.h>
#else
   #include <GL/glut.h>    //use this on Linux & PC
#endif

//------------- precision: real is an alias for double or float:

typedef double real;

//------------- initial size of the visible region:

const real vsize[] = { 1.0, 1.0 };

//------------- size of the display window:

GLint viewport[4];

//------------- delay for the Timer function in milli-seconds:

int delay = 30;

//------------- user adjustable zoom:

real zoom = 0.8;

//------------- variables for the mouse driven zoom:

int mouse_action;
real zoom_save, mouse_zoom_scalar;

//position of click:
real click[] = { 0, 0, 0 };
//size of point drawn:
int point_size = 1;

//------------- matrices for the unprojection by GLU:

GLdouble modelViewMatrix[16];
GLdouble projectionMatrix[16];

//------------- function to set the OpenGL transformation:

void setModelView(int redisplay = 0)
{
    glMatrixMode(GL_MODELVIEW);
 
    //set the matrix with a simple zoom:
    glLoadIdentity();
    glScaled( zoom, zoom, zoom );
    
    //get the modelView matrix:
    glGetDoublev(GL_MODELVIEW_MATRIX, modelViewMatrix);

    if ( redisplay ) glutPostRedisplay();
}


//--------------- standard function to check for OpenGL errors:

void glPrintErrors(char * where)
{
  GLenum glError;
  do {
    glError = glGetError();
    switch( glError ) {
      case GL_NO_ERROR:  
        return;
      case GL_INVALID_ENUM:
        printf("OpenGL error GL_INVALID_ENUM in %s\n", where); 
        break;
      case GL_INVALID_VALUE: 
        printf("OpenGL error GL_INVALID_VALUE in %s\n", where); 
        break;
      case GL_INVALID_OPERATION: 
        printf("OpenGL error GL_INVALID_OPERATION in %s\n", where);
        break;
      case GL_STACK_OVERFLOW:
        printf("OpenGL error GL_STACK_OVERFLOW in %s\n", where); 
        break;
      case GL_STACK_UNDERFLOW:
        printf("OpenGL error GL_STACK_UNDERFLOW in %s\n", where);
        break;
      case GL_OUT_OF_MEMORY:
        printf("OpenGL error GL_OUT_OF_MEMORY in %s\n", where); 
        break;
    }
  } while ( glError != GL_NO_ERROR );
}


//----------------------------- KEYS --------------------------------

void processNormalKey(unsigned char c, int mouse_x, int mouse_y)
{
  switch (c) {
	
  case 27:
  case 'q':
    exit(EXIT_SUCCESS);
    
  case ' ':
    zoom = 1;
    break;   //we use break to call the setModelView() below
    
  default:
    printf("unknown normal key %c\n", c);
  }
  setModelView(1);
}


void processInputKey(int c, int mouse_x, int mouse_y)
//these are special keys: arrows, ctrl, etc.
{
  printf("unknown special key %c\n", c);
}


//----------------------------- MENUS --------------------------------

enum MENUS_ID { MENU_QUIT, MENU_RESETVIEW };

void processMenu(int item)
{
  switch ( item ) {
  case MENU_QUIT: 
    exit( EXIT_SUCCESS );
  case MENU_RESETVIEW:
    zoom = 1;
    setModelView(1);
    break;
  }
}

void initMenus()
{
    // --- one menu with two entries:
  glutCreateMenu(processMenu);
  glutAddMenuEntry("Reset", MENU_RESETVIEW);
  glutAddMenuEntry("Quit", MENU_QUIT);
  glutAttachMenu(GLUT_RIGHT_BUTTON);
}

//----------------------------- MOUSE -------------------------------

//----------------list of different actions control by the mouse

enum MOUSE_ACTION { MOUSE_NOTHING, MOUSE_ZOOM, MOUSE_CLICK };


void processMouse(int button, int state, int x, int y)
  //this is called when the mouse button is pressed or released:
{
  mouse_action = MOUSE_NOTHING; 
    
  //for button release, we do nothing:
  if ( state != GLUT_DOWN ) return;   
  
  //choice of action depending on the button pressed:
  switch( button ) {
  case GLUT_LEFT_BUTTON:
    mouse_action = MOUSE_CLICK;
    break;
  case GLUT_MIDDLE_BUTTON:
    mouse_action = MOUSE_ZOOM;
    break;
  }
  
  //depending on the action...
  switch( mouse_action ) {
    
  case MOUSE_CLICK: {
    // --- we unproject the mouse to get the natural coordinates:
    gluUnProject(x, viewport[3]-y, 0, modelViewMatrix, projectionMatrix,
		 viewport, click, click+1, click+2);
    point_size = 0;
    setModelView(1);
  } break;
  
  case MOUSE_ZOOM: {
    real xx = x - viewport[2]/2.0;
    real yy = y - viewport[3]/2.0;
    mouse_zoom_scalar = sqrt( xx*xx + yy*yy );
    if ( mouse_zoom_scalar > 0 ) 
      mouse_zoom_scalar = 1.0 / mouse_zoom_scalar;
    zoom_save = zoom;
  } break;
  
  case MOUSE_NOTHING:
    return;
  }
}


void processMotion(int x, int y)
{
  switch( mouse_action ) {
    
  case MOUSE_CLICK: {
    // --- unproject the mouse, to get the natural coordinates:
    gluUnProject(x, viewport[3]-y, 0, modelViewMatrix, projectionMatrix, 
		 viewport, click, click+1, click+2);
  } break;
  
  case MOUSE_ZOOM: {
    // --- we set the zoom from how far the mouse is from the window center
    real xx = x - viewport[2]/2.0;
    real yy = y - viewport[3]/2.0;
    real Z = mouse_zoom_scalar * sqrt( xx*xx + yy*yy );
    if ( Z <= 0 ) return;
    zoom = zoom_save * Z;
  } break;
  
  case MOUSE_NOTHING:
    return;
  }
  setModelView(1);
}

//----------------------------- RESHAPE --------------------------------


void reshapeWindow(int w, int h)
{
  glViewport(0, 0, w, h);
  
  // --- get the viewport, needed for the unprojection:
  glGetIntegerv(GL_VIEWPORT, viewport);
  
  // --- set-up the projection matrix:
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  real ratio = w * vsize[1] / double( vsize[0] * h );
  
  if ( ratio > 1 )
    glOrtho(-vsize[0], vsize[0], -vsize[1]/ratio, vsize[1]/ratio, 0, 1 );
  else
    glOrtho(-vsize[0]*ratio, vsize[0]*ratio, -vsize[1], vsize[1], 0, 1 );    
  
  // --- get the projection matrix needed for the unprojection
  glGetDoublev(GL_PROJECTION_MATRIX, projectionMatrix);
}

//----------------------------- DISPLAY --------------------------------

void display()
{
  // --- clear window to the current clearing color:
  glClear( GL_COLOR_BUFFER_BIT );
  
  // --- set line width to 1 and color to white:
  glColor3f(1.0, 1.0, 1.0);
  glLineWidth(1.0);
  
  // --- draw a wireframe triangle:
  glBegin(GL_LINE_LOOP);
  glVertex2f( 1.0, -1.0);
  glVertex2f(-1.0, -1.0);
  glVertex2f( 0.0,  1.0);
  glEnd();
  
  // --- a point of variable size at the last clicked position:
  glColor3f(1.0, 1.0, 0.0);
  glPointSize( 1 + point_size );
  glBegin(GL_POINTS);
  glVertex2f(click[0], click[1]);
  glEnd();
  
  glFinish();
  glutSwapBuffers();
  // --- check for any OpenGL error:
  glPrintErrors("display()");  
}

//------------------------- TIMER FUNCTION -----------------------------
void timerFunction(int value)
{
  point_size = ( point_size + 1 ) % 16;
  glutPostRedisplay();
  //register another timer call back in PP.delay milli-sec:
  glutTimerFunc(delay, timerFunction, 1);
}


//----------------------------- INIT GL --------------------------------
void initGL()
{
  // --- choose the clearing color: black
  glClearColor(0.0, 0.0, 0.0, 0.0);
  
  //--- hints for OpenGL rendering:
  glEnable(GL_BLEND);
  glBlendFunc( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA );
  
  glEnable(GL_POINT_SMOOTH);
  glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);
  
  glEnable(GL_LINE_SMOOTH);
  glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
  
  initMenus();
  setModelView();
}


//-----------------------------   MAIN  --------------------------------
int main(int argc, char *argv[])
{
  
  // --- initialization of GLUT:
  glutInit(&argc, argv);
  glutInitDisplayMode( GLUT_RGBA | GLUT_DOUBLE );
  glutInitWindowSize(400, 400); 
  glutInitWindowPosition(50, 50);
  glutCreateWindow(argv[0]);
  
  // --- further initialization of OpenGL:
  initGL();
  
  // --- register all the necessary functions:
  glutDisplayFunc(display);
  glutReshapeFunc(reshapeWindow);
  glutMouseFunc(processMouse);
  glutMotionFunc(processMotion);
  glutSpecialFunc(processInputKey);
  glutKeyboardFunc(processNormalKey);
  glutTimerFunc(50, timerFunction, 0);
  
  // --- starts the event loop, which will never return:
  glutMainLoop();
  
  return EXIT_SUCCESS;
}
