//RCS: $Id: glut_base.cc,v 2.5 2005/03/30 21:04:07 nedelec Exp $
//--------------------------------------------------------------------------
//                               glut_base.cc
//      mouse driven zoom and rotation with quaternions and GLUT unproject
//
//   Francois Nedelec nedelec@embl.de  October 2002, modified Nov. 2003

#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "glut.h"
#include "vecteur3.h"
#include "quaternion.h"

//an amplification factor to make the rotation more responsive
const real mouse_amplification = 3;

//size of viewing box in the real-world units:
Vecteur3 visibleRegion( 1.0, 1.0, 1.0 );

//the translation of camera
Vecteur3 cameraTranslation( 0, 0, -3.0 );

//point of focus
Vecteur3 focus_save, focus( 0, 0, 0 );

//the zoom factor, and saved value:
real zoom_save, zoom = 1;

//the quaternion which defines the rotation, and a save value:
//quat is initialized to Id (no rotation)
Quaternion quat_save, quat(1);

//----------------list of different actions controled by the mouse

enum MOUSE_ACTION { MOUSE_NOTHING, MOUSE_TRANSLATE, MOUSE_TRANSLATE_AXIS,
                    MOUSE_ZOOM,    MOUSE_ROTATE,    MOUSE_SPIN };

//----------------the attribution of mouse actions
const int MENU_BUTTON     = GLUT_RIGHT_BUTTON;

//----------------attributions of the actions with mouse + option-keys
MOUSE_ACTION mouseAction( int button, int modifier ) {
  switch( button ) {
  case GLUT_LEFT_BUTTON:
    switch( modifier ) {
    case GLUT_ACTIVE_CTRL:
      return MOUSE_SPIN;
    case GLUT_ACTIVE_SHIFT:
      return MOUSE_TRANSLATE;
    case GLUT_ACTIVE_SHIFT + GLUT_ACTIVE_CTRL:
      return MOUSE_TRANSLATE_AXIS;
    case 0:
      return MOUSE_ROTATE;
    }
  case GLUT_MIDDLE_BUTTON:
    return MOUSE_ZOOM;
  }
  return MOUSE_NOTHING;
}

//------------ variables used for the mouse zoom/rotation
int mouse_action;

Vecteur3 mouse_origin, mouse_normal, mouse_axis;
real mouse_zoom_scalar;

// copies of the OpenGL transformation matrices used in glutUnproject()
GLdouble modelviewMat[16];
GLdouble projectionMat[16];
GLint viewport[4];


//----------------------------------------------------------------------------
// OpenGL debug function:
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
    default:
      printf("OpentGL unknow error in %s\n", where);
      break;
    }
  } while ( glError != GL_NO_ERROR );
}


//----------------------------------------------------------------------------
// set the Model-view transformation
void setModelView()
{
  //the opengl matrix of model-view transformation
  static GLdouble mat[16];

  glMatrixMode(GL_MODELVIEW);

  quat.setThisMatrix16( mat, cameraTranslation );
  glLoadMatrixd( mat );
  glScaled( zoom, zoom, zoom );
  glTranslated( -focus[0], -focus[1], -focus[2]);

  glutPostRedisplay();
}


//----------------------------------------------------------------------------
void resetModelView(real zoomset = 1.0)
{
  quat.set(1);
  focus.set(0,0,0);
  zoom = zoomset;
  setModelView();
}

//----------------------------------------------------------------------------
// window reshape callback-function
void reshapeWindow(int w, int h)
{
  glViewport(0, 0, w, h);
  glGetIntegerv(GL_VIEWPORT, viewport);

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
  
  glGetDoublev(GL_PROJECTION_MATRIX, projectionMat);

#ifdef USE_OPENGL_FOG
  glFogf(GL_FOG_DENSITY, 0.5/visibleRegion[2]);
  glFogf(GL_FOG_START, znear);
  glFogf(GL_FOG_END,   zfar);
#endif
  glMatrixMode(GL_MODELVIEW);
}

//----------------------------------------------------------------------------

void processMouse(const int button, const int state, const int x, const int y)
{
  if ( state != GLUT_DOWN ) {
    mouse_action = MOUSE_NOTHING; 
    return; 
  }

  mouse_action = mouseAction( button, glutGetModifiers() );
  if ( mouse_action == MOUSE_NOTHING ) 
    return;

  glGetDoublev(GL_MODELVIEW_MATRIX, modelviewMat);

  GLdouble dx, dy, dz;
  gluUnProject(x, viewport[3]-y, 0, modelviewMat,
			   projectionMat, viewport, &dx, &dy, &dz);

  switch( mouse_action ) {

  case MOUSE_TRANSLATE: {
    mouse_origin.set( dx, dy, dz );
    focus_save = focus;
  } break;
  
  case MOUSE_TRANSLATE_AXIS: {
    mouse_origin.set( dx, dy, dz );
    GLdouble cx, cy, cz;
    gluUnProject(viewport[2]/2.0, viewport[3]/2.0, 0, modelviewMat,
		 projectionMat, viewport, &cx, &cy, &cz);
    mouse_normal = ( Vecteur3( cx, cy, cz ) - focus ).normalized();
    gluUnProject(viewport[2]/2.0, viewport[3], 0, modelviewMat,
		 projectionMat, viewport, &dx, &dy, &dz);
    mouse_axis = Vecteur3( dx-cx, dy-cy, dz-cz ).normalized();
    focus_save = focus;
  } break;
  
  case MOUSE_ROTATE: {
    mouse_origin.set( dx, dy, dz );
    mouse_normal = mouse_origin - focus;
    if ( mouse_normal.normSquare() )
      mouse_normal *= mouse_amplification / mouse_normal.normSquare();
    quat_save = quat;
  } break;
  
  case MOUSE_SPIN: {
    GLdouble cx, cy, cz;
    gluUnProject(viewport[2]/2.0, viewport[3]/2.0, 0, modelviewMat, 
		 projectionMat, viewport, &cx, &cy, &cz);
    mouse_origin.set( cx, cy, cz );
    mouse_axis = ( mouse_origin - focus ).normalized();
    mouse_normal = Vecteur3( dx, dy, dz ) - mouse_origin;
    quat_save = quat;
  } break;
  
  case MOUSE_ZOOM: {
    real xx = x - viewport[2]/2.0;
    real yy = y - viewport[3]/2.0;
    mouse_zoom_scalar = sqrt( xx*xx + yy*yy );
    if ( mouse_zoom_scalar > 0 ) 
      mouse_zoom_scalar = 1.0 / mouse_zoom_scalar;
    zoom_save = zoom;
  } break;
  }
}


//----------------------------------------------------------------------------
void processMotion(int x, int y)
{
  if ( mouse_action == MOUSE_NOTHING ) 
    return;

  GLdouble dx, dy, dz;
  gluUnProject(x, viewport[3]-y, 0, modelviewMat,
			   projectionMat, viewport, &dx, &dy, &dz);
  Vecteur3 mouse_drag = Vecteur3( dx, dy, dz ) - mouse_origin;
  
  switch( mouse_action ) {

  case MOUSE_ROTATE: {
    quat.setFromAxisAngle( mouse_normal ^ mouse_drag );
    quat.leftMult( quat_save );
    setModelView();
  } break;
  
  case MOUSE_SPIN: {
    real cos = mouse_normal * mouse_drag;
    real sin = ( mouse_normal ^ mouse_drag ) * mouse_axis;
    quat.setFromAxisAngle( mouse_axis, atan2( sin, cos ) ); 
    quat.leftMult( quat_save ); 
    setModelView();
  } break;
  
  case MOUSE_TRANSLATE: {
    focus = focus_save - mouse_drag;
    setModelView();
  } break;
  
  case MOUSE_TRANSLATE_AXIS: {
    real S = mouse_drag * mouse_axis;
    focus = focus_save + S * ( mouse_normal + mouse_axis ) - mouse_drag;
    setModelView();
  } break;
  
  case MOUSE_ZOOM: {
    real xx = x - viewport[2]/2.0;
    real yy = y - viewport[3]/2.0;
    real Z = mouse_zoom_scalar * sqrt( xx*xx + yy*yy );
    if ( Z > 0 ) {
      zoom = zoom_save * Z;
      setModelView();
    }
  } break;
  }
}

//---------------------------------------------------------------------------
void displayText(const char text[], const int position = 0 )
{
  int line = 0;
  GLint shift, raster_position[2];
    
  switch( position ) {
    case 0: //bottom-left
      raster_position[0] = 10;
      raster_position[1] = 10;
      shift = +13;
      break;
    case 1: //bottom-right
      raster_position[0] = viewport[2] - 8*strlen(text) -10;
      raster_position[1] = 10;
      shift = +13;
      break;
    case 2: //top-right
      raster_position[0] = viewport[2] - 8*strlen(text) -10;
      raster_position[1] = viewport[3] - 10;
      shift = -13;
      break;
    default:
    case 3: //top-left
      raster_position[0] = 10;
      raster_position[1] = viewport[3] - 10;
      shift = -13;
      break;
  }

  //set the matrices
  glMatrixMode(GL_MODELVIEW);
  glPushMatrix();
  glLoadIdentity();
  glMatrixMode(GL_PROJECTION);
  glPushMatrix();
  glLoadIdentity();
  glOrtho( 0, viewport[2], 0, viewport[3], -1, 1 );
  
  // set color and position of text
  glColor3f(1.0, 1.0, 1.0);
  glRasterPos2i(raster_position[0], raster_position[1] + shift * line);

  // draw the string character per character in a small font:
  for(const char * p = text; *p; ++p) {
    if ( *p == '\n' )
      glRasterPos2i(raster_position[0], raster_position[1] + shift * (++line));
    else
      glutBitmapCharacter(GLUT_BITMAP_8_BY_13, *p);
  }
  glPopMatrix();
  glMatrixMode(GL_MODELVIEW);
  glPopMatrix();
}


//----------------------------------------------------------------------------
#ifdef DEFINE_MENUS

enum MENUS_ID { MENU_QUIT, MENU_RESETVIEW };

void processMenu(int item)
{
  switch( item ) {
    case MENU_QUIT:
      exit( EXIT_SUCCESS );
    case MENU_RESETVIEW:    
      resetModelView();
      break;
  }
}

void initMenus()
{
  glutCreateMenu(processMenu);
  glutAddMenuEntry("Reset", MENU_RESETVIEW);
  glutAddMenuEntry("Quit", MENU_QUIT);
  glutAttachMenu(MENU_BUTTON);
}

#endif

//----------------------------------------------------------------------------
void initGLUT_default()
{
  glClearColor(0.0, 0.0, 0.0, 0.0);
  
  glEnable(GL_BLEND);
  glBlendFunc( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA );
  
  glEnable(GL_POINT_SMOOTH);
  glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);
  
  glEnable(GL_LINE_SMOOTH);
  glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
    
#ifdef USE_OPENGL_FOG
  glEnable(GL_FOG);
  glFogi(GL_FOG_MODE, GL_LINEAR);
  GLfloat fogColor[] = { 0.0, 0.0, 0.0, 1.0 };
  glFogfv(GL_FOG_COLOR, fogColor);
#endif
  
#ifdef DEFINE_MENUS
  initMenus();
#endif
  
}


