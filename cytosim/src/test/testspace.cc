//RCS: $Id: testspace.cc,v 2.15 2005/04/08 10:09:39 foethke Exp $
// Visual test for the rasterizer used in attachment algorithm of Cytosim
// Francois Nedelec, nedelec@embl.de, October 2002

#if ( DIM == 3 )
  #define USE_OPENGL_FOG
#endif
#include "glut_base.cc"

#include <ctime>

#include "types.h"
#include "smath.h"
#include "random.h"
#include "space.h"
#include "vecteur.h"

#include "glextensions.h"
using namespace glExtensions;

//the space to be tested:
Space * space;

//number of points
const int maxpts = 256000;
      int nbpts  = 1000;

//INFLATION of the rectangle containing point to be projected
const real INFLATION = 10;

//
real boxsize_set[5] = { 5, 10, 5, 2, 2 }; 


//inflation factor
int  inflation_mode = 0;
real inflation[] = { 0, 1, 2, 5, 10, 0, -1, -2, -5, -10 };
real inflation_radius = 0;

//shape for the space:
int spaceShape = SHAPE_SQUARE;

//coordinates of the points:
Vecteur      points[ maxpts ];

//true if inside
int          inside[ maxpts ];

//coordinates of the projections
Vecteur projections[ maxpts ];

//coordinates of the projections of the projections
Vecteur projections2[ maxpts ];

//max distance from projection to second projection
real  max_error_projection;

//slicing parameters
const real sliceStep = 0.2;
      real thickness = sliceStep;
      real slicePos  = 0;

//slicing is 1 for slicing in x,
//           2 for slicing in y,
//           4 for slicing in z
//or combinations of these
int  slicing   = false;

//show or hide points in or outside
int showInside  = true;
int showOutside = true;
int showLines   = true;
int showProject = true;

//use timer function on or off
int timer = false;
int delay = 50;

//display parameter for OpenGL
int line_width = 2;

//amount of white added to colors
const float COL = 0.8;
//----------------------------------------------------------------------------
void setPoints() {
  max_error_projection = 0;
  for( int ii = 0; ii < nbpts; ++ii ) {
    
    points[ ii ] = space->getBoundingRect();
    
    //we inflate the bounding rect in every direction
    points[ ii ].XX += INFLATION;
    points[ ii ].YY += INFLATION;
    points[ ii ].ZZ += INFLATION;
    
    //pick a random point withing this space:
    points[ ii ].XX *= RNG.sreal();
    points[ ii ].YY *= RNG.sreal();
    points[ ii ].ZZ *= RNG.sreal();

    //see if space finds it inside:
    inside[ ii ] = space->isInside( points[ii] );
    
    //calculate the projection:
    space->project( points[ii], projections[ii] );
    
    //calculate the projection of the projection:
    space->project( projections[ii], projections2[ii] );
    
    real d = (projections[ii] - projections2[ii]).normSquare();
    if ( d > max_error_projection ) max_error_projection = d;
  }
  max_error_projection = sqrt( max_error_projection );
}

//----------------------------------------------------------------------------
void timerFunction(int value)
{
  if( timer ) {
    setPoints();
    glutPostRedisplay();
  }
  
  glutTimerFunc( delay, timerFunction, timer );
}

//----------------------------------------------------------------------------
//set the space
void setSpace(int shape = SHAPE_SQUARE) 
{
  spaceShape = shape;
  real boxsize[8];
  
  switch( shape ) {
    case SHAPE_SQUARE:
    case SHAPE_PERIODIC:
    case SHAPE_ELLIPSE:
    case SHAPE_STRIP:
    case SHAPE_ROUND_SQUARE:
    case SHAPE_OVAL:
    case SHAPE_SPHERE:
    case SHAPE_CYLINDER:
    case SHAPE_CYLINDER_Z:
      for( int ii = 0; ii < 5; ++ii )
        boxsize[ii] = boxsize_set[ii];
      break;
      
    case SHAPE_OUTSIDE_SPHERE:
      boxsize[0] = 9;
      boxsize[1] = 10;
      boxsize[2] = 11;
      boxsize[4] = 7;
      break;
    case SHAPE_BANANA:
      boxsize[0] = 15;
      boxsize[1] = 4;
      boxsize[2] = 1.25;
      break;
    case SHAPE_BOOMERANG:
      boxsize[0] = 15;
      boxsize[1] = 4;
      boxsize[2] = 0.5;
      break;
    case SHAPE_TEE:
      boxsize[0] = 15;
      boxsize[1] = 4;
      boxsize[2] = 5;
      boxsize[3] = 8;
      break;
  }
  
  space =  newSpace( spaceShape, boxsize, inflation_radius);
  setPoints();
  glutPostRedisplay();
}

//----------------------------------------------------------------------------
enum MENUS_ID { MENU_QUIT = 102, MENU_RESETVIEW = 103,
                MENU_INSIDE = 104, MENU_OUTSIDE = 105, MENU_LINES = 106,
                MENU_XSLICING = 107, MENU_YSLICING = 108, MENU_ZSLICING = 109,
                MENU_INFLATION_0 = 110 };

void processMenu(int item)
{
  switch( item ) {
    case MENU_INFLATION_0:
      inflation_radius = 0;
      setSpace( spaceShape );
      break;
      
    case MENU_QUIT:
      exit( EXIT_SUCCESS );
    case MENU_RESETVIEW:  
      resetModelView(0.04);
      break;
    case MENU_INSIDE:
      showInside = ! showInside;
      break;
    case MENU_OUTSIDE:
      showOutside = ! showOutside;
      break;
    case MENU_LINES:
      showLines = ! showLines;
      break;
    case MENU_XSLICING:
      if( slicing & 0x1 ) slicing -= 1;
      else slicing += 1;    
      break;
    case MENU_YSLICING:
      if( slicing & 0x2 ) slicing -= 2;
      else slicing += 2;    
      break;
    case MENU_ZSLICING:
      if( slicing & 0x4 ) slicing -= 4;
      else slicing += 4;    
      break;
    default:
      setSpace( item );
      break;
  }
  glutPostRedisplay();
}


void initMenus()
{
  glutCreateMenu(processMenu);

  glutAddMenuEntry("SQUARE",         SHAPE_SQUARE);
  glutAddMenuEntry("SPHERE",         SHAPE_SPHERE);
  glutAddMenuEntry("OVAL",           SHAPE_OVAL);
  glutAddMenuEntry("ELLIPSE",        SHAPE_ELLIPSE);
  glutAddMenuEntry("OUTSIDE_SPHERE", SHAPE_OUTSIDE_SPHERE);
  glutAddMenuEntry("STRIP",          SHAPE_STRIP);
  glutAddMenuEntry("PERIODIC",       SHAPE_PERIODIC);
  glutAddMenuEntry("ROUND_SQUARE",   SHAPE_ROUND_SQUARE);
  glutAddMenuEntry("BANANA",         SHAPE_BANANA);
  glutAddMenuEntry("BOOMERANG",      SHAPE_BOOMERANG);
  glutAddMenuEntry("TEE",            SHAPE_TEE);
#if ( DIM == 3 )
  //cylinder and cylinderZ are only implemented in 3D
  //but we could give them an implementation...to be fool-proof
  glutAddMenuEntry("CYLINDER",       SHAPE_CYLINDER);
  glutAddMenuEntry("CYLINDER_Z",     SHAPE_CYLINDER_Z);
#endif
  
  glutAddMenuEntry("-",              0);
  glutAddMenuEntry("Toggle inside  (i)",   MENU_INSIDE);
  glutAddMenuEntry("Toggle outside (o)",   MENU_OUTSIDE);
  glutAddMenuEntry("Toggle lines   (l)",   MENU_LINES);
  
  glutAddMenuEntry("Toggle x-slicing (x)", MENU_XSLICING);
  glutAddMenuEntry("Toggle y-slicing (y)", MENU_YSLICING);
  glutAddMenuEntry("Toggle z-slicing (z)", MENU_ZSLICING);
  glutAddMenuEntry("Inflation =  0 (s)",   MENU_INFLATION_0);
  glutAddMenuEntry("Reset",                MENU_RESETVIEW);
  glutAddMenuEntry("Quit",                 MENU_QUIT);
  
  glutAttachMenu(MENU_BUTTON);
}

//----------------------------------------------------------------------------
void processSpecialKey(int key, int x=0, int y=0)
{

  switch (key) {
    case GLUT_KEY_LEFT:
      slicePos -= 0.2;
      break;
    case GLUT_KEY_RIGHT:
      slicePos += 0.2;
      break;
    case GLUT_KEY_UP:
      thickness += sliceStep;
      break;
    case GLUT_KEY_DOWN:
      thickness -= sliceStep;
      if( thickness < sliceStep ) thickness = sliceStep;
      break;
    default:
      break;
  }
  
  glutPostRedisplay();
}

void processNormalKey(unsigned char c, int x=0, int y=0)
{
  
  switch (c) {
  case 27:
  case 'q':
    exit(EXIT_SUCCESS);
    
  case ' ':
    setPoints();
    break;
    
  case '0':
    resetModelView(0.04);
    break;
    
  case '+': 
    nbpts*=2; 
    if ( nbpts >= maxpts )
      nbpts = maxpts-1;
    setPoints();
    break;
    
  case '-': 
    nbpts/=2; 
    if ( nbpts < 1 )
      nbpts = 1;
    break;
  
  case 'x':
    if( slicing & 0x1 ) slicing -= 1;
    else slicing += 1;    
    break;
  
  case 'y':
    if( slicing & 0x2 ) slicing -= 2;
    else slicing += 2;
    break;
  
  case 'z':
    if( slicing & 0x4 ) slicing -= 4;
    else slicing += 4;
    break;
  
  case 'i':
    showInside = ! showInside;
    break;
  
  case 'o':
    showOutside = ! showOutside;
    break;
    
  case 'l':
    showLines = ! showLines;
    break;
    
  case 'a':
    for(int ii = 0; ii < 5; ++ii )
      boxsize_set[ii] = RNG.choose( 1, 2, 3, 5, 7, 9 );
    setSpace( spaceShape );
    break;
    
  case 's':
    inflation_mode = ( inflation_mode + 1 ) % ( sizeof(inflation)/sizeof(real) );
    inflation_radius = inflation[ inflation_mode ];
    setSpace( spaceShape );
    break;
  
  case 't':
    timer = ! timer;
    break;
    
  }
  
  glutPostRedisplay();
}

//----------------------------------------------------------------------------

bool showPoint( const int ii )
{
  return ( (!slicing) || \
           ( (slicing & 0x01) && (projections[ii].XX > (slicePos-0.5*thickness)) && (projections[ii].XX < (slicePos+0.5*thickness))) ||
           ( (slicing & 0x02) && (projections[ii].YY > (slicePos-0.5*thickness)) && (projections[ii].YY < (slicePos+0.5*thickness))) ||
           ( (slicing & 0x04) && (projections[ii].ZZ > (slicePos-0.5*thickness)) && (projections[ii].ZZ < (slicePos+0.5*thickness))) );
}


void display()
{
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

  //a yellow dot at the center
  glPointSize(10.0);
  glColor3f( 1.0, 1.0, 0 );  
  glBegin(GL_POINTS);
  glVertex3f( 0., 0., 0. );
  glEnd();
  
  //plot a gren dot for points inside, a red dot for point outside:
  glPointSize(2.0);
  glBegin(GL_POINTS);
  for( int ii = 0; ii < nbpts; ++ii ) {
    if( showPoint(ii) ) {
      if ( inside[ ii ] ) {
        if ( showInside ) {
          glColor3f( 0.0, COL, 0.0 ); glVertex( points[ii] );
          glColor3f( COL, 1.0, COL ); glVertex( projections[ii] );
        }
      } else
        if ( showOutside ) {
          glColor3f( 0.0, 0.0, COL ); glVertex( points[ii] );
          glColor3f( COL, COL, 1.0 ); glVertex( projections[ii] );
        }
    }
  }
  glEnd();

  if ( showLines ) {
    //plot a blue line from the point to its projection:
    glLineWidth(line_width);
    glBegin(GL_LINES);
    for( int ii = 0; ii < nbpts; ++ii ) {
      if( showPoint(ii) ) {
        if ( inside[ ii ] ) {
          if ( showInside ) {
            glColor3f( 0.0, COL, 0.0 ); glVertex( points[ii] );
            glColor3f( COL, 1.0, COL ); glVertex( projections[ii] );
          }
        }
        else {
          if ( showOutside ) {
            glColor3f( 0.0, 0.0, COL ); glVertex( points[ii] );
            glColor3f( COL, COL, 1.0 ); glVertex( projections[ii] );
          }
        }
      }
    }
    glEnd();
  }
  
  if ( showProject ) {
    glLineWidth(line_width);
    glBegin(GL_LINES);
    for( int ii = 0; ii < nbpts; ++ii ) {
      if( showPoint(ii) ) {
        glColor3f( COL, 0.0, 0.0 );
        glVertex( projections[ii] );
        glColor3f( 1.0, COL, COL );
        glVertex( projections2[ii] );
      }
      glEnd();
    }
  }
  
  char text[128];
  snprintf(text, sizeof(text), "inflation=%5.1f, error=%5.2e (press s)", 
           inflation_radius, max_error_projection);
  displayText(text, 1);
  
  glFlush();
}

//----------------------------------------------------------------------------
void initGLUT()
{
  initGLUT_default();
  initMenus();
  
  if ( DIM == 3 )
    glEnable(GL_DEPTH_TEST);

//  glDisable(GL_BLEND);
  
  setSpace();
  setPoints();
  
  glutTimerFunc(delay, timerFunction, timer);
  
  resetModelView(0.04);
}


//----------------------------------------------------------------------------
int main(int argc, char *argv[])
{
  glutInit(&argc, argv);

  glutInitDisplayMode( GLUT_RGBA | GLUT_DEPTH );
  glutInitWindowSize(600, 600); 
  glutInitWindowPosition(50, 50);
  glutCreateWindow (argv[0]);

  initGLUT();
  
  glutDisplayFunc(display);
  glutReshapeFunc(reshapeWindow);
  glutMouseFunc(processMouse);
  glutMotionFunc(processMotion);
  glutSpecialFunc(processSpecialKey);
  glutKeyboardFunc(processNormalKey);
  
  glutMainLoop();

  return EXIT_SUCCESS;
}


