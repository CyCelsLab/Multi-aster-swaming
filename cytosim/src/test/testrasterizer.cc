//RCS: $Id: testrasterizer.cc,v 2.1 2005/01/10 17:27:47 foethke Exp $
// Visual test for the rasterizer used in attachment algorithm of Cytosim
// Francois Nedelec, nedelec@embl.de, October 2002

#define DEFINE_MENUS
//#define USE_OPENGL_FOG
#define DIM 3

#include <ctime>
#include "glut_base.cc"
#include "types.h"
#include "smath.h"
#include "random.h"
#include "rasterizer.h"

//------------- delay for the Timer function:

int delay = 0;


const int size = 20;
const int maxpts = 10;

int nbpts = 2;
real pts[ 3 * maxpts ];
real width = 5;


void manyTests();

void processNormalKey(unsigned char c, int x=0, int y=0)
{
  switch (c) {
  case 27:
  case 'q':
    exit(EXIT_SUCCESS);
  case ' ':
    for( int i = 0; i < 3*maxpts ; ++i )
      pts[i] = (size-1) * RNG.sreal();
    break;
  case '0':
    resetModelView();
    break;
  case '+': ++nbpts; break;
  case '-': --nbpts; break;
  case 'r':
      manyTests();
      break;
  default:
    printf("short manual:\n\
            space : draw new distribution\n\
            +     : increase number of points\n\
            -     : decrease number of points\n\
            r     : many tests as fast as possible\n");
  }
  glutPostRedisplay();
}

//===================================================================
//===================================================================

real min[]   = {0, 0, 0};
real delta[] = {1, 1, 1};

void dummyPaint(int x, int y, int z) 
{
}

void glPaint(int x, int y, int z)
{
  glPointSize(1.0);
  glBegin(GL_POINTS);
  glColor3f(1.0, 1.0, 1.0);

#if ( DIM == 2 )
  glVertex2i( x, y );
#else
  glVertex3i( x, y, z );
#endif

  glEnd();
}


void display2D()
{
  int i,j;
  //--------------draw a grid in gray:
  glPointSize(1.0);
  glBegin(GL_POINTS);
  glColor3f(0.3, 0.3, 0.3);
  for( i = -size-5; i <= size+5; i += 1)
  for( j = -size-5; j <= size+5; j += 1)
    glVertex2i(i, j);
  glEnd();


  glPointSize(6.0);
  glBegin(GL_POINTS);
  glColor3f(0, 0, 1.0);
  for( i = 0; i < nbpts ; ++i )
    glVertex2f( pts[2*i], pts[2*i+1] );
  glEnd();


  if ( nbpts == 2 ) {
      Rasterizer::paintFatLine2D( glPaint, pts, pts+DIM, width, min, delta );
    return;
  }

  int nb = nbpts;
  Rasterizer::convexHull2D( pts, &nb );
   Rasterizer::paintPolygon2D(glPaint, nb, pts);  

  glBegin(GL_LINE_LOOP);
  glColor3f(0, 1.0, 0);
  for( i = 0; i < nb ; ++i )
    glVertex2f( pts[2*i], pts[2*i+1] );
  glEnd();
}


void display3D()
{
  /*
   //--------------draw a grid:
   glPointSize(1.0);
   glBegin(GL_POINTS);
   glColor3f(1.0, 0, 0);
   for( i = -size; i <= size; i += 1)
   for( j = -size; j <= size; j += 1)
   //for( k = -size; k <= size; k += 1)
   glVertex3i(i, j, k);
   glEnd();
   */


  glPointSize(6.0);
  glBegin(GL_POINTS);
  glColor3f(0, 0, 1.0);
  for(int i = 0; i < nbpts ; ++i )
    glVertex3f( pts[3*i], pts[3*i+1], pts[3*i+2] );
  glEnd();
  
  if ( nbpts == 2 )
      Rasterizer::paintFatLine3D(glPaint, pts, pts+DIM, width, min, delta);
  else
      Rasterizer::paintPolygon3D(glPaint, nbpts, pts);
}



void display()
{
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

  //------------draw the origin:
  glPointSize(2.0);
  glBegin(GL_POINTS);
  glColor3f(1.0, 0.0, 1.0);
  glVertex2i(0, 0);
  glEnd();

  //a vecteur for the Z direction:
  glBegin(GL_LINES);
  glColor3f(1.0, 0.0, 1.0);
  glVertex3i(0, 0, 0);
  glVertex3f(0, 0, 1.5);
  glEnd();

#if (DIM == 2)
  display2D();
#else
  display3D();
#endif

  glFlush();
  glutSwapBuffers();
}

void manyTests()
{
    //perform 100 visual tests:
    for( int test = 0; test < 100; ++test ) {
        for( int i = 0; i < 3*maxpts ; ++i )
            pts[i] = (size-1) * RNG.sreal();
        display();
    }
    
    int NB = 10000;
    clock_t start = clock();
    //performs NB tests with the dummy paint function:
    for( int test = 0; test < NB; ++test ) {
        for( int i = 0; i < 3*maxpts ; ++i )
            pts[i] = (size-1) * RNG.sreal();

        if ( nbpts == 2 )
            Rasterizer::paintFatLine3D(dummyPaint, pts, pts+DIM, width, min, delta);
        else
            Rasterizer::paintPolygon3D(dummyPaint, nbpts, pts);
        
    }
    printf("%i tests in %lu clock ticks\n", NB, clock()-start);
}


void initGLUT()
{
  initGLUT_default();
  
  glEnable(GL_DEPTH_TEST);

  resetModelView(0.04);
  processNormalKey(' ');
}


int main(int argc, char *argv[])
{
  glutInit(&argc, argv);

  glutInitDisplayMode( GLUT_DOUBLE | GLUT_RGBA ); //| GLUT_DEPTH );
  glutInitWindowSize(600, 600); 
  glutInitWindowPosition(50, 50);
  glutCreateWindow (argv[0]);

  initGLUT();

  glutDisplayFunc(display);
  glutReshapeFunc(reshapeWindow);
  glutMouseFunc(processMouse);
  glutMotionFunc(processMotion);
  //glutSpecialFunc(processInputKey);
  glutKeyboardFunc(processNormalKey);

  glutMainLoop();

  return EXIT_SUCCESS;
}


