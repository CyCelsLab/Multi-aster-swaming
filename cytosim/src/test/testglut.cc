//RCS: $Id: testglut.cc,v 2.2 2005/04/10 15:03:58 nedelec Exp $
//--------------------------------------------------------------------------
//                               testglut.cc
//      mouse driven zoom and rotation with quaternions and GLUT unproject
//
//   Francois Nedelec nedelec@embl.de  October 2002, modified Nov. 2003

#define DEFINE_MENUS
#define USE_OPENGL_FOG
#include "glut_base.cc"

#include "pointsonsphere.h"

int nbpts = 30;
PointsOnSphere S(30);


//----------------------------------------------------------------------------
void processNormalKey(unsigned char c, int x, int y)
{
  int force = 0;
  switch (c) {
  case 27:
  case 'q':
    exit(EXIT_SUCCESS);
  case '0':
    resetModelView();
    break;  

  case 'u': nbpts-=10; break;
  case 'i': nbpts-=1; break;
  case 'o': nbpts+=1; break;
  case 'p': nbpts+=10; break;
  case '[': nbpts+=100; break;
  case ']': nbpts+=1000; break;
  case ' ': force = 1; break;

  default:
    printf("normal key %c %i %i\n", c, x, y);
  }

  if ( force || ( nbpts != S.nbPoints() )) {
    S.distributePoints( nbpts );
    S.reportConvergence();
  }

  glutPostRedisplay();
}

//----------------------------------------------------------------------------
void display()
{
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

  glColor3f(1.0, 0.0, 0.0);
  glLineWidth(2.0);
  glutWireCube(2.0);

  glPointSize(2.0);
  glBegin(GL_POINTS);
  glColor3f(1.0, 1.0, 1.0);
  glVertex3f(0.0, 0.0, 0.0);
  glColor3f(1.0, 0.0, 0.0);
  glVertex3f(0.5, 0.0, 0.0);
  glColor3f(0.0, 1.0, 0.0);
  glVertex3f(0.0, 0.5, 0.0);
  glColor3f(0.0, 0.0, 1.0);
  glVertex3f(0.0, 0.0, 0.5);
  glEnd();

  glPointSize(8.0);
  glBegin(GL_POINTS);
  glColor4f(0.0, 0.0, 1.0, 0.8);
  for( int ii=0; ii < S.nbPoints(); ++ii )
#ifdef SINGLE
    glVertex3fv( S.addressForThisPoint(ii) );
#else
    glVertex3dv( S.addressForThisPoint(ii) );
#endif
  glEnd();

  char displayed_text[128];
  snprintf(displayed_text, sizeof(displayed_text), "%i points", S.nbPoints());
  displayText(displayed_text );
    
  glFinish();
  glutSwapBuffers();
  glPrintErrors("display");
}

//----------------------------------------------------------------------------
void initGLUT()
{
  initGLUT_default();
  glEnable(GL_DEPTH_TEST);
  resetModelView(0.5);
}

//----------------------------------------------------------------------------
void glStat()
{
  printf("has keyboard %i\n", glutDeviceGet(GLUT_HAS_KEYBOARD));
  printf("has mouse %i\n", glutDeviceGet(GLUT_HAS_MOUSE));
  printf("nb of mouse buttons %i\n", glutDeviceGet(GLUT_NUM_MOUSE_BUTTONS));
  printf("color bit depth %i\n", glutGet(GLUT_WINDOW_BUFFER_SIZE));
  printf("alpha bit depth %i\n", glutGet(GLUT_WINDOW_ALPHA_SIZE));
  printf("Current display is RGBA %i\n", glutGet(GLUT_WINDOW_RGBA));
  printf("Current display mode possible %i\n", glutGet(GLUT_DISPLAY_MODE_POSSIBLE));

  //anti-aliasing of points and lines:
  printf("GL_POINT_SMOOTH enabled %i\n", glIsEnabled(GL_POINT_SMOOTH));
  GLfloat s[2];
  glGetFloatv(GL_SMOOTH_POINT_SIZE_RANGE, s);
  printf("GL_SMOOTH_POINT_SIZE_RANGE %f %f\n", s[0],s[1]);

  printf("GL_LINE_SMOOTH enabled %i\n", glIsEnabled(GL_LINE_SMOOTH));
  glGetFloatv(GL_SMOOTH_LINE_WIDTH_RANGE, s);
  printf("GL_SMOOTH_LINE_WIDTH_RANGE %f %f\n", s[0], s[1]);
  
  GLint numBuffers;
  glGetIntegerv(GL_AUX_BUFFERS, &numBuffers);
  printf("GL_AUX_BUFFERS %i\n", int(numBuffers));
}

//----------------------------------------------------------------------------
int main(int argc, char *argv[])
{
  S.reportConvergence();
  glutInit(&argc, argv);
	
  glutInitDisplayMode( GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH );

  glutInitWindowSize(400, 400); 
  glutInitWindowPosition(50, 50);
  glutCreateWindow(argv[0]);

  initGLUT();

  glutDisplayFunc(display);
  glutReshapeFunc(reshapeWindow);
  glutMouseFunc(processMouse);
  glutMotionFunc(processMotion);
  glutKeyboardFunc(processNormalKey);

  //glStat();
  glutMainLoop();

  return EXIT_SUCCESS;
}
