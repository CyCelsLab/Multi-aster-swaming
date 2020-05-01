//RCS: $Id: testsolve.cc,v 2.2 2005/01/10 17:27:47 foethke Exp $
//--------------------------------------------------------------------------
//
// nedelec@embl-heidelberg.de  April 2003

#include <cstdlib>
#include <cstdio>
#include "glut.h"

#include "smath.h"
#include "vecteur2.h"
#include "matrix2.h"
#include "random.h"

//a point in space:
GLdouble Gx=1, Gy=1, Gz;

Vecteur2 pg(5,0,0), px(1,0,0), pf, pn, pc, pp;

//rate of calculation:
int delay = 50;
int normalize = 1;
int correct = 1;
int mobile = 1;

real radius = 1;
real dt = 1./16;
real km = 1;
real noise = 1;

void timerFunction(int value)
{
  real s, h = km * dt;
  Matrix22 I, P, D, C;
  I.setToOne();
  real xn = px.norm();
  Vecteur2 pxn = px.normalized();
  
  P(0,0) = 1.0 - pxn[0] * pxn[0];
  P(1,0) =     - pxn[1] * pxn[0];
  P(0,1) =     - pxn[0] * pxn[1];
  P(1,1) = 1.0 - pxn[1] * pxn[1];

  Vecteur2 rhs, f, fx;
  Vecteur2 random_force( RNG.gauss()*noise/km, RNG.gauss()*noise/km, 0);
  

  printf("\nkm*dt %2f noise %f correct %i norm %i :\n",
	 km*dt, noise, correct, normalize);
  printf("px    : "); px.println();
  //P.println();

  //free point obtained without the constraints:
  pf = ( px + h * ( pg + random_force ) ) / ( 1.0 + h );
  printf("free  : "); pf.println();

  //projection of the free point on the constraints
  pp = radius * pf.normalized();
  printf("proj  : "); pp.println();

  //adding the projector in the dynamic matrix
  rhs = px + h * ( P * ( pg + random_force ));
  pc = (I + h * P).inverse() * rhs;
  printf("const : "); pc.println();

  f = ( pg - px );
  s = ( pxn * f ) / xn;
  fx = 2.0 * s * pxn - f / xn;

  //projector with its corrections due to the derivative of the constraints:
  C(0,0) = pxn[0] * fx[0] - s;
  C(1,0) = pxn[1] * fx[0];
  C(0,1) = pxn[0] * fx[1];
  C(1,1) = pxn[1] * fx[1] - s;

  //C.println();
  rhs = px + h * ( P * ( pg + random_force - C * px ));
  D = I + h * ( P * ( I - C ));
  //D.println();
  //D.inverseit();

  pn = D.inverse() * rhs;
  printf("new   : ");   pn.println();

  if ( mobile )    px = correct ? pn : pc;
  if ( normalize ) px = radius * px.normalized();
  
  glutPostRedisplay();
  glutTimerFunc( delay, timerFunction, 1);
}


void display()
{
  glClear(GL_COLOR_BUFFER_BIT);

  glColor3f(0,0,1.0);
  glBegin(GL_LINE_LOOP);
  for( real ii = 0 ; ii < 6.28; ii+=0.0314 )
    glVertex2f( radius*cos(ii), radius*sin(ii) );
  glEnd();

  glPointSize(7.0);
  glBegin(GL_POINTS);

  glColor3f(1.0, 1.0, 1.0);
  glVertex2f(px[0], px[1]);

  glColor3f(0.0, 0.0, 1.0);
  glVertex2f(pg[0], pg[1]);

  glColor3f(1.0, 0.0, 0.0);
  glVertex2f(pf[0], pf[1]);

  glColor3f(0.0, 1.0, 0.0);
  glVertex2f(pc[0], pc[1]);
  
  glColor3f(0.0, 0.0, 1.0);
  glVertex2f(pn[0], pn[1]);

  glColor3f(1.0, 0.0, 1.0);
  glVertex2f(pp[0], pp[1]);

  glEnd();

  glColor3f(1.0,1.0,1.0);
  glBegin(GL_LINE_LOOP);
  glVertex2f(px[0], px[1]);
  glVertex2f(pg[0], pg[1]);
  glEnd();

  
  glFlush();
}


//----------------------------------------------------------------------------

//size of viewing box:
GLdouble vsize[] = { 10, 10 };

//the zoom factor:
real zoom_save, zoom = 1;
int mouse_action, mouse_x, mouse_y;
GLint viewport[4];


//----------------matrices to compute the inverse projection of mouse locations

GLdouble mat_model[16];
GLdouble mat_proj[16];

#define MOUSE_ZOOM    GLUT_RIGHT_BUTTON
#define MOUSE_SET     GLUT_LEFT_BUTTON
#define MENU_BUTTON   GLUT_MIDDLE_BUTTON

void setModelView()
{
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  glScaled( zoom, zoom, zoom );
  glutPostRedisplay();
}


void processNormalKey(unsigned char c, int x, int y)
{
  switch (c) {
  case 27:
  case 'q':
    exit(EXIT_SUCCESS);

  case 'o':
    delay *= 2; break;
  case 'p':
    if ( delay > 1 ) delay /= 2; break;

  case 'K':
    km *= 2; break;
  case 'J':
    km /= 2; break;
  case 'k':
    km += 1; break;
  case 'j':
    km -= 1; break;
    
  case 'i':
    noise *= 2; break;
  case 'u':
    noise /= 2; break;
    
  case 'n':
    normalize = !normalize;
    break;    

  case 'c':
    correct = !correct;
    break;
    
  case 'm':
    mobile = !mobile;
    break;
    
  case 'z':
    px.set( 1, 0, 0);
    pg.set( 5, 0, 0 );
    break;
    
  default:
    printf("normal key %c %i %i\n", c, x, y);
  }
}

enum MENUS_ID { MENU_QUIT };

void processMenu(int item)
{
  if ( item == MENU_QUIT ) exit( EXIT_SUCCESS );
}

void reshapeWindow(int w, int h)
{
  glViewport(0, 0, w, h);
  glGetIntegerv(GL_VIEWPORT, viewport);

  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  double ratio = w * vsize[1] / double( vsize[0] * h );
  
  if ( ratio > 1 )
    glOrtho(-vsize[0], vsize[0], -vsize[1]/ratio, vsize[1]/ratio, 0, 1);
  else
    glOrtho(-vsize[0]*ratio, vsize[0]*ratio, -vsize[1], vsize[1], 0, 1);
}

//----------------------------------------------------------------------------

void processMouse(int button, int state, int x, int y)
{
//  printf("button %i %i %i %i\n", button, state, x, y);
  if ( state != GLUT_DOWN ) return;
  mouse_action = button;

  mouse_x = x;
  mouse_y = y;

  switch( mouse_action ) {

  case MOUSE_ZOOM:
    zoom_save = zoom;
    break;

  case MOUSE_SET:
    glGetDoublev(GL_MODELVIEW_MATRIX, mat_model);
    glGetDoublev(GL_PROJECTION_MATRIX, mat_proj);
    gluUnProject(x, viewport[3]-y, 0, mat_model, mat_proj, viewport,
		 &Gx, &Gy, &Gz);
    pg.set(Gx, Gy, Gz);
    glutPostRedisplay();
    break;
  }
}


void processMotion(int x, int y)
{
  real d;
  switch( mouse_action ) {

  case MOUSE_ZOOM:
    
    d = 1.0 + 4 * real( x - mouse_x ) / viewport[2];
    if ( d > 0 ) zoom = zoom_save * d;
    setModelView();
    break;
    
  case MOUSE_SET:
    glGetDoublev(GL_MODELVIEW_MATRIX, mat_model);
    glGetDoublev(GL_PROJECTION_MATRIX, mat_proj);
    gluUnProject(x, viewport[3]-y, 0, mat_model, mat_proj, viewport,
		 &Gx, &Gy, &Gz);
    pg.set(Gx, Gy, Gz);
    glutPostRedisplay();
    break;
  }
}


void initGLUT()
{
  glClearColor(0.0, 0.0, 0.0, 0.0);
  glEnable(GL_POINT_SMOOTH);
  glEnable(GL_LINE_SMOOTH);
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glHint(GL_POINT_SMOOTH_HINT, GL_FASTEST);
  glHint(GL_LINE_SMOOTH_HINT, GL_FASTEST);

  glutCreateMenu(processMenu);
  glutAddMenuEntry("Quit", MENU_QUIT);
  glutAttachMenu(MENU_BUTTON);

  setModelView(); 
  glutTimerFunc( delay, timerFunction, 1);
}

int main(int argc, char *argv[])
{
  glutInit(&argc, argv);
	
  glutInitDisplayMode( GLUT_SINGLE | GLUT_RGBA );
  glutInitWindowSize(400, 400); 
  glutInitWindowPosition(50, 50);
  glutCreateWindow (argv[0]);

  initGLUT();

  glutDisplayFunc(display);
  glutReshapeFunc(reshapeWindow);
  glutMouseFunc(processMouse);
  glutMotionFunc(processMotion);
  glutKeyboardFunc(processNormalKey);

  glutMainLoop();

  return EXIT_SUCCESS;
}
