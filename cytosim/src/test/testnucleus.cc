// RCS: $Id: testnucleus.cc,v 2.5 2005/04/10 15:06:52 nedelec Exp $
//--------------------------------------------------------------------------
//
// foethke@embl-heidelberg.de July 2003

//#include <iostream>
#include <cstdlib>
#include <cstdio>
#include "smath.h"
//#include <GL/gl.h>
//#include <GL/glext.h>
//#include <GL/glut.h>
#include "glut.h"

#include "types.h"
#include "main.h"
#include "random.h"
#include "cblas.h"
#include "clapack.h"
#include "vecteur.h"
#include "smath.h"
#include "gmres.h"
#include "conjgradient.h"
#include "nucleus.h"

//#define PI 3.1415926535898

#ifndef CALLBACK
#define CALLBACK
#endif

//using namespace std;

//-------------------------variables for the nucleus----------------------------

Nucleus * nuc;

//rate of calculation:
int  delay     = 50;
int  normalize = 1;
bool correct   = 1;
bool ppprime   = 1;
bool explCalc  = 0;

//real noise     = 1.;
real thermalNoise = 0;

real *pg            = 0;
real *noiseVec      = 0;
real *axis          = 0;
real *E             = 0;
real *dummy1        = 0;
real *dummy2        = 0;

real *force         = 0;

real *b             = 0;
real *j             = 0;
real *jjt           = 0;
real *jjtInvj       = 0;
real *p             = 0;

real *dx            = 0;
real *work          = 0;

//-----------------------------------------------------------------------------
// These are new variables and functions that are needed for compatibility with
// pointset.h of cytosim

int  allocatedP     = 0;


//--------------------------variables for OpenGL--------------------------------

const int screenWidth    = 800;
const int screenHeight   = 800;
GLdouble ratio             = (GLdouble)screenWidth/(GLdouble)screenHeight;
const GLdouble glWindowHeight = 5.0;
const GLdouble glWindowWidth  = glWindowHeight*ratio;

GLdouble cameraPos[3];                        // position of the camera
GLdouble lookAt[3];                           // where to look at
GLdouble thisSideUp[3];                       // this is up
GLdouble translation[3] = {0.0, 0.0, 0.0};    // total translation of the center
GLdouble transSave[3]   = {0.0, 0.0, 0.0};    // backup of the translation
GLdouble rotMat[16];                          // total rotation around the center
GLdouble rotMatSave[16];                      // backup of the rotation matrix
GLdouble zoom = 1.0;                          // zoom-factor
GLdouble zoomSave;                            // backup of the zoom-factor
GLdouble mouseX;                              // x-coordinate of the mouseclick
GLdouble mvMatrix[16], projMatrix[16];
GLint    viewPort[4];
GLdouble clickPointVec[3];
GLdouble nclickPointVec[3];
GLuint*  selBuf;                              // the selection buffer
GLint    selection;                           // the name of the selection
GLint    hits;                                // the number of selection hits
GLUquadricObj* qobj;                          // a pointer to quadric objects
GLdouble sphereRotMat[16];                    // rotation Matrix for the sphere
GLdouble opaque = 0.6;                        // opaqueness of the nucleus

//----------------list of different actions control by the mouse----------------

enum MOUSE_ACTION { MOUSE_NOTHING, MOUSE_TRANSLATE, MOUSE_ZOOM,
                    MOUSE_ROTATE, MOUSE_SPIN, MOUSE_GRAB };

MOUSE_ACTION mouse_action;

MOUSE_ACTION mouseAction( int button, int modifier ) {

  switch( button ) {
  case GLUT_LEFT_BUTTON:
    switch( modifier ) {
    case GLUT_ACTIVE_CTRL:
      return MOUSE_GRAB;
    case GLUT_ACTIVE_SHIFT:
      return MOUSE_TRANSLATE;
    case GLUT_ACTIVE_ALT:
      return MOUSE_NOTHING;
    default:
      return MOUSE_ROTATE;
    }
  case GLUT_MIDDLE_BUTTON:
    return MOUSE_ZOOM;
    //case GLUT_RIGHT_BUTTON:    
  }
  return MOUSE_NOTHING;
}

//------------functions for the calculation of the rotation matrix--------------

GLdouble norm(const GLdouble *vector)
{
  return(sqrt(vector[0]*vector[0] + vector[1]*vector[1] + vector[2]*vector[2]));
}

void createRotMat(const GLdouble *v, GLdouble *mat)
{
  // creates a rotation around the axis v

  GLdouble angle  = norm(v);
  GLdouble u[3]   = {v[0]/angle, v[1]/angle, v[2]/angle};
  GLdouble sinFac = sin(0.5*angle);
  GLdouble cosFac = 1.0 - cos(0.5*angle);

  mat[0]          = cosFac*u[0]*u[0] + cos(0.5*angle);
  mat[1]          = cosFac*u[1]*u[0] + sinFac*u[2];
  mat[2]          = cosFac*u[2]*u[0] - sinFac*u[1];
  mat[3]          = 0.0;
  mat[4]          = cosFac*u[0]*u[1] - sinFac*u[2];
  mat[5]          = cosFac*u[1]*u[1] + cos(0.5*angle);
  mat[6]          = cosFac*u[2]*u[1] + sinFac*u[0];
  mat[7]          = 0.0;
  mat[8]          = cosFac*u[0]*u[2] + sinFac*u[1];
  mat[9]          = cosFac*u[1]*u[2] - sinFac*u[0];
  mat[10]         = cosFac*u[2]*u[2]+cos(0.5*angle);
  mat[11]         = 0.0;
  mat[12]         = 0.0;
  mat[13]         = 0.0;
  mat[14]         = 0.0;
  mat[15]         = 1.0;
}

void initRotMat(GLdouble *mat)
{
  for(int i = 0; i < 16; i++) mat[i] = 0.0;
  for(int i = 0; i < 16; i += 5) mat[i] = 1.0;
}

void copyMat(const GLdouble *source, GLdouble *dest)
{
  // copy matrix source to matrix dest
  
  for(int i = 0; i < 16; i++) dest[i] = source[i];
}

void multMat(const GLdouble *leftMat, GLdouble *rightMat)
{
  // calculates leftMat.rightMat, the result is written in rightMat
 
  GLdouble result[16];
  for(int i = 0; i < 16; i++) result[i] = 0.0;  

  for(int i = 0; i < 4; i++) {
    for(int k = 0; k < 4; k++) {
      for(int j = 0; j < 4; j++) {
	result[i+k*4] += leftMat[i+j*4]*rightMat[j+k*4];
      }
    }
  }
  
  copyMat(result, rightMat);
}

//------------------------------------------------------------------------------

void initialise() {

  // build dimensions-array
  real dimensions[3] = {glWindowWidth, glWindowHeight, glWindowWidth};

  // create nucleus
  Nucleus::space = newSpace(1, dimensions);

  MP.parseFile("test.in");
  nuc = new Nucleus(0);
  
  // theses arrays are needed for all calculation methods
  // no handle for the reference points
  pg = new real[DIM*(nuc->nbPoints()-nuc->nbrefpts)];
  
  // initialise points
  if(DIM == 2) {
    
    pg[0] = 0.;
    pg[1] = 0.;
    //pg[2] = nuc->whereP(0)[0] + nuc->radius() + 3.;
    //pg[3] = nuc->whereP(0)[1];
    //if(nuc->nbSegments() > 1) {
    //  pg[4] = nuc->whereP(0)[0] - nuc->radius() - 10.;
    //  pg[5] = nuc->whereP(0)[1];
   // }
    
    // no handle for the reference points
    for(int i = 1; i < nuc->nbPoints()-nuc->nbrefpts; i++) {
      //randAnglePhi     = 2.*PI*RNG.preal();
      pg[i*DIM]        = nuc->whereP(i)[0]*(1. + 3./nuc->radius());
      pg[i*DIM + 1]    = nuc->whereP(i)[1]*(1. + 3./nuc->radius());
    }
  }
  else if(DIM == 3) {    
    pg[0] = 0.;
    pg[1] = 0.;
    pg[2] = 0.;
    //pg[3] = nuc->whereP(0)[0] + nuc->radius() - 3.;
    //pg[4] = nuc->whereP(0)[1];
    //pg[5] = nuc->whereP(0)[2];
    
    // no handle for the reference points
    for(int i = 1; i < nuc->nbPoints()-nuc->nbrefpts; i++) {
      pg[i*DIM]        = nuc->whereP(i)[0]*(1. + 3./nuc->radius());
      pg[i*DIM + 1]    = nuc->whereP(i)[1]*(1. + 3./nuc->radius());
      pg[i*DIM + 2]    = nuc->whereP(i)[2]*(1. + 3./nuc->radius());
    }
  }
  printf("reached end of initialise!\n");

}

void cleanUp()
{
  printf("Cleaning up...\n");
  
  // theses arrays are needed for all calculation methods
  delete[] pg;
  delete[] noiseVec;
  delete[] force;
  delete[] dummy1;
  delete[] dummy2;

  // these arrays are only needed for implicit
  delete[] b;
  delete[] j;
  delete[] jjt;
  delete[] jjtInvj;
  delete[] p;
  
  delete[] dx;
}

/*****************************************************************************/

void allocate()
{
  int nbMTOCs = nuc->nbSegments();
  int vecDim  = DIM*nuc->nbPoints();

  if (allocatedP < vecDim) {
    printf("Allocating more memory in testnucleus!\n");
    if( noiseVec )      delete[] noiseVec;
    if( force )         delete[] force;
    if( dummy1 )        delete[] dummy1;
    if( dummy2 )        delete[] dummy2;

    if( b )             delete[] b;
    if( j )             delete[] j;
    if( jjt )           delete[] jjt;
    if( jjtInvj )       delete[] jjtInvj;
    if( p )             delete[] p;

    if( dx )            delete[] dx;
    if( work )          delete[] work;

    allocatedP    = vecDim + 2*DIM;
    // theses arrays are needed for all calculation methods
    noiseVec      = new real[allocatedP];
    force         = new real[allocatedP];
    dummy1        = new real[allocatedP];
    dummy2        = new real[allocatedP];
        
    // these arrays are only needed for implicit
    b             = new real[allocatedP];
    j             = new real[nbMTOCs*allocatedP];
    jjt           = new real[nbMTOCs*nbMTOCs];
    jjtInvj       = new real[nbMTOCs*allocatedP];
    p             = new real[allocatedP*allocatedP];
        
    // these are working-arrays for the solvers
    dx            = new real[allocatedP];
    work          = new real[allocatedP];
  }
}


void calculateExtForce(real* extForce)
{
  // no forces on the center point
  for(int i = 0; i < DIM; i++)
    extForce[i] = 0.0;
  
  for(int i = DIM; i < DIM*(nuc->nbPoints()-nuc->nbrefpts); i++)
    extForce[i] = MP.km*(pg[i] - nuc->coord(i));

  // no forces on the reference points
  for(int i = DIM*(nuc->nbPoints()-nuc->nbrefpts); i < DIM*nuc->nbPoints(); i++)
    extForce[i] = 0.0;
}

/*****************************************************************************/

// these functions are for the iterative methods

void calculateMatVec(const real *x, real *Ax)
{
  int vecDim = DIM*nuc->nbPoints();
    
  // calculate Y = A*X
  // no forces on the center and the reference points
  for(int i = 0; i < DIM; i++)
    dummy1[i] = 0.;
  for(int i = DIM; i < vecDim-nuc->nbrefpts*DIM; i++)
    dummy1[i] = -MP.km*x[i];
  for(int i = vecDim-nuc->nbrefpts*DIM; i < vecDim; i++)
    dummy1[i] = 0.;
  
  // add the P'-term if we calculate with correction
  if( correct )
    nuc->setSpeedsWithCorrection(x, dummy1, Ax, -MP.dt);
    //nuc->addProjectionDiff(x, dummy1);
  else {
    // project and scale with dt
    nuc->setSpeedsFromForces(dummy1, Ax, -MP.dt);
    //nuc->setProjectedForces(dummy1, Ax);
    //blas_xscal(vecDim,-MP.dt*nuc->nurotmobil,Ax,1);
  }
  
  // add unit*x
  blas_xaxpy(vecDim,1.,x,1,Ax,1);  
}

/*****************************************************************************/

void timerFunction(int value)
{
  //  printf("reached timerFunction!\n");
  
//  int nbMTOCs = nuc->nbSegments();
  int vecDim  = DIM*nuc->nbPoints();

  allocate();
  nuc->prepareProjection();
  
  //thermalNoise = MP.km*MP.dt*nuc->setBrownian(noiseVec);
  //thermalNoise = MP.km*MP.dt*nuc->setBrownian(dummy1);

  if( explCalc ) {
    calculateExtForce(force);
    nuc->setSpeedsFromForces(force, dummy1, MP.dt);
    //nuc->setProjectedForces(force, dummy1);
    //blas_xscal(vecDim,MP.dt*nuc->nurotmobil,dummy1,1);
            
    // add noise
    nuc->addBrownianMotion(dummy1);
    
    blas_xaxpy(vecDim,1.0,dummy1,1,nuc->pts(),1);
  }
  else {
    int nIter = vecDim;
    real tolerance = 1e-2;
    
//    for(int d = 0; d < 10000; d++) {
//      printf("slow down!\n");
//    }
    
    // calculate the right hand side of Ax = b
    calculateExtForce(force);
    nuc->setSpeedsFromForces(force, b, MP.dt);
    //nuc->setProjectedForces(force, b);
    //blas_xscal(vecDim,MP.dt*nuc->nurotmobil,b,1);
    
    if( correct ) nuc->prepareProjectionDiff(force);
    
    // add noise
    thermalNoise = nuc->addBrownianMotion(b);
    //printf("thermalNoise: %f\n", thermalNoise);
                
    // call the solver
    for(int ii=0; ii < vecDim; ii++)
      dx[ii] = 0.;
    //Solver::GMRESHH(vecDim, b, dx, calculateMatVec, nIter, tolerance, thermalNoise);
    Solver::biCGstab(vecDim, b, dx, calculateMatVec, nIter, tolerance, thermalNoise);
    //printf("converged after %d steps of %d, tolerance: %e\n", nIter, vecDim, tolerance);
    //solverGMRES(vecDim, b, dx, calculateMatVec, nIter);
    //solverFOM(vecDim, b, dx, calculateMatVec, nIter);
    
    for(int i = 0; i < vecDim; i++) {
      nuc->pts(i) += dx[i];
    }
  }
    
  // get rid of finite-step errors but conserve the shape
  // and the center of gravity for the nucleus
  nuc->reshape();
  //nuc->reshapeWithConservedCenter();
    
  glutPostRedisplay();
  glutTimerFunc( delay, timerFunction, 1);
}

//-----------------------------OpenGL functions---------------------------------

void initGlut2D()
{
  int nbMTOCs = nuc->nbSegments();
  
  // set the viewpoint
  cameraPos[0]   = 0.0;   // x-position of the camera
  cameraPos[1]   = 0.0;   // y-position of the camera
  cameraPos[2]   = 10.0;  // z-position of the camera
  lookAt[0]      = 0.0;   // where to look at
  lookAt[1]      = 0.0;   // where to look at
  lookAt[2]      = 0.0;   // where to look at
  thisSideUp[0]  = 0.0;   // this is up
  thisSideUp[1]  = 1.0;   // this is up
  thisSideUp[2]  = 0.0;   // this is up
  
  glLoadIdentity();
  gluLookAt(cameraPos[0], cameraPos[1], cameraPos[2],
	    lookAt[0], lookAt[1], lookAt[2],
	    thisSideUp[0], thisSideUp[1], thisSideUp[2]);

  // set lighting parameters

  GLfloat matAmbDiff[]    = { 0.0, 0.5, 1.0, 1.0 };
  GLfloat matSpecular[]   = { 1.0, 1.0, 1.0, 1.0 };
  GLfloat matShininess[]  = { 50.0 };
  GLfloat light0Pos[]     = { 2.0, -3.0,  2.0, 0.0 };
  GLfloat light1Pos[]     = {-2.0, -3.0, -2.0, 0.0 };
  GLfloat light0Diffuse[] = { 1.0, 1.0, 1.0, 1.0 };
  GLfloat light0Spec[]    = { 1.0, 1.0, 1.0, 1.0 };
  GLfloat lModelAmbient[] = { 0.2, 0.2, 0.2, 1.0 };
  
  glClearColor(0.0, 0.0, 0.0, 0.0);
  glShadeModel(GL_SMOOTH);

  glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, matAmbDiff);
  glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR,            matSpecular);
  glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS,           matShininess);

  glLightfv(GL_LIGHT0, GL_POSITION, light0Pos);
  glLightfv(GL_LIGHT0, GL_DIFFUSE,  light0Diffuse);
  glLightfv(GL_LIGHT0, GL_SPECULAR, light0Spec);

  glLightfv(GL_LIGHT1, GL_POSITION, light1Pos);
  glLightfv(GL_LIGHT1, GL_DIFFUSE,  light0Diffuse);
  glLightfv(GL_LIGHT1, GL_SPECULAR, light0Spec);

  glLightModelfv(GL_LIGHT_MODEL_AMBIENT, lModelAmbient);

  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);
  //glEnable(GL_LIGHT1);
  //glEnable(GL_RESCALE_NORMAL);
  glEnable(GL_NORMALIZE);
  glEnable(GL_DEPTH_TEST);

  // set selection parameters
  
  selBuf = new GLuint[(3*nbMTOCs+1)*10];
  glSelectBuffer((3*nbMTOCs+1)*10, selBuf);
  
  // some nice parameters
  
  glEnable(GL_POINT_SMOOTH);
  //glEnable(GL_LINE_SMOOTH);
  glPointSize(7.0);

  // the timer function that calculates the new positions

  glutTimerFunc(delay, timerFunction, 1);

  qobj = gluNewQuadric();                
}

void initGlut3D()
{
  int nbMTOCs = nuc->nbSegments();
  
  // set the viewpoint
  cameraPos[0]   = 0.0;   // x-position of the camera
  cameraPos[1]   = -10.0; // y-position of the camera
  cameraPos[2]   = 0.0;   // z-position of the camera
  lookAt[0]      = 0.0;   // where to look at
  lookAt[1]      = 0.0;   // where to look at
  lookAt[2]      = 0.0;   // where to look at
  thisSideUp[0]  = 0.0;   // this is up
  thisSideUp[1]  = 0.0;   // this is up
  thisSideUp[2]  = 1.0;   // this is up
  
  glLoadIdentity();
  gluLookAt(cameraPos[0], cameraPos[1], cameraPos[2],
	    lookAt[0], lookAt[1], lookAt[2],
	    thisSideUp[0], thisSideUp[1], thisSideUp[2]);

  // set lighting parameters

  GLfloat matAmbDiff[]    = { 0.0, 0.5, 1.0, 1.0 };
  GLfloat matSpecular[]   = { 1.0, 1.0, 1.0, 1.0 };
  GLfloat matShininess[]  = { 50.0 };
  GLfloat light0Pos[]     = { 2.0, -3.0,  2.0, 0.0 };
  GLfloat light1Pos[]     = {-2.0, -3.0, -2.0, 0.0 };
  GLfloat light0Diffuse[] = { 1.0, 1.0, 1.0, 1.0 };
  GLfloat light0Spec[]    = { 1.0, 1.0, 1.0, 1.0 };
  GLfloat lModelAmbient[] = { 0.2, 0.2, 0.2, 1.0 };
  
  glClearColor(0.0, 0.0, 0.0, 0.0);
  glShadeModel(GL_SMOOTH);

  glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, matAmbDiff);
  glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR,            matSpecular);
  glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS,           matShininess);

  glLightfv(GL_LIGHT0, GL_POSITION, light0Pos);
  glLightfv(GL_LIGHT0, GL_DIFFUSE,  light0Diffuse);
  glLightfv(GL_LIGHT0, GL_SPECULAR, light0Spec);

  glLightfv(GL_LIGHT1, GL_POSITION, light1Pos);
  glLightfv(GL_LIGHT1, GL_DIFFUSE,  light0Diffuse);
  glLightfv(GL_LIGHT1, GL_SPECULAR, light0Spec);

  glLightModelfv(GL_LIGHT_MODEL_AMBIENT, lModelAmbient);
  glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);

  glEnable(GL_LIGHTING);
  glEnable(GL_LIGHT0);
  //glEnable(GL_LIGHT1);
  //glEnable(GL_RESCALE_NORMAL);
  glEnable(GL_NORMALIZE);
  glEnable(GL_DEPTH_TEST);

  // set selection parameters
  
  selBuf = new GLuint[(3*nbMTOCs+1)*10];
  glSelectBuffer((3*nbMTOCs+1)*10, selBuf);
  
  // some nice parameters
  
  glEnable(GL_POINT_SMOOTH);
  //glEnable(GL_LINE_SMOOTH);
  glPointSize(7.0);

  // the timer function that calculates the new positions

  glutTimerFunc(delay, timerFunction, 1);

  qobj = gluNewQuadric();                
}

GLvoid CALLBACK gluErrorCallback()
{
  fprintf(stderr, "A GLU error occured!\n");
  exit(0);
}

GLvoid CALLBACK gluErrorCallback(...)
{
  fprintf(stderr, "A GLU error occured!\n");
  exit(0);
}

void drawTube(GLdouble length)
{
  //gluQuadricCallback(qobj, GLU_ERROR, gluErrorCallback);
  gluQuadricDrawStyle(qobj, GLU_FILL);
  gluQuadricOrientation(qobj, GLU_OUTSIDE);
  gluQuadricNormals(qobj, GLU_SMOOTH);
  gluCylinder(qobj, 0.03, 0.03, length, 5, 1);
  //gluDeleteQuadric(qobj);
}

void drawNucleus()
{
  // draw the nucleus

  glPushMatrix();
  glPushName(0);
  if(DIM == 2) {
    glTranslated((GLfloat)nuc->coord(0), (GLfloat)nuc->coord(1), 0.0);
    gluQuadricDrawStyle(qobj, GLU_SILHOUETTE);
    gluQuadricNormals(qobj, GLU_SMOOTH);
    gluDisk(qobj, 0.0, (GLdouble)nuc->radius(), 100, 1);
  }
  else if(DIM ==3) {
    glTranslated(nuc->coord(0), nuc->coord(1), nuc->coord(2));
    
    int   refpoint1    = DIM*(nuc->nbPoints()-3);
    int   refpoint2    = DIM*(nuc->nbPoints()-2);
    int   refpoint3    = DIM*(nuc->nbPoints()-1);
    
    initRotMat(sphereRotMat);
    for( int ii = 0; ii < 3; ii++) {
      sphereRotMat[ii]   = (nuc->coord(refpoint1+ii) - nuc->coord(ii))/nuc->radius();
      sphereRotMat[ii+4] = (nuc->coord(refpoint2+ii) - nuc->coord(ii))/nuc->radius();
      sphereRotMat[ii+8] = (nuc->coord(refpoint3+ii) - nuc->coord(ii))/nuc->radius();
    }
    
    glMultMatrixd( sphereRotMat );
    glRotated(90, 1, 0, 0);
    glutSolidSphere(nuc->radius(), 20, 20);
    //glutWireSphere(nuc->radius(), 20, 20);
  }
  glPopMatrix();
}

void drawGrabbing()
{
  // draw the grabbing points
  
  // no grabbing point for the reference points
  for(int i = 1; i < nuc->nbPoints()-nuc->nbrefpts; i++) {
    glPushMatrix();
    if(DIM == 2) {
      glTranslated(pg[i*DIM], pg[i*DIM + 1], 0.0);
    }
    else if(DIM ==3) {
      glTranslated((GLfloat)pg[i*DIM], (GLfloat)pg[i*DIM + 1], (GLfloat)pg[i*DIM + 2]);
    }
    glLoadName(i);
    glutSolidSphere(0.1, 10, 10);
    glPopMatrix();
  }
}

void display()
{
  static GLfloat matAmbDiffNuc[]       = { 0.0, 0.5, 1.0, opaque };
  static GLfloat matAmbDiffNucInside[] = { 0.0, 0.0, 0.0, 0.0 };
  static GLfloat matAmbDiffMTOCs[]     = { 1.0, 1.0, 1.0, 1.0 };
  static GLfloat matAmbDiffRef[]       = { 0.0, 0.8, 1.0, 1.0 };
  static GLfloat matAmbDiffHands[]     = { 1.0, 1.0, 0.0, 1.0 };
  matAmbDiffNuc[3] = opaque;
  
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glColor3f(1.0, 1.0, 1.0);
  
  // draw the grabbing points
  glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, matAmbDiffHands);
  drawGrabbing();

  if(DIM == 2) {

    // draw the MTOCs on the nucleus
    
    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, matAmbDiffMTOCs);
    for(int i = 1; i < nuc->nbPoints()-nuc->nbrefpts; i++) {
      glPushMatrix();
      glTranslated(nuc->coord(i*DIM), nuc->coord(i*DIM + 1), 0.0);
      glutSolidSphere(0.1, 10, 10);
      glPopMatrix();
    }
    // a different color for the reference points
    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, matAmbDiffRef);
    for(int i = nuc->nbPoints()-nuc->nbrefpts; i < nuc->nbPoints(); i++) {
      glPushMatrix();
      glTranslated(nuc->coord(i*DIM), nuc->coord(i*DIM + 1), 0.0);
      glutSolidSphere(0.1, 10, 10);
      glPopMatrix();
    }
    
    // draw lines between MTOCs and grabbing points
    
    GLdouble pxg[2];
    
    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, matAmbDiffHands);
    // no line for the reference points
    for(int i = 1; i < nuc->nbPoints()-nuc->nbrefpts; i++) {
      pxg[0] = (GLdouble)(pg[i*DIM] - nuc->coord(i*DIM));
      pxg[1] = (GLdouble)(pg[i*DIM + 1] - nuc->coord(i*DIM + 1));
    
      glPushMatrix();
      glTranslated(nuc->coord(i*DIM), nuc->coord(i*DIM + 1), 0.0);
      glRotated(90.0, -pxg[1], pxg[0], 0.0);
      drawTube((GLdouble)sqrt(pxg[0]*pxg[0] + pxg[1]*pxg[1]));
      glPopMatrix();
    }
  }
  else if(DIM == 3) {
    
    // draw the MTOCs on the nucleus
    
    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, matAmbDiffMTOCs);
    for(int i = 1; i < nuc->nbPoints()-nuc->nbrefpts; i++) {
      glPushMatrix();
      glTranslated(nuc->coord(i*DIM), nuc->coord(i*DIM + 1), nuc->coord(i*DIM + 2));
      glutSolidSphere(0.1, 10, 10);
      glPopMatrix();      
    }
    // a different color for the reference points
    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, matAmbDiffRef);
    for(int i = nuc->nbPoints()-nuc->nbrefpts; i < nuc->nbPoints(); i++) {
      glPushMatrix();
      glTranslated(nuc->coord(i*DIM), nuc->coord(i*DIM + 1), nuc->coord(i*DIM + 2));
      glutSolidSphere(0.1, 10, 10);
      glPopMatrix();      
    }
    
    // draw lines between MTOCs and grabbing points
    
    GLdouble pxg[3];
    GLdouble length;
    
    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, matAmbDiffHands);
    // no line for the reference points
    for(int i = 1; i < nuc->nbPoints()-nuc->nbrefpts; i++) {
      pxg[0] = (GLdouble)(pg[i*DIM] - nuc->coord(i*DIM));
      pxg[1] = (GLdouble)(pg[i*DIM + 1] - nuc->coord(i*DIM + 1));
      pxg[2] = (GLdouble)(pg[i*DIM + 2] - nuc->coord(i*DIM + 2));
      length = norm(pxg);
	
      glPushMatrix();
      glTranslated(nuc->coord(i*DIM), nuc->coord(i*DIM + 1), nuc->coord(i*DIM + 2));
      glRotated(acos(pxg[2]/length)/PI*180.0, -pxg[1], pxg[0], 0.0);
      drawTube(length);
      glPopMatrix();
    }
  }
  
  // draw the nucleus  
  glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, matAmbDiffNuc);
  glMaterialfv(GL_BACK, GL_AMBIENT_AND_DIFFUSE, matAmbDiffNucInside);
  glEnable(GL_BLEND);
  glDepthMask(GL_FALSE);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  //glBlendFunc(GL_SRC_ALPHA, GL_ONE);
  drawNucleus();
  glDepthMask(GL_TRUE);
  glDisable(GL_BLEND);
    
  glFlush();
  glutSwapBuffers();
}

void reshapeWindow(int w, int h)
{
  glViewport(0, 0, (GLsizei) w, (GLsizei) h);
  glGetIntegerv(GL_VIEWPORT, viewPort);

  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
//   glOrtho(-1.0*(GLdouble)glWindowWidth*(GLdouble)w/(GLdouble)screenWidth,
//                (GLdouble)glWindowWidth*(GLdouble)w/(GLdouble)screenWidth,
//  	  -1.0*(GLdouble)glWindowHeight*(GLdouble)h/(GLdouble)screenHeight,
// 	       (GLdouble)glWindowHeight*(GLdouble)h/(GLdouble)screenHeight,
//  	  0.0, 100.0);
  glFrustum(-1.0*(GLdouble)glWindowWidth*(GLdouble)w/(GLdouble)screenWidth,
                 (GLdouble)glWindowWidth*(GLdouble)w/(GLdouble)screenWidth,
            -1.0*(GLdouble)glWindowHeight*(GLdouble)h/(GLdouble)screenHeight,
                 (GLdouble)glWindowHeight*(GLdouble)h/(GLdouble)screenHeight,
             5.0, 15.0);
  glGetDoublev(GL_PROJECTION_MATRIX, projMatrix);
  glMatrixMode(GL_MODELVIEW);
}

void processKey(unsigned char key, int x, int y)
{
  switch(key) {
  case 'q':
    exit(0);
  case 'r':
    glMatrixMode(GL_MODELVIEW);
    translation[0] = 0.0;
    translation[1] = 0.0;
    translation[2] = 0.0;
    zoom           = 1.0;
    initRotMat(rotMat);
    glLoadIdentity();
    gluLookAt(cameraPos[0], cameraPos[1], cameraPos[2],
	      lookAt[0], lookAt[1], lookAt[2],
	      thisSideUp[0], thisSideUp[1], thisSideUp[2]);
    glutTimerFunc(delay, timerFunction, 1);
    glutPostRedisplay();
    break;
  case 'o':
    delay *= 2; break;
  case 'p':
    if ( delay > 1 ) delay /= 2; break;
  case 'K':
    MP.km *= 2;
    printf("set km = %f\n", MP.km);
    break;
  case 'J':
    MP.km /= 2;
    printf("set km = %f\n", MP.km);
    break;
  case 'j':
    MP.km -= 1;
    printf("set km = %f\n", MP.km);
    break;
  case 'k':
    MP.km += 1;
    printf("set km = %f\n", MP.km);
    break;
/*  case 'u':
    noise *= 2;
    printf("set noise = %f\n", noise);
    break;
  case 'y':
    noise /= 2;
    printf("set noise = %f\n", noise);
    break;*/
  case 'e':
    explCalc = !explCalc;
    if(explCalc == 1) printf("switching to explicit calculation.\n");
    else printf("switching to implicit calculation.\n");
    break;
  case 'n':
    normalize = !normalize;
    break;
  case 'd':
    opaque += 0.1;
    if(opaque > 1.0) opaque = 1.0;
    break;
  case 'D':
    opaque -= 0.1;
    if(opaque < 0.0) opaque = 0.0;
    break;
  case 'c':
    correct = !correct;
    if(correct == 1) printf("switching to correct linearization (only implicit mode).\n");
    else printf("switching to bad linearization     (only implicit mode).\n");
    break;
  case 't':
    ppprime = !ppprime;
    if(ppprime == 1) printf("using p*(1+p')                     (only implicit mode).\n");
    else printf("using p+p'                         (only implicit mode).\n");
    break;
  default:
    printf("ignored key %c %i %i\n", key, x, y);
    break;
  }
}

void processHits (GLint hits, GLuint buffer[])
{
   unsigned int i, j;
   GLuint names, *ptr;

   printf ("hits = %d\n", hits);
   ptr = (GLuint *) buffer;
   for (i = 0; i < (unsigned int)hits; i++) {	/*  for each hit  */
      names = *ptr;
      printf (" number of names for hit = %d\n", names); ptr++;
      printf("  z1 is %g;", (float) *ptr/0x7fffffff); ptr++;
      printf(" z2 is %g\n", (float) *ptr/0x7fffffff); ptr++;
      printf ("   the names on the stack are ");
      for (j = 0; j < names; j++) {	/*  for each name */
         printf ("%d ", *ptr); ptr++;
      }
      printf ("\n");
   }
}

void processMouse(int button, int state, int x, int y)
{
  GLdouble clickNorm;
  
  if(state != GLUT_DOWN) { return; }

  mouse_action = mouseAction( button, glutGetModifiers() );
  if ( mouse_action == MOUSE_NOTHING ) return;
  
  glGetDoublev(GL_MODELVIEW_MATRIX, mvMatrix);
  gluUnProject((GLdouble)x, ((GLdouble)viewPort[3] - (GLdouble)y - 1), 0.0,
  	       mvMatrix, projMatrix, viewPort,
  	       clickPointVec, clickPointVec+1, clickPointVec+2);
  
  switch(mouse_action) {
  case MOUSE_ROTATE:
    nclickPointVec[0] = clickPointVec[0] - translation[0];
    nclickPointVec[1] = clickPointVec[1] - translation[1];
    nclickPointVec[2] = clickPointVec[2] - translation[2];
    clickNorm = norm(nclickPointVec);
    nclickPointVec[0] /= clickNorm;
    nclickPointVec[1] /= clickNorm;
    nclickPointVec[2] /= clickNorm;
    copyMat(rotMat, rotMatSave);
    break;
  case MOUSE_SPIN:
    printf("Spinning is not yet implemented!\n");
    break;
  case MOUSE_TRANSLATE:
    transSave[0] = translation[0];
    transSave[1] = translation[1];
    transSave[2] = translation[2];
    break;
  case MOUSE_ZOOM:
    zoomSave = zoom;
    mouseX   = x;
    break;
  case MOUSE_GRAB:
    GLdouble clickZMinusOne[3];
    GLdouble intoTheScreen[3];
    GLdouble factor;
    
    // this is only needed to do 2d-grabbing in 3d
    gluUnProject((GLdouble)x, ((GLdouble)viewPort[3] - (GLdouble)y - 1), -1.0,
		 mvMatrix, projMatrix, viewPort,
		 clickZMinusOne, clickZMinusOne+1, clickZMinusOne+2);
    intoTheScreen[0] = clickPointVec[0] - clickZMinusOne[0];
    intoTheScreen[1] = clickPointVec[1] - clickZMinusOne[1];
    intoTheScreen[2] = clickPointVec[2] - clickZMinusOne[2];
    if( intoTheScreen[2] != 0. ) {
      factor = -clickPointVec[2]/intoTheScreen[2];
    }
    else {
      factor = 0.0;
      printf("division by zero!\n");
      //exit(0);
    }
    //

    glGetIntegerv(GL_VIEWPORT, viewPort);

    glRenderMode(GL_SELECT);
    glInitNames();

    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();
    gluPickMatrix((GLdouble)x, (GLdouble)(viewPort[3] - y), 3.0, 3.0, viewPort);
    glFrustum(-1.0*(GLdouble)glWindowWidth*(GLdouble)viewPort[2]/(GLdouble)screenWidth,
	           (GLdouble)glWindowWidth*(GLdouble)viewPort[2]/(GLdouble)screenWidth,
	      -1.0*(GLdouble)glWindowHeight*(GLdouble)viewPort[3]/(GLdouble)screenHeight,
	           (GLdouble)glWindowHeight*(GLdouble)viewPort[3]/(GLdouble)screenHeight,
	      5.0, 15.0);

    // switch to modelview for drawing the objects that are to be selected
    glMatrixMode(GL_MODELVIEW);
    drawNucleus();
    drawGrabbing();

    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
    glFlush();

    hits = glRenderMode(GL_RENDER);

    // set the new position of the selected object
    if(hits) {
      selection = selBuf[3];
      //printf("Number: %d\n", selection);
      
      if(selection == 0) {
	if(DIM == 2) {
	  nuc->pts(0) = clickPointVec[0];
	  nuc->pts(1) = clickPointVec[1];
	}
	else if(DIM == 3) {
	  nuc->pts(0) = clickPointVec[0];
	  nuc->pts(1) = clickPointVec[1];
	  nuc->pts(2) = clickPointVec[2];
	}
      }
      else {
	if(DIM == 2) {
	  pg[DIM*selection] = clickPointVec[0] + factor*intoTheScreen[0];
	  pg[DIM*selection + 1] = clickPointVec[1] + factor*intoTheScreen[1];
	}
	else if(DIM == 3) {
	  pg[DIM*selection] = clickPointVec[0];
	  pg[DIM*selection + 1] = clickPointVec[1];
	  pg[DIM*selection + 2] = clickPointVec[2];
	}
      }
    }
    //processHits(hits, selBuf);
    glMatrixMode(GL_MODELVIEW);

    break;
  default:
    break;
  }
}

void processMouseMotion(int x, int y)
{
  GLdouble w[3];
  //GLdouble axis[3];
  //GLdouble angle;

  switch(mouse_action) {
  case MOUSE_ROTATE:
    GLdouble motionVec[3];
    GLdouble axis[3];
    //GLdouble axisNorm;
    
    gluUnProject((GLdouble)x, ((GLdouble)viewPort[3] - (GLdouble)y - 1), 0.0,
		 mvMatrix, projMatrix, viewPort, w, w+1, w+2);
    motionVec[0] = w[0] - clickPointVec[0];
    motionVec[1] = w[1] - clickPointVec[1];
    motionVec[2] = w[2] - clickPointVec[2];
    axis[0]  = nclickPointVec[1] * motionVec[2] - nclickPointVec[2] * motionVec[1];
    axis[1]  = nclickPointVec[2] * motionVec[0] - nclickPointVec[0] * motionVec[2];
    axis[2]  = nclickPointVec[0] * motionVec[1] - nclickPointVec[1] * motionVec[0];
    //axisNorm = norm(axis);
    //angle = 0.5*axisNorm/PI*180;
    // Obviously we don't have to normalize the axis for glRotate
    //   axis[0] /= axisNorm;
    //   axis[1] /= axisNorm;
    //   axis[2] /= axisNorm;

    createRotMat(axis, rotMat);
    multMat(rotMatSave, rotMat);
    break;
  case MOUSE_SPIN:
    break;
  case MOUSE_TRANSLATE:
    gluUnProject((GLdouble)x, ((GLdouble)viewPort[3] - (GLdouble)y - 1), 0.0,
		 mvMatrix, projMatrix, viewPort, w, w+1, w+2);
    translation[0] = transSave[0] + clickPointVec[0] - w[0];
    translation[1] = transSave[1] + clickPointVec[1] - w[1];
    translation[2] = transSave[2] + clickPointVec[2] - w[2];
    break;
  case MOUSE_ZOOM:
    GLdouble zFactor;
    zFactor = 1.0 + (GLdouble)(x - mouseX)/(GLdouble)viewPort[2];
    if( zFactor > 0 ) zoom = zoomSave * zFactor;
    break;
  case MOUSE_GRAB:
    if(hits) {
      GLdouble clickZMinusOne[3];
      GLdouble intoTheScreen[3];
      GLdouble factor;
      
      gluUnProject((GLdouble)x, ((GLdouble)viewPort[3] - (GLdouble)y - 1), 0.0,
		   mvMatrix, projMatrix, viewPort, w, w+1, w+2);
      
      // this is only needed to do 2d-grabbing in 3d
      gluUnProject((GLdouble)x, ((GLdouble)viewPort[3] - (GLdouble)y - 1), -1.0,
		   mvMatrix, projMatrix, viewPort,
		   clickZMinusOne, clickZMinusOne+1, clickZMinusOne+2);
      intoTheScreen[0] = w[0] - clickZMinusOne[0];
      intoTheScreen[1] = w[1] - clickZMinusOne[1];
      intoTheScreen[2] = w[2] - clickZMinusOne[2];
      if( intoTheScreen[2] != 0. ) {
	factor = -w[2]/intoTheScreen[2];
      }
      else {
	factor = 0.0;
	printf("division by zero!\n");
      }
      //
      
      if(selection == 0) {
	if(DIM == 2) {
	  nuc->pts(0) = w[0];
	  nuc->pts(1) = w[1];
	}
	else if(DIM == 3) {
	  nuc->pts(0) = w[0];
	  nuc->pts(1) = w[1];
	  nuc->pts(2) = w[2];
	}
      }
      else {
	if(DIM == 2) {
	  pg[DIM*selection] = w[0] + factor*intoTheScreen[0];
	  pg[DIM*selection + 1] = w[1] + factor*intoTheScreen[1];
	}
	else if(DIM == 3) {
	  pg[DIM*selection] = w[0];
	  pg[DIM*selection + 1] = w[1];
	  pg[DIM*selection + 2] = w[2];
	}
      }
    }
    break;
  default:
    break;
  }
  
  glLoadIdentity();
  gluLookAt(cameraPos[0], cameraPos[1], cameraPos[2],
	    lookAt[0], lookAt[1], lookAt[2],
	    thisSideUp[0], thisSideUp[1], thisSideUp[2]);
  glMultMatrixd( rotMat );
  glScaled(zoom, zoom, zoom);
  glTranslated(-translation[0], -translation[1], -translation[2]);
  glutPostRedisplay();
  
}

int main(int argc, char** argv)
{
  // first do the GLUT initialisation
  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
  glutInitWindowSize(screenWidth, screenHeight);
  //glutInitWindowPosition(200, 200);
  glutCreateWindow(argv[0]);

  // initialize the nucleus...
  initialise();
    
  // initialise GLUT...
  switch( DIM ) {
  case 2:
    initGlut2D();
    break;
  default:
    initGlut3D();
  }
  initRotMat(rotMat);

  glutDisplayFunc(display);
  glutReshapeFunc(reshapeWindow);
  glutKeyboardFunc(processKey);
  glutMouseFunc(processMouse);
  glutMotionFunc(processMouseMotion);

  // ...and enter the main loop
  glutMainLoop();
  cleanUp(); // this is actually never called

  return 0;
}
