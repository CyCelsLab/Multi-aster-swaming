//RCS: $Id: testsolved.cc,v 2.1 2005/01/10 15:59:59 foethke Exp $
//--------------------------------------------------------------------------
//
// nedelec@embl-heidelberg.de  April 2003

#include <stdlib.h>
#include <stdio.h>
#include "glut.h"

#include "types.h"
#include "random.h"
#include "cblas.h"
#include "clapack.h"
#include "smath.h"

//a point in space:
GLdouble Gx=1, Gy=1, Gz;

real  radius = 2;                    // the radius of the circle
int   nbp    = 3;                   // number of points
int   dim    = 2;
int   vecDim = dim*(nbp + 1);
//real*    pf;
//real*    pn;
//real*    pc;

real* px;
real* pg;

//rate of calculation:
int  delay     = 50;
int  normalize = 1;
bool correct   = 1;
bool ppprime   = 1;
bool explCalc  = 1;

real dt = 1./16.;
real mu = 1.;
real km = 1.;
real noise = 1.;

void initialize() {

  real randAngle;

  px = new real[vecDim];
  pg = new real[vecDim];

  // initialize the points on the circle and the hands
  px[0] = 0.;
  px[1] = 0.;
  px[2] = px[0] + radius;
  px[3] = px[1];

  pg[0] = 0.;
  pg[1] = 0.;
  pg[2] = px[0] + radius + 3.;
  pg[3] = px[1];

//   px[4] = px[0];
//   px[5] = px[1] + radius;
//   px[6] = px[0] - radius;
//   px[7] = px[1];
//   pg[4] = px[0];
//   pg[5] = px[1] + radius + 3.;
//   pg[6] = px[0] - radius - 3.;
//   pg[7] = px[1];

  for(int i = 2; i <= nbp; i++) {
    randAngle     = 2.*PI*RNG.preal();
    px[i*dim]     = radius*cos(randAngle);
    px[i*dim + 1] = radius*sin(randAngle);
    pg[i*dim]     = (radius + 3.)*cos(randAngle);
    pg[i*dim + 1] = (radius + 3.)*sin(randAngle);
  }

}

void timerFunction(int value)
{

  real noiseVec[vecDim];
  real axis[dim];
  real dummy[vecDim];
  real j[nbp*vecDim];
  real jjt[nbp*nbp];
  real jjtInvj[nbp*vecDim];
  real p[vecDim*vecDim];
  real norm;
  real* unit = new real(1.);
  //  int* info  = new int;
  int info;

  // set the jacobian matrix j
  for(int i = 0; i < (nbp*vecDim); i++) {
    j[i] = 0.;
  }
  for(int i = 1; i <= nbp; i++) {  
    for(int k = 0; k < dim; k++) {
      j[(i-1) + k*nbp] = -2.*(px[dim*i+k] - px[k]);
    }
  }
  for(int i = 1; i <= nbp; i++) {  
    for(int k = 0; k < dim; k++) {
      j[(i-1) + (i*dim+k)*nbp] =  2.*(px[dim*i+k] - px[k]);
    }
  }

  // set (jj^T)
  for(int k = 0; k < nbp; k++) {
    for(int i = 0; i <= k; i++) {
      jjt[i+k*nbp] = 0.;
    }
  }
  for(int k = 0; k < nbp; k++) {
    for(int i = 0; i <= k; i++) {
      for(int l = 0; l < dim; l++) {
	jjt[i+k*nbp] += 4.*(px[(i+1)*dim + l] - px[l])*(px[(k+1)*dim + l] - px[l]);	
      }
    }
  }
  for(int i = 0; i < nbp; i++) {
    jjt[i+i*nbp] *= 2.;
  }
  
  // calculate jjt^(-1) from jjt
  lapack_xpotrf('u',nbp,jjt,nbp,&info);
  lapack_xpotri('u',nbp,jjt,nbp,&info);
  
  // calculate the projector P
  blas_xsymm('l','u',nbp,vecDim,1.,jjt,nbp,j,nbp,0,jjtInvj,nbp);
  blas_xgemm('t','n',vecDim,vecDim,nbp,1.,j,nbp,jjtInvj,nbp,0,p,vecDim);

  // calculate 1-P
  blas_xscal(vecDim*vecDim,-1.,p,1);
  blas_xaxpy(vecDim,1.,unit,0,p,vecDim+1);
  
  // initialize noise
  for(int i = 0; i < vecDim; i++) {
    noiseVec[i] = RNG.gauss()*noise;
    //noiseVec[i] = 0.;
  }

  if( explCalc ) {

    real fEx[vecDim];
    real fRes[vecDim];
    real a[vecDim*vecDim];
    
    // initialize Matrix a
    for(int i = 0; i < (vecDim*vecDim); i++) {
      a[i] = 0.;
    }
    // starts at dim because we don't have forces on the center point
    for(int i = dim; i < vecDim; i++) {
      a[i + i*vecDim] = - km;
    }
        
    // calculate external forces
    blas_xgemv('n',vecDim,vecDim,1.,a,vecDim,px,1,0.,dummy,1);
    for(int i = 0; i < vecDim; i++) {  
      fEx[i] = dummy[i] + km*pg[i] + noiseVec[i];
    }
    
    // calculate legal forces
    blas_xgemv('n',vecDim,vecDim,1.,p,vecDim,fEx,1,0.,fRes,1);
    
  
    // set the new positions
    for(int i = 0; i < vecDim; i++) {
      px[i] += fRes[i]*dt;
    }

  }
  else {
    
    int  ipiv[vecDim + 1];
    //int* lwork = new int;
    //*lwork   = vecDim*32;
    int  lwork = vecDim*32;
    real work[lwork];


    if( correct ) {

      real jtjjtInvt[vecDim*nbp];
      real djdxcp0[nbp*vecDim];
      real djdxip0[nbp*vecDim];
      real djtdxcjjtInvj[vecDim*vecDim];
      real djtdxijjtInvj[vecDim*vecDim];
      real pPrimeSlice[vecDim*vecDim];
      real pPrimeB[vecDim*vecDim];
      real pPrimeBt[vecDim*vecDim];

      //blas_xscal(vecDim,km,pg,1);
      blas_xaxpy(vecDim,km,pg,1,noiseVec,1);
    
      // set jtjjtInvt
      blas_xsymm('l','u',nbp,vecDim,1.,jjt,nbp,j,nbp,0,jtjjtInvt,nbp);
      
      for(int d = 0; d < dim; d++) {

	// set djdxcp0
	for(int i = 0; i < nbp; i++) {
	  for(int k = 0; k < vecDim; k++) {
	    djdxcp0[i + k*nbp] = 2.*(p[d + k*vecDim] - p[(i+1)*dim+d + k*vecDim]);
	  }
	}

	// set djtdxcjjtInvj
	for(int i = 0; i < vecDim; i++) {
	  for(int k = 0; k < vecDim; k++) {
	    djtdxcjjtInvj[i + k*vecDim] = 0.;
	  }
	}
	for(int i = 0; i < nbp; i++) {
	  for(int k = 0; k < vecDim; k++) {
	    djtdxcjjtInvj[(i+1)*dim+d + k*vecDim]  = -2.*jjtInvj[i + k*nbp];
	    djtdxcjjtInvj[d           + k*vecDim] -= djtdxcjjtInvj[(i+1)*dim+d + k*vecDim];
	  }
	}
	
	// calculate jtjjtInv * djdxcp0
	blas_xgemm('t','n',vecDim,vecDim,nbp,1.,jtjjtInvt,nbp,djdxcp0,nbp,0,pPrimeSlice,vecDim);
	
	// calculate slice of pPrime = p0 * djtdxcjjtInvj + jtjjtInv * djdxcp0
	blas_xgemm('n','n',vecDim,vecDim,vecDim,1.,p,vecDim,djtdxcjjtInvj,vecDim,1.,pPrimeSlice,vecDim);
	
	// calculate pPrimeSlice * B
	blas_xgemv('n',vecDim,vecDim,1.,pPrimeSlice,vecDim,noiseVec,1,0.,(pPrimeB + d*vecDim),1);

	// now do it for the other slices that derive from the center slice
	for(int i = 0; i < nbp; i++) {

	  // set djdxip0
	  for(int l = 0; l < nbp; l++) {
	    for(int k = 0; k < vecDim; k++) {
	      djdxip0[l + k*nbp] = 0.;
	    }
	  }
	  for(int k = 0; k < vecDim; k++) {
	    djdxip0[i + k*nbp] = -1.*djdxcp0[i + k*nbp];
	  }
	  	
	  // set djtdxijjtInvj
	  for(int l = 0; l < vecDim; l++) {
	    for(int k = 0; k < vecDim; k++) {
	      djtdxijjtInvj[l + k*vecDim] = 0.;
	  }
	  }
	  for(int k = 0; k < vecDim; k++) {
	    djtdxijjtInvj[d           + k*vecDim] =     djtdxcjjtInvj[(i+1)*dim+d + k*vecDim];
	    djtdxijjtInvj[(i+1)*dim+d + k*vecDim] = -1.*djtdxcjjtInvj[(i+1)*dim+d + k*vecDim];
	  }
	
	  // calculate jtjjtInv * djdxip0
	  blas_xgemm('t','n',vecDim,vecDim,nbp,1.,jtjjtInvt,nbp,djdxip0,nbp,0,pPrimeSlice,vecDim);
	  
	  // calculate slice of pPrime = p0 * djtdxijjtInvj + jtjjtInv * djdxip0
	  blas_xgemm('n','n',vecDim,vecDim,vecDim,1.,p,vecDim,djtdxijjtInvj,vecDim,1.,pPrimeSlice,vecDim);
	  
	  // calculate pPrimeSlice * B
	  blas_xgemv('n',vecDim,vecDim,1.,pPrimeSlice,vecDim,noiseVec,1,0.,(pPrimeB+((i+1)*dim+d)*vecDim),1);
	  
	}
      }

      if( ppprime ) {
	blas_xgemm('n','t',vecDim,vecDim,vecDim,dt*mu,p,vecDim,pPrimeB,vecDim,0.,pPrimeBt,vecDim);
      }
      else {

	// beware, pPrimeB is actually pPrimeB transposed!!
	for(int i = 0; i < vecDim; i++) {
	  for(int k = 0; k < vecDim; k++) {
	    pPrimeBt[i + k*vecDim] = dt*mu*pPrimeB[k + i*vecDim];
	  }
	}
      }
      
      blas_xgemv('n',vecDim,vecDim,dt*mu,p,vecDim,noiseVec,1,0.,dummy,1);
      blas_xaxpy(vecDim,1.,unit,0,pPrimeBt,vecDim+1);
      blas_xgemv('n',vecDim,vecDim,1.,pPrimeBt,vecDim,px,1,1.,dummy,1);
      
      blas_xscal(vecDim*vecDim,dt*mu*km,p,1);
      // clean up forces on the center point
      for(int k = 0; k < dim; k++) {
	for(int i = 0; i < vecDim; i++) {
	  p[i + k*vecDim] = 0.;
	}
      }
      
      blas_xaxpy(vecDim*vecDim,1.,pPrimeBt,1,p,1);
      lapack_xgetrf(vecDim,vecDim,p,vecDim,ipiv,&info);
      lapack_xgetri(vecDim,p,vecDim,ipiv,work,lwork,&info);
      
      blas_xgemv('n',vecDim,vecDim,1.,p,vecDim,dummy,1,0,px,1);
      
    }
    
    else {
      
      blas_xaxpy(vecDim,km,pg,1,noiseVec,1);
      blas_xgemv('n',vecDim,vecDim,dt*mu,p,vecDim,noiseVec,1,1.,px,1);
      
      blas_xscal(vecDim*vecDim,dt*mu*km,p,1);
      // clean up forces on the center point
      for(int k = 0; k < dim; k++) {
	for(int i = 0; i < vecDim; i++) {
	  p[i + k*vecDim] = 0.;
	}
      }

      blas_xaxpy(vecDim,1.,unit,0,p,vecDim+1);
      lapack_xgetrf(vecDim,vecDim,p,vecDim,ipiv,&info);
      lapack_xgetri(vecDim,p,vecDim,ipiv,work,lwork,&info);
      
      blas_xgemv('n',vecDim,vecDim,1.,p,vecDim,px,1,0,dummy,1);
      for(int i = 0; i < vecDim; i++) {
	px[i] = dummy[i];
      }    
    }
  }
    
  // calculate center of gravity ...
  int nucWeight = 500;
  real* cgrav = new real[dim];
  cgrav[0] = nucWeight*px[0];
  cgrav[1] = nucWeight*px[1];
  for(int j = 1; j <= nbp; j++) {
    for(int i = 0; i < dim; i++) {
      cgrav[i] += px[i+j*dim];
    }
  }
  for(int i = 0; i < dim; i++) {
    cgrav[i] /= (nbp+nucWeight);
  }

//   printf("px[0]:    %f, px[1]:    %f\n", px[0], px[1]);
//   printf("px[2]:    %f, px[3]:    %f\n", px[2], px[3]);  
//   printf("cgrav[0]: %f, cgrav[1]: %f\n\n", cgrav[0], cgrav[1]);

  for(int i = 0; i < dim; i++) {
    px[i] = cgrav[i];
  }

  // ... and ged rid of finite-step errors
  for(int j = 1; j <= nbp; j++) {
    for(int i = 0; i < dim; i++) {
      axis[i] = px[i+j*dim] - px[i];
    }
    norm = blas_xnrm2(dim,axis,1);
    for(int i = 0; i < dim; i++) {
      axis[i]    /= norm;
      px[i+j*dim] = px[i] + axis[i]*radius;
    }
  }

  glutPostRedisplay();
  glutTimerFunc( delay, timerFunction, 1);
}


void display()
{
  glClear(GL_COLOR_BUFFER_BIT);

  // the circle
  glColor3f(0,0,1.0);
  glBegin(GL_LINE_LOOP);
  for( real ii = 0 ; ii < 6.28; ii+=0.0314 )
    glVertex2f( px[0] + radius*cos(ii), px[1] + radius*sin(ii) );
  glEnd();


  glPointSize(7.0);
  glBegin(GL_POINTS);

  // the points on the circle
  glColor3f(1.0, 1.0, 1.0);
  for(int i = 1; i <= nbp; i ++) {
    glVertex2f(px[i*dim], px[i*dim + 1]);
  }

  // the grab points
  glColor3f(0.0, 0.0, 1.0);
  for(int i = 1; i <= nbp; i ++) {
    glVertex2f(pg[i*dim], pg[i*dim + 1]);
  }

  //glColor3f(1.0, 0.0, 0.0);
  //glVertex2f(pf[0], pf[1]);

  //glColor3f(0.0, 1.0, 0.0);
  //glVertex2f(pc[0], pc[1]);
  
  //glColor3f(0.0, 0.0, 1.0);
  //glVertex2f(pn[0], pn[1]);
  glEnd();


  // lines between circle- and grab-points
  glColor3f(1.0,1.0,1.0);
  for(int i = 1; i <= nbp; i ++) {
    glBegin(GL_LINE_LOOP);
    glVertex2f(px[i*dim], px[i*dim + 1]);
    glVertex2f(pg[i*dim], pg[i*dim + 1]);
    glEnd();
  }

  // the radius
  glColor3f(0.0,1.0,1.0);
  glBegin(GL_LINE_LOOP);
  glVertex2f(px[2], px[3]);
  glVertex2f(px[0], px[1]);
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
    km *= 2;
    printf("set km = %f\n", km);
    break;
  case 'J':
    km /= 2;
    printf("set km = %f\n", km);
    break;
  case 'k':
    km += 1;
    printf("set km = %f\n", km);
    break;
  case 'j':
    km -= 1;
    printf("set km = %f\n", km);
    break;
    
  case 'i':
    noise *= 2;
    printf("set noise = %f\n", noise);
    break;
  case 'u':
    noise /= 2;
    printf("set noise = %f\n", noise);
    break;
    
  case 'e':
    explCalc = !explCalc;
    if(explCalc == 1) printf("switching to explicit calculation.\n");
    else printf("switching to implicit calculation.\n");
    break;

  case 'n':
    normalize = !normalize;
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
    
//   case 'z':
//     px.set( 1, 0, 0 );
//     pg.set( 5, 0, 0 );
//     break;
    
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
    //pg.set(Gx, Gy, Gz);
    pg[2] = Gx;
    pg[3] = Gy;
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
    //pg.set(Gx, Gy, Gz);
    pg[2] = Gx;
    pg[3] = Gy;    
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
  initialize();

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
