//RCS: $Id: play_draw3.cc,v 2.19 2005/04/22 12:28:11 nedelec Exp $


#include "globjects3d.h"


using namespace glObjects3d;


//the "real" radius of mts is ca 12nm
//to go from mtwidth and mtsize (2 pixels) to 12nm we devide by 167
//const real sizeFactor = 167.;
const real sizeFactor = 80.;

// pointers to quadric objects
GLUquadricObj* qobj;
GLUExtQuadric* extqobj;

// width of tubes used to draw the box and the nucleus
real lineWidth = 0.04;


//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//functions for translation and rotation
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

real glRelocateTube(const real* v1, const real* v2)
// translate and rotate the world in order to draw a tube
// that extends from v1 to v2, returns the length of the tube
{
  static Vecteur tubeDirection;
  static real    tubeLength;

  //calculate the rotation angle for the cylinders
  tubeDirection.XX = v2[0] - v1[0];
  #if( DIM > 1 )
  tubeDirection.YY = v2[1] - v1[1];
  #endif
  #if( DIM > 2 )
  tubeDirection.ZZ = v2[2] - v1[2];
  #endif
  tubeLength       = tubeDirection.norm();

  glTranslatev( v1 );
  #if( DIM == 2 )
  glRotated(90,-tubeDirection.YY,tubeDirection.XX,0);
  //this second rotation is necessary in 2D to make the front of the tube
  //actually face the front of the screen, otherwise we see a part of the back
  if(tubeDirection.YY > 0)
    glRotated( acos(tubeDirection.XX/tubeLength)/PI*180., 0.0, 0.0, 1.0 );
  else
    glRotated( -acos(tubeDirection.XX/tubeLength)/PI*180., 0.0, 0.0, 1.0 );
  #elif( DIM == 3 )
  glRotated(acos(tubeDirection.ZZ/tubeLength)/PI*180.,-tubeDirection.YY,tubeDirection.XX,0);
  #endif

  return(tubeLength);
}


//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//functions to draw specific objects like the bounding box, a half-sphere,
//an oval, an enhanced cylinder etc.
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------


void displayBoxSection3(const Space * mSpace, int dim, float pos, float step)
// This is a direct adaptation of displayBoxSection1 in play_draw1.cc.
// Instead of lines we simply draw little cylinders. Drawing a lot of
// cylinders is slow and this function should only be used if there's
// no better way to draw the box.
{
  Vecteur q, p( pos, pos, pos );
  Vecteur boxStart;
  Vecteur tubeStart;
  real    tubeLength;
  bool    odd = 0;
  int     xx = ( dim + 1 ) % DIM;
  int     yy = ( xx + 1 ) % DIM;
  real    xs = mSpace->getBoundingRect()[xx];
  real    ys = mSpace->getBoundingRect()[yy];

  p[yy] = ys;
  for( real a = -xs; a < xs; a += step ) {
    p[xx] = a;
    mSpace->project(p, q);
    if ( odd ) {
      glPushMatrix();
      tubeLength = glRelocateTube( tubeStart, q );
      gluCylinder(qobj, lineWidth, lineWidth, tubeLength, 5, 5);
      glPopMatrix();
    }
    else
      boxStart = q;
    tubeStart = q;
    odd = 1;
  };
  p[xx] = xs;
  for( real a = -ys; a < ys; a += step ) {
    p[yy] = -a;
    mSpace->project(p, q);
    if ( odd ) {
      glPushMatrix();
      tubeLength = glRelocateTube( tubeStart, q );
      gluCylinder(qobj, lineWidth, lineWidth, tubeLength, 5, 5);
      glPopMatrix();
    }
    else
      boxStart = q;
    tubeStart = q;
    odd = 1;
  };
  p[yy] = -ys;
  for( real a = -xs; a < xs; a += step ) {
    p[xx] = -a;
    mSpace->project(p, q);
    if ( odd ) {
      glPushMatrix();
      tubeLength = glRelocateTube( tubeStart, q );
      gluCylinder(qobj, lineWidth, lineWidth, tubeLength, 5, 5);
      glPopMatrix();
    }
    else
      boxStart = q;
    tubeStart = q;
    odd = 1;
  };
  p[xx] = -xs;
  for( real a = -ys; a < ys; a += step ) {
    p[yy] = a;
    mSpace->project(p, q);
    if ( odd ) {
      glPushMatrix();
      tubeLength = glRelocateTube( tubeStart, q );
      gluCylinder(qobj, lineWidth, lineWidth, tubeLength, 5, 5);
      glPopMatrix();
    }
    else
      boxStart = q;
    tubeStart = q;
    odd = 1;
  };
  //close the loop:
  if ( odd ) {
    glPushMatrix();
    tubeLength = glRelocateTube( tubeStart, boxStart );
    gluCylinder(qobj, lineWidth, lineWidth, tubeLength, 5, 5);
    glPopMatrix();
  }
}


//-----------------------------------------------------------------------------
void displayBox3(const Space * space)
// This is a copy of displayBox1 in play_draw1.cc.
// Normal glColor commands are ignored when lightning is enabled,
// so they are useless here!
///\todo exchange glColor with glMaterialProp
///\todo exchange points and lines with real 3D objects
{
  glColor( PP.boxcolor );
  glLineWidth( PP.boxwidth );

  switch( space->getShape() ) {

  case SHAPE_SQUARE:
  case SHAPE_PERIODIC:
  case SHAPE_STRIP:
    glBegin(GL_LINE_LOOP);
    glVertex3(  MP.boxsize[0],  MP.boxsize[1], MP.boxsize[2] );
    glVertex3(  MP.boxsize[0], -MP.boxsize[1], MP.boxsize[2] );
    glVertex3( -MP.boxsize[0], -MP.boxsize[1], MP.boxsize[2] );
    glVertex3( -MP.boxsize[0],  MP.boxsize[1], MP.boxsize[2] );
    glEnd();
#if ( DIM > 2 )
    glBegin(GL_LINE_LOOP);
    glVertex3(  MP.boxsize[0],  MP.boxsize[1], -MP.boxsize[2] );
    glVertex3(  MP.boxsize[0], -MP.boxsize[1], -MP.boxsize[2] );
    glVertex3( -MP.boxsize[0], -MP.boxsize[1], -MP.boxsize[2] );
    glVertex3( -MP.boxsize[0],  MP.boxsize[1], -MP.boxsize[2] );
    glEnd();
#endif
    break;

  case SHAPE_SPHERE:
    glBegin(GL_LINE_LOOP);
    for( real aa = 0; aa < 6.28; aa += 0.05 )
      glVertex3( MP.boxsize[0]*cos(aa), MP.boxsize[0]*sin(aa), 0 );
    glEnd();
#if (DIM > 2)
    glBegin(GL_LINE_LOOP);
    for( real aa = 0; aa < 6.28; aa += 0.1 )
      glVertex3( 0, MP.boxsize[0]*cos(aa), MP.boxsize[0]*sin(aa) );
    glEnd();
    glBegin(GL_LINE_LOOP);
    for( real aa = 0; aa < 6.28; aa += 0.1 )
      glVertex3( MP.boxsize[0]*cos(aa), 0, MP.boxsize[0]*sin(aa) );
    glEnd();
#endif
    break;

  case SHAPE_OVAL:
  case SHAPE_CYLINDER:
    displayBoxSection3( space, 2, 0, 0.1 );
    if ( DIM == 3 ) {
      for( real Z = +0; Z <= +MP.boxsize[0]; Z += 1.0 )
        displayBoxSection3( space, 0, Z, 0.1 );
      for( real Z = -1; Z >= -MP.boxsize[0]; Z -= 1.0 )
        displayBoxSection3( space, 0, Z, 0.1 );
      displayBoxSection3( space, 1, 0, 0.1 );
    }
  break;

  case SHAPE_BANANA: {
    real center[2];
    real angle, rad;
    real bRadius = MP.boxsize[0]/MP.boxsize[2];

    //the rightmost point on the backbone of the banana
    //(for the leftmost point we use -center[0] below)
    center[0] = bRadius*sin(MP.boxsize[2]);
    center[1] = bRadius*cos(MP.boxsize[2]) - bRadius;

    glBegin(GL_LINE_LOOP);
    //the upper bow from right to left (+x to -x)
    for( real aa = -MP.boxsize[2]; aa <=  MP.boxsize[2]; aa += 0.01 )
      glVertex3( -(bRadius+MP.boxsize[1])*sin(aa), (bRadius+MP.boxsize[1])*cos(aa)-bRadius, 0. );
    //the left cap
    for( real aa = MP.boxsize[2];  aa <= MP.boxsize[2]+PI; aa += 0.01 )
      glVertex3( -center[0] - MP.boxsize[1]*sin(aa), center[1] + MP.boxsize[1]*cos(aa), 0. );
    //the lowe bow from left to right
    for( real aa = MP.boxsize[2];  aa >= -MP.boxsize[2]; aa -= 0.01 )
      glVertex3( -(bRadius-MP.boxsize[1])*sin(aa), (bRadius-MP.boxsize[1])*cos(aa)-bRadius, 0. );
    //the right cap
    for( real aa = (PI-MP.boxsize[2]);  aa <= (PI-MP.boxsize[2])+PI; aa += 0.01 )
      glVertex3( center[0] - MP.boxsize[1]*sin(aa), center[1] + MP.boxsize[1]*cos(aa), 0. );
    glEnd();

    if ( DIM == 3 ) {

      real xxmax, yymax;

      //display rings ~ every 1.2 um
      int nbSegments = (int)floor( 2.*MP.boxsize[0] / 1.2 );
      for( int ii = 0; ii <= nbSegments; ii++ ) {
        angle = (-MP.boxsize[0] + 2.*MP.boxsize[0]/nbSegments*ii)/bRadius;
        glBegin(GL_LINE_LOOP);
        for( real aa = 0;  aa < 2.*PI; aa += 0.01 ) {
          rad = bRadius + MP.boxsize[1]*cos(aa);
          glVertex3( rad*sin(angle), rad*cos(angle) - bRadius, MP.boxsize[1]*sin(aa) );
        }
        glEnd();
      }

      //display borderline
      //xxmax and yymax are the maximum x an y extensions of the tilted caps
      xxmax = MP.boxsize[1]*cos(PI-MP.boxsize[2]);
      yymax = MP.boxsize[1]*sin(PI-MP.boxsize[2]);
      glBegin(GL_LINE_LOOP);
      //the bow in the front (+z)
      for( real aa = -MP.boxsize[2]; aa <=  MP.boxsize[2]; aa += 0.01 )
        glVertex3( -bRadius*sin(aa), bRadius*cos(aa)-bRadius, MP.boxsize[1] );
      //the left cap
      for( real aa = 0;  aa <= PI; aa += 0.01 ) {
        glVertex3( -center[0] + xxmax*sin(aa), center[1] - yymax*sin(aa), MP.boxsize[1]*cos(aa) );
      }
      //the bow in the back (-z)
      for( real aa = MP.boxsize[2]; aa >=  -MP.boxsize[2]; aa -= 0.01 )
        glVertex3( -bRadius*sin(aa), bRadius*cos(aa)-bRadius, -MP.boxsize[1] );
      //the right cap
      for( real aa = -PI;  aa <= 0; aa += 0.01 ) {
        glVertex3( center[0] + xxmax*sin(aa), center[1] + yymax*sin(aa), MP.boxsize[1]*cos(aa) );
      }
      glEnd();
    }
    break;
  }

  case SHAPE_ROUND_SQUARE:
    if ( DIM == 3 ) {
      displayBoxSection3( space, 0, MP.boxsize[3]-MP.boxsize[0], 0.1 );
      displayBoxSection3( space, 0, MP.boxsize[0]-MP.boxsize[3], 0.1 );
      displayBoxSection3( space, 1, MP.boxsize[3]-MP.boxsize[1], 0.1 );
      displayBoxSection3( space, 1, MP.boxsize[1]-MP.boxsize[3], 0.1 );
      displayBoxSection3( space, 2, MP.boxsize[2]-MP.boxsize[3], 0.1 );
      displayBoxSection3( space, 2, MP.boxsize[3]-MP.boxsize[2], 0.1 );
    } else {
      displayBoxSection3( space, 2, 0, 0.1 );
    }
    break;

  default:
    displayBoxSection3( space, 2, 0, 0.1 );
    if ( DIM == 3 ) {
      displayBoxSection3( space, 0, 0, 0.1 );
      displayBoxSection3( space, 1, 0, 0.1 );
    }
  break;
  }
}


//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
//some initialisation functions for OpenGL
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------


void Player::initGL3()
{
  gleClearColor( PP.bgcolor );

//  if ( PP.blend ) {
//  glEnable(GL_BLEND);
//  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
//  }

  // set lighting parameters

  //GLfloat matAmbDiff[]    = { 0.0, 0.5, 1.0, 1.0 };
  GLfloat matSpecular[]   = { 1.0, 1.0, 1.0, 1.0 };
  GLfloat matShininess[]  = { 50.0 };
  GLfloat light0Pos[]     = { 5.0, -3.0,  3.0, 0.0 };
  GLfloat light1Pos[]     = {-5.0,  3.0,  3.0, 0.0 };
  GLfloat light0Diffuse[] = { 1.0, 1.0, 1.0, 1.0 };
  GLfloat light0Spec[]    = { 1.0, 1.0, 1.0, 1.0 };
  GLfloat lModelAmbient[] = { 0.2, 0.2, 0.2, 1.0 };

  //glClearColor(0.0, 0.0, 0.0, 0.0);
  glShadeModel(GL_SMOOTH);

  //glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, matAmbDiff);
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
  glEnable(GL_LIGHT1);
  //glEnable(GL_RESCALE_NORMAL);
  glEnable(GL_NORMALIZE);

  glEnable(GL_POINT_SMOOTH);
  glEnable(GL_LINE_SMOOTH);
  glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);
  glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);

#if ( DIM == 3 )
  glEnable(GL_DEPTH_TEST);
  if ( PP.fog ) {
    glEnable(GL_FOG);
    glFogi(GL_FOG_MODE, GL_LINEAR);
    glFogf(GL_FOG_START, -cameraTranslation[2] -0.5*visibleRegion[2] );  // mid-front
    glFogf(GL_FOG_END,   -cameraTranslation[2] +    visibleRegion[2] );  // far back
    GLfloat fogColor[] = { 0.0, 0.0, 0.0, 1.0 };
    glFogfv(GL_FOG_COLOR, fogColor);
  }
#endif

  //openGL GLU Quadric initialisation
  qobj    = gluNewQuadric();
  extqobj = gluNewExtQuadric();
  gluQuadricDrawStyle(qobj, GLU_FILL);
  gluQuadricNormals(qobj, GLU_SMOOTH);
  gluExtQuadricDrawStyle(extqobj, GLU_FILL);
  gluExtQuadricNormals(extqobj, GLU_SMOOTH);
}

//CHAITANYA: display gradient of catastrophe and rescue in color using openGL commands

/*void Player::displayGradient( int gradTyp ){
	glColor3f( 1.0/2, 0.0, 0.0);
	int maxC=sim.resMap.nbCells();
	//printf("num. cells: %i\n", sim.resMap.nbCells() );
	real x_cor = -MP.boxsize[0], y_cor = -MP.boxsize[1];
	//for ( int n=0; n< maxC; n++ ){
		glRectf(x_cor, y_cor, (x_cor + MP.ficellsize), (y_cor + MP.ficellsize) ); //lowercorner_x, lowercorner_y, width, height
	//	x_cor = x_cor + MP.ficellsize;
	//	y_cor = y_cor + MP.ficellsize;
	//}
    glEnd();
	}
*/

//-----------------------------------------------------------------------------


void Player::drawFrame3()
{
  int ii;
  real ab;
  Vecteur xx, yy;
  real    tubeLength;
  //printf("disp %i\n", sim.frameInBuffer() );

  //----------------------------generate display lists-------------------------

  static bool initBox        = 0;
  static real oldCxSize      = 0.;
  static real oldMtSize      = 0.;
  static GLuint boxDispList  = glGenLists(1);
  static GLuint cxDispList   = glGenLists(1);
  static GLuint mtPtDispList = glGenLists(1);
  static GLuint rodDispList  = glGenLists(1);

  //create the display list for the bounding box
  if ( !initBox ) {
    glNewList(boxDispList, GL_COMPILE);
    if ( PP.boxwidth && Microtub::space ) {
      glLightningColorOI( PP.boxcolor, 0x00000000 );

      // if we are in 3D, check if we have a nice way to draw the box
      if ( DIM == 3 ) {
        switch( Microtub::space->getShape() ) {
          case SHAPE_SPHERE:
            gluSphere(qobj, MP.boxsize[0], 50, 50);
            break;
          case SHAPE_OVAL:
          case SHAPE_CYLINDER:
            oval3D(qobj, MP.boxsize[0], MP.boxsize[1], 5, 10);
            break;
          case SHAPE_BANANA:
            banana3D(extqobj, MP.boxsize[0], MP.boxsize[1], MP.boxsize[2], 5, 30);
            break;
          case SHAPE_TEE:
            tee3D(extqobj, qobj, MP.boxsize[0], MP.boxsize[1], MP.boxsize[2], MP.boxsize[3], MP.boxinflate, 3, 10);
            break;
          default:
            // no nice box, try it the normal way
            displayBox3( Microtub::space );
            break;
        }
      }
      // in 2D always draw the normal box
      else
      displayBox3( Microtub::space );

	  //--------------Display gradient
	  if ( PP.showGrad && MP.fieldTyp[0]){
	  	displayGradient( PP.showGrad );
	  }


      glEndList();
      initBox = 1;
    }
  }

  //if cxsize changed, regenerate complex display list
  if ( PP.cxsize != oldCxSize ) {
    //create the display list for the complexes
    glNewList(cxDispList, GL_COMPILE);
    glutSolidSphere(PP.cxsize/sizeFactor, 10, 10);
    glEndList();
    oldCxSize = PP.cxsize;
  }

  //if mtsize changed, regenerate mt point display list
  if ( PP.mtsize != oldMtSize ) {
    //create the display list for the complexes
    glNewList(mtPtDispList, GL_COMPILE);
    glutSolidSphere(PP.mtsize/sizeFactor, 10, 10);
    glEndList();
    oldMtSize = PP.mtsize;
  }


  //---------------------------------------------------------------------------
  //draw opaque objects first
  //depth buffer remains read/write
  //---------------------------------------------------------------------------


  //---------------------------- microtubules ----------------------
  ///\todo display arrows on microtubules

  Microtub* mt;
  glLightningColor( GL_FRONT_AND_BACK, PP.mtcolor[0] );

  //display lines only:
  if( PP.mtwidth && PP.mtlines && !PP.mtratchets && !PP.mtrainbow ) {
    for( mt = sim.firstMicrotub(); mt; mt = mt->next())
      if ( mt->mtColor >= 0) {
        glLightningColor( GL_FRONT_AND_BACK, PP.mtcolor[ mt->mtColor ] );
        //create a display list for the rods of this tube
        glNewList(rodDispList, GL_COMPILE);
        gluCylinder( qobj, PP.mtwidth/sizeFactor, PP.mtwidth/sizeFactor, mt->rodLength(), 10, 1 );
        glEndList();
        for( ii = 0; ii < mt->nbPoints()-1; ++ii ) {
          glPushMatrix();
          glRelocateTube( mt->pointAddr(ii), mt->pointAddr(ii+1) );
          glCallList(rodDispList);
          glPopMatrix();
        }
      }
  }

  //display segments which gradient of color indicating polarity
  if( PP.mtwidth && PP.mtratchets && !PP.mtrainbow ) {
    for( mt = sim.firstMicrotub(); mt; mt = mt->next())
      if ( mt->mtColor >= 0) {
        //create a display list for the rods
        glNewList(rodDispList, GL_COMPILE);
        gluEnhancedCylinder( extqobj, PP.mtwidth/sizeFactor, PP.mtwidth/sizeFactor, mt->rodLength(), 10, 1,
                             PP.bgcolor, PP.mtcolor[ mt->mtColor ], FRONT_AND_BACK );
        glEndList();
        for( ii = 0; ii < mt->nbSegments(); ++ii ) {
          //calculate the rotation angle for the cylinders
          glPushMatrix();
          glRelocateTube( mt->pointAddr(ii), mt->pointAddr(ii+1) );
          glCallList(rodDispList);
          glPopMatrix();
        }
      }
  }

  //display fix points of the tubes (every micro-meter):
  if( PP.mtspeckles ) {
    for( mt = sim.firstMicrotub(); mt; mt = mt->next())
      if ( mt->mtColor >= 0) {
        glLightningColor( GL_FRONT_AND_BACK, PP.mtcolor[ mt->mtColor ] );
        ab = PP.mtspeckledist * ceil( mt->abscissa(MT_MINUS_END) / PP.mtspeckledist );
        while( ab <= mt->abscissa(MT_PLUS_END) ) {
          glPushMatrix();
          glTranslatev( mt->where(ab, MT_ORIGIN) );
          glCallList(mtPtDispList);
          glPopMatrix();
          ab += PP.mtspeckledist;
        }
      }
  }

  //display simulated points of the tubes:
  if( PP.mtpoints ) {

    for( mt = sim.firstMicrotub(); mt; mt = mt->next())
      if ( mt->mtColor >= 0) {
        glLightningColor( GL_FRONT_AND_BACK, PP.mtcolor[ mt->mtColor ] );
        for( ii = 0; ii < mt->nbPoints(); ++ii ) {
          glPushMatrix();
          glTranslatev( mt->pointAddr(ii) );
          glCallList(mtPtDispList);
          glPopMatrix();
        }
      }
  }

  //display points at minus-ends:
  if( PP.mtends[MT_MINUS_END] ) {
    for( mt = sim.firstMicrotub(); mt; mt = mt->next())
      if ( mt->mtColor >= 0) {
        glPushMatrix();
        glLightningColor( GL_FRONT_AND_BACK, PP.mtcolor[ mt->getDynamicState(MT_MINUS_END) ] );
        glTranslatev( mt->pointAddr(0) );
        glutSolidSphere(2.*PP.mtsize/sizeFactor, 10, 10);
        glPopMatrix();
      }
  }

  //display points at plus-ends:
  if( PP.mtends[MT_PLUS_END] ) {
    for( mt = sim.firstMicrotub(); mt; mt = mt->next())
      if ( mt->mtColor >= 0) {
        glPushMatrix();
        glLightningColor( GL_FRONT_AND_BACK, PP.mtcolor[ mt->getDynamicState(MT_PLUS_END) ] );
        glTranslatev( mt->pointAddr(mt->nbSegments()) );
        glutSolidSphere(2.*PP.mtsize/sizeFactor, 10, 10);
        glPopMatrix();
      }
  }

  //display forces acting on the tubes:
  if( PP.mtwidth && PP.mtforces ) {
    glLightningColor( GL_FRONT_AND_BACK, PP.mtfcolor );
    for( mt = sim.firstMicrotub(); mt; mt = mt->next())
      if (( mt->mtColor >= 0 ) && ( mt->getForcesAddr() )) {
        for( ii = 0; ii < mt->nbPoints(); ++ii ) {
          xx = mt->whereP(ii) + PP.mtforces * mt->forcesP(ii);
          glPushMatrix();
          tubeLength = glRelocateTube( mt->pointAddr(ii), xx );
          gluCylinder( qobj, PP.mtwidth/sizeFactor, PP.mtwidth/sizeFactor, tubeLength, 10, 1 );
          glPopMatrix();
        }
      }
  }

  //display projection of the points on the box:
  if ( PP.mtwidth && PP.mtproject ) {
    glLightningColor( GL_FRONT_AND_BACK, PP.mtfcolor );
    for( mt = sim.firstMicrotub(); mt; mt = mt->next())
      if ( mt->mtColor >= 0 )
        for( ii = 0; ii < mt->nbPoints(); ++ii ) {
          xx = mt->whereP( ii );
          if ( Microtub::space->isOutside( xx ) ) {
            Microtub::space->project( xx, yy );
            glPushMatrix();
            tubeLength = glRelocateTube( xx, yy );
            gluCylinder( qobj, PP.mtwidth/sizeFactor, PP.mtwidth/sizeFactor, tubeLength, 10, 1 );
            glPopMatrix();

            glPushMatrix();
            glTranslatev( yy );
            glutSolidSphere( PP.mtsize/sizeFactor, 10, 10 );
            glPopMatrix();
          }
        }
  }

  //display the internal tensions of the microtubules:
  if ( PP.width && PP.mtrainbow ) {
    for( mt = sim.firstMicrotub(); mt; mt = mt->next())
      if ( mt->mtColor >= 0 ) {
        for( ii = 0; ii < mt->nbSegments(); ++ii ) {
          glLightningColorRainbow( 0.2 - mt->getLagrange(ii) * PP.mtrainbow );
          glPushMatrix();
          glRelocateTube( mt->pointAddr(ii), mt->pointAddr(ii+1) );
          gluCylinder( qobj, PP.mtwidth/sizeFactor, PP.mtwidth/sizeFactor, mt->rodLength(), 10, 1 );
          glPopMatrix();
        }
      }
  }


  //---------------------------- solids ----------------------------

  //display points:
  if ( PP.sopoints ) {
    glPointSize(PP.sosize);
    glBegin(GL_POINTS);
    glColor(PP.socolor);
    for(Solid * so=sim.firstSolid(); so ; so=so->next()) {
      for( ii = 0; ii < so->nbPoints(); ++ii )
        glVertex( so->pointAddr( ii ) );
    }
    glEnd();
  }

  //display a close loop joinning all the points
  if ( PP.solinks ) {
    glLineWidth( PP.sosize );
    glColor(PP.socolor);
    for(Solid * so=sim.firstSolid(); so ; so=so->next()) {
      glBegin(GL_LINE_LOOP);
      for( ii = 0; ii < so->nbPoints(); ++ii )
        glVertex( so->pointAddr( ii ) );
      glEnd();
    }
  }


  //---------------------------- nuclei ----------------------------

  //display the nuclear envelope
  if ( PP.nuenvelope ) {
    glLineWidth( PP.nuenvwidth );
    glLightningColorOI( PP.nuenvcolor, 0x00000000 );
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    for(Nucleus * nu=sim.firstNucleus(); nu ; nu=nu->next()) {
      glPushMatrix();
      if(DIM == 2) {
        glTranslate(nu->coord(0), nu->coord(1), 0.0);
        glutSolidTorus( lineWidth, (GLdouble)nu->radius(), 10, 50 );
      }
      else if(DIM == 3) {
        glTranslate(nu->coord(0), nu->coord(1), nu->coord(2));
        // rotate the nucleus
        int refpoint1 = DIM*(nu->nbPoints()-3);
        int refpoint2 = DIM*(nu->nbPoints()-2);
        int refpoint3 = DIM*(nu->nbPoints()-1);
        static GLdouble sphereRotMat[16];

        for( int ii = 0; ii < 3; ii++) {
          sphereRotMat[ii]     = (nu->coord(refpoint1+ii) - nu->coord(ii))/nu->radius();
          sphereRotMat[ii+4]   = (nu->coord(refpoint2+ii) - nu->coord(ii))/nu->radius();
          sphereRotMat[ii+8]   = (nu->coord(refpoint3+ii) - nu->coord(ii))/nu->radius();
          sphereRotMat[ii+12]  = 0.;
          sphereRotMat[ii*4+3] = 0.;
        }
        sphereRotMat[15] = 1.;
        glMultMatrixd( sphereRotMat );

        glutWireSphere(nu->radius(), 16, 16);
      }
      glPopMatrix();
    }
    glDisable(GL_BLEND);
  }

  //display the points on the nucleus (without the center-point)
  if ( PP.nupoints ) {
    glLightningColor( GL_FRONT_AND_BACK, PP.nuptcolor );
    for(Nucleus * nu=sim.firstNucleus(); nu ; nu=nu->next())
      for( ii = 1; ii < nu->nbPoints()-nu->nbrefpts; ++ii ) {
        glPushMatrix();
        glTranslatev( nu->pointAddr( ii ) );
        glutSolidSphere(PP.nuptsize/sizeFactor, 10, 10);
        glPopMatrix();
      }
  }

  //display the reference points
  if ( PP.nurefpts ) {
    glLightningColor( GL_FRONT_AND_BACK, PP.nurefcolor );
    for(Nucleus * nu=sim.firstNucleus(); nu ; nu=nu->next())
      for( ii = nu->nbPoints()-nu->nbrefpts; ii < nu->nbPoints(); ++ii ) {
        glPushMatrix();
        glTranslatev( nu->pointAddr( ii ) );
        glutSolidSphere(PP.nuptsize/sizeFactor, 10, 10);
        glPopMatrix();
      }
  }

  //display the mtoc-microtubule links
  if ( PP.nulinks ) {
    PointExact pte1, pte2;
    //glLightningColor( GL_FRONT_AND_BACK, PP.nulinkcolor );
    for(Nucleus * nu=sim.firstNucleus(); nu ; nu=nu->next())
      for( int nuindx = 0 ; nuindx < nu->maxNbMicrotub(); nuindx++ ) {
	if ( nu -> setPointClamp( &pte1, &pte2, nuindx )) {
	  glPushMatrix();
	  tubeLength = glRelocateTube( pte1.where(), pte2.where() );
	  //gluCylinder(qobj, PP.mtwidth/sizeFactor, PP.mtwidth/sizeFactor, tubeLength, 10, 1);
	  gluEnhancedCylinder(extqobj, PP.mtwidth/sizeFactor, PP.mtwidth/sizeFactor, tubeLength, 10, 1,
			      PP.nuptcolor, PP.nulinkcolor, FRONT_AND_BACK );
	  glPopMatrix();
	}
      }
  }

  //--------------------------- grafted ----------------------------

  Grafted * gh;

  // display the grafted positions
  if ( PP.ghhands && ( PP.ghselect & 2 )) {
    for(gh=sim.firstBoundGrafted(); gh ; gh=gh->next()) {
      if ( PP.ghrainbow )
        glLightningColorRainbow( gh->getStretch().norm() * PP.ghrainbow );
      else
        glLightningColor( GL_FRONT_AND_BACK, PP.hacolor[ gh->getType() ] );
      glPushMatrix();
      glTranslatev( gh->whereHand() );
      glutSolidSphere( PP.ghsize/sizeFactor, 10, 10 );
      glPopMatrix();
    }
  }


  //------------------------- complexes ---------------------------


  Complex * cx;


  // display bridging complexes
  if ( PP.cxhands && ( PP.cxselect & 4 )) {
    for(cx=sim.firstBridgeComplex(); cx ; cx=cx->next()) {
      glPushMatrix();
      glTranslatev( cx->where1() );
      if ( PP.cxrainbow ) {
        glLightningColorRainbow( cx->getStretch().norm() * PP.cxrainbow );
      } else {
        glLightningColor( GL_FRONT_AND_BACK, PP.hacolor[ cx->getType1() ] );
      }
      glCallList(cxDispList);
      glPopMatrix();

      glPushMatrix();
      glTranslatev( cx->where2() );
      if ( PP.cxrainbow ) {
        glLightningColorRainbow( cx->getStretch().norm() * PP.cxrainbow );
      } else {
        glLightningColor( GL_FRONT_AND_BACK, PP.hacolor[ cx->getType2() ] );
      }
      glCallList(cxDispList);
      glPopMatrix();
    }
  }

  // display the link for bridging complexes
  if ( PP.cxlinks && ( PP.cxselect & 4 )) {
    for(cx=sim.firstBridgeComplex(); cx ; cx=cx->next()) {
      xx = cx->where1();
      yy = cx->where2();
      Space::moduloNear( yy, xx );

      glPushMatrix();
      tubeLength = glRelocateTube( xx, yy );
      if ( PP.cxrainbow ) {
        glLightningColorRainbow( cx->getStretch().norm() * PP.cxrainbow );
        gluCylinder(qobj, PP.cxwidth/sizeFactor, PP.cxwidth/sizeFactor, tubeLength, 10, 1);
      } else {
        gluEnhancedCylinder(extqobj, PP.cxwidth/sizeFactor, PP.cxwidth/sizeFactor, tubeLength, 10, 1,
                            PP.hacolor[ cx->getType1() ], PP.hacolor[ cx->getType2() ], FRONT_AND_BACK );
      }
      glPopMatrix();
    }
  }


  //--------------------------------------------------------------------------
  //draw translucent objects
  //depth buffer is now read only, so order of drawing is important here!
  //--------------------------------------------------------------------------


  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glDepthMask(GL_FALSE);


  //--------------------------- grafted ----------------------------

  //display the attached position:
  if ( PP.ghhands && ( PP.ghselect & 1 )) {
    glLightningColor( GL_BACK, 0x00000000 );
    for(gh=sim.firstFreeGrafted(); gh ; gh=gh->next()) {
      glLightningColor( GL_FRONT_AND_BACK, PP.hacolor[ gh->getType() ], 2 );
      glPushMatrix();
      glTranslatev( gh->whereGrafted() );
      glutSolidSphere( PP.ghsize/sizeFactor, 10, 10 );
      glPopMatrix();
    }
  }

  // display the grafted link or forces
  if (( PP.ghlinks || PP.ghrainbow ) && ( PP.ghselect & 2 )) {
    // make the inside of the tubes invisible
    glLightningColor( GL_BACK, 0x00000000 );
    for(gh=sim.firstBoundGrafted(); gh ; gh=gh->next()) {
      xx = gh->whereGrafted();
      yy = gh->whereHand();
      Space::moduloNear( yy, xx );

      glPushMatrix();
      tubeLength = glRelocateTube( xx, yy );
      if ( PP.ghrainbow ) {
        glLightningColorRainbow( gh->getStretch().norm() * PP.ghrainbow );
        gluCylinder( qobj, PP.ghwidth/sizeFactor, PP.ghwidth/sizeFactor, tubeLength, 10, 1 );
      } else {
        // the tubes should have a gradient in the alpha value
        // devide given alpha value by 4:
        long int alpha     = (long int)(((unsigned char)(PP.hacolor[ gh->getType() ] & 0x000000FF)) >> 2);
        // set new alpha value:
        long int tubeColor = ( (PP.hacolor[ gh->getType() ] & 0xFFFFFF00) + alpha );
        //printf("tubeColor: %x, hacolor: %x\n", tubeColor, PP.hacolor[ gh->getType() ]);
        gluEnhancedCylinder( extqobj, PP.ghwidth/sizeFactor, PP.ghwidth/sizeFactor, tubeLength, 10, 1,
                             tubeColor, PP.hacolor[ gh->getType() ] );
      }
      glPopMatrix();
    }
  }

  // display the grafted displacements
  if ( PP.ghdisp && ( PP.ghselect & 1 )) {
    glLightningColor( GL_BACK, 0x00000000 );
    //give the backgroundcolor a zero alpha value to make the tubes vanish
    long int bgcolor = ((PP.bgcolor | 0x000000FF) & 0xFFFFFF00);

    for(gh=sim.firstFreeGrafted(); gh ; gh=gh->next()) {
      glPushMatrix();
      tubeLength = glRelocateTube( gh->getPosition(), gh->getPositionOld() );
      gluEnhancedCylinder( extqobj, PP.ghwidth/sizeFactor, PP.ghwidth/sizeFactor, tubeLength, 10, 1,
                           PP.hacolor[ gh->getType() ], bgcolor );
      glPopMatrix();
      gh->setPositionOld();
    }
  }

  if ( PP.ghdisp && ( PP.ghselect & 2 )) {
    glLightningColor( GL_BACK, 0x00000000 );
    //give the backgroundcolor a zero alpha value to make the tubes vanish
    long int bgcolor = ((PP.bgcolor | 0x000000FF) & 0xFFFFFF00);

    for(gh=sim.firstBoundGrafted(); gh ; gh=gh->next()) {
      glPushMatrix();
      tubeLength = glRelocateTube( gh->getPosition(), gh->getPositionOld() );
      gluEnhancedCylinder( extqobj, PP.ghwidth/sizeFactor, PP.ghwidth/sizeFactor, tubeLength, 10, 1,
                           PP.hacolor[ gh->getType() ], bgcolor );
      glPopMatrix();
      gh->setPositionOld();
    }
  }


  //------------------------- complexes ---------------------------

  // display free complexes
  if ( PP.cxhands && ( PP.cxselect & 1 )) {
    glLightningColor( GL_BACK, 0x00000000 );
    for(cx=sim.firstFreeComplex(); cx ; cx=cx->next()) {
      glLightningColor( GL_FRONT, PP.hacolor[ cx->getType1() ], 2 );
//      glLightningColor( GL_FRONT_AND_BACK, PP.hacolor[ cx->getType1() ]);
      glPushMatrix();
      glTranslatev( cx->getPosition() );
      glCallList(cxDispList);
      glPopMatrix();
      glLightningColor( GL_FRONT, PP.hacolor[ cx->getType2() ], 2 );
      glPushMatrix();
      glTranslatev( cx->getPosition() );
      glCallList(cxDispList);
      glPopMatrix();
    }
  }


  // display bound complexes
  if ( PP.cxhands && ( PP.cxselect & 2 )) {
    glLightningColor( GL_BACK, 0x00000000 );
    for(cx=sim.firstBoundComplex(); cx ; cx=cx->next()) {
      if ( cx->isAttached1() )
        glLightningColor( GL_FRONT, PP.hacolor[ cx->getType1() ], 1 );
      else
        glLightningColor( GL_FRONT, PP.hacolor[ cx->getType2() ], 1 );
      glPushMatrix();
      glTranslatev( cx->calculatePosition() );
      glCallList(cxDispList);
      glPopMatrix();
    }
  }


  // display the displacements of the complexes
  if ( PP.cxdisp ) {
    glLightningColor( GL_BACK, 0x00000000 );
    //give the backgroundcolor a zero alpha value to make the tubes vanish
    long int bgcolor = ((PP.bgcolor | 0x000000FF) & 0xFFFFFF00);

    for(cx=sim.firstFreeComplex(); cx ; cx=cx->next()) {
      if ( PP.cxselect & 1 ) {
        glPushMatrix();
        tubeLength = glRelocateTube( cx->getPosition(), cx->getPositionOld() );
        gluEnhancedCylinder( extqobj, PP.cxwidth/sizeFactor, PP.cxwidth/sizeFactor, tubeLength, 10, 1,
                             PP.hacolor[ cx->getType1() ], bgcolor );
        glPopMatrix();
      }
      cx->setPositionOld();
    }


    for(cx=sim.firstBoundComplex(); cx ; cx=cx->next()) {
      cx->getPosition();
      if ( PP.cxselect & 2 ) {
        glPushMatrix();
        tubeLength = glRelocateTube( cx->getPosition(), cx->getPositionOld() );
        gluEnhancedCylinder( extqobj, PP.cxwidth/sizeFactor, PP.cxwidth/sizeFactor, tubeLength, 10, 1,
                             PP.hacolor[ cx->getType1() ], bgcolor );
        glPopMatrix();
      }
      cx->setPositionOld();
    }

    for(cx=sim.firstBridgeComplex(); cx ; cx=cx->next()) {
      cx->getPosition();
      if ( PP.cxselect & 4 ) {
        glPushMatrix();
        tubeLength = glRelocateTube( cx->getPosition(), cx->getPositionOld() );
        gluEnhancedCylinder( extqobj, PP.cxwidth/sizeFactor, PP.cxwidth/sizeFactor, tubeLength, 10, 1,
                             PP.hacolor[ cx->getType1() ], bgcolor );
        glPopMatrix();
      }
      cx->setPositionOld();
    }
  }


  //----------------- draw a frame around the box:
  glCallList( boxDispList );


  glDepthMask(GL_TRUE);
  glDisable(GL_BLEND);


  //set color to white to display help text
  glLightningColor( GL_FRONT_AND_BACK, 0xFFFFFFFF );
}
