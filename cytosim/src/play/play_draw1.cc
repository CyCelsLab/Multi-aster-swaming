//RCS: $Id: play_draw1.cc,v 2.30 2005/04/22 18:53:25 nedelec Exp $

// drawFrame1() implements a wireframe display using Lines and Points.
// This works well if anti-aliasing is enabled (PC).

// Nk: 31 - 03 -2016, definining a color scheme for complexes - different from the hand color scheme
// Neha Khetan, 10 July 2018 - Modifying the color of links of the complex so that it is same as the color of complex - Changed hacolor to cxcolor

void Player::initGL1()
{
  gleClearColor( PP.bgcolor );
  // blending is enabled by default in styles 1 and 2
  if ( PP.blend ) {
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glEnable(GL_POINT_SMOOTH);
    glEnable(GL_LINE_SMOOTH);
    glHint(GL_POINT_SMOOTH_HINT, GL_FASTEST);
    glHint(GL_LINE_SMOOTH_HINT, GL_FASTEST);
  }

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
}


//---------------------------------------------------------------------------
void displayBoxSection1(const Space * mSpace, int dim, float pos, float step)
{
  Vecteur q, p( pos, pos, pos );
  int xx = ( dim + 1 ) % DIM;
  int yy = ( xx + 1 ) % DIM;
  real xs = mSpace->getBoundingRect()[xx];
  real ys = mSpace->getBoundingRect()[yy];

  glBegin(GL_LINE_LOOP);
  p[yy] = ys;
  for( real a = -xs; a < xs; a += step ) {
    p[xx] = a;
    mSpace->project(p, q);
    glVertex( q );
  };
  p[xx] = xs;
  for( real a = -ys; a < ys; a += step ) {
    p[yy] = -a;
    mSpace->project(p, q);
    glVertex( q );
  };
  p[yy] = -ys;
  for( real a = -xs; a < xs; a += step ) {
    p[xx] = -a;
    mSpace->project(p, q);
    glVertex( q );
  };
  p[xx] = -xs;
  for( real a = -ys; a < ys; a += step ) {
    p[yy] = a;
    mSpace->project(p, q);
    glVertex( q );
  };
  glEnd();

}

//we define glVertexPseudo(ARG) to avoid using the variable mt->pseudoY
//(and being obliged to allocate it) if PSEUDO_YZ is not enabled
#ifdef PSEUDO_YZ
  #define glVertexPseudo(ARG) glVertex(ARG, mt->pseudoY);
#else
  #define glVertexPseudo(ARG) glVertex(ARG);
#endif
//---------------------------------------------------------------------------
void displayBox1(const Space * space)
{
  glColor( PP.boxcolor );
  glLineWidth( PP.boxwidth );

  switch( space->getShape() ) {

  case SHAPE_PERIODIC:
  case SHAPE_SQUARE:
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
    displayBoxSection1( space, 2, 0, 0.1 );
    if ( DIM == 3 ) {
      for( real Z = +0; Z <= +MP.boxsize[0]; Z += 1.0 )
        displayBoxSection1( space, 0, Z, 0.1 );
      for( real Z = -1; Z >= -MP.boxsize[0]; Z -= 1.0 )
        displayBoxSection1( space, 0, Z, 0.1 );
      displayBoxSection1( space, 1, 0, 0.1 );
    }
    break;

  case SHAPE_BANANA: {
    real center[2];
    real bRadius = MP.boxsize[0]/MP.boxsize[2];

    //the rightmost point on the backbone of the banana
    //(for the leftmost point we use -center[0] below)
    center[0] = bRadius*sin(MP.boxsize[2]);
    center[1] = bRadius*cos(MP.boxsize[2]) - bRadius;

    glColor( PP.boxcolor );
    glLineWidth( PP.boxwidth );

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

    #if ( DIM == 3 )
      real angle, rad;
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
    #endif

    break;
  }

  case SHAPE_TEE: {

    const int  pisteps = 50; //how many lines to draw for a circle of pi
                             //this should be a multiple of 2 (we devide by 2 below)
    const real pifrac  = PI/(real)pisteps;
    const real inflRad = MP.boxsize[1]+MP.boxinflate;
    const real sqrttwo = sqrt(2.);

    glColor( PP.boxcolor );
    glLineWidth( PP.boxwidth );

    glBegin(GL_LINE_LOOP);
    //the upper side from the tJunction to the left
    glVertex3(  MP.boxsize[2]-MP.boxsize[1], inflRad, 0. );
    glVertex3( -MP.boxsize[0],               inflRad, 0. );
    //the left cap
    for( int aa = 1; aa < pisteps; aa++ )
      glVertex3( -MP.boxsize[0]-inflRad*sin(pifrac*aa), inflRad*cos(pifrac*aa), 0. );
    //the lower side from left to right
    glVertex3( -MP.boxsize[0], -inflRad, 0. );
    glVertex3(  MP.boxsize[0], -inflRad, 0. );
    //the right cap
    for( int aa = 1; aa < pisteps; aa++ )
      glVertex3( MP.boxsize[0]+inflRad*sin(pifrac*aa), -inflRad*cos(pifrac*aa), 0. );
    //the upper side from the right to the tJunction
    glVertex3( MP.boxsize[0],               inflRad, 0. );
    glVertex3( MP.boxsize[2]+MP.boxsize[1], inflRad, 0. );
    //the right bow if the tee is inflated
    if( MP.boxinflate ) {
      for( int aa = 1; aa <= pisteps/2; aa++ )
        glVertex3( MP.boxsize[2]+MP.boxsize[1]+MP.boxinflate*sin(pifrac*aa), \
                   MP.boxsize[1]+MP.boxinflate*cos(pifrac*aa), 0. );
    }
    //the right side of the arm
    glVertex3( MP.boxsize[2]+inflRad, MP.boxsize[1]+MP.boxsize[3], 0. );
    //the cap of the arm
    for( int aa = 1; aa < pisteps; aa++ )
      glVertex3( MP.boxsize[2]+inflRad*cos(pifrac*aa), MP.boxsize[3]+MP.boxsize[1]+inflRad*sin(pifrac*aa), 0. );
    //the left side of the arm
    glVertex3( MP.boxsize[2]-inflRad, MP.boxsize[1]+MP.boxsize[3], 0. );
    //the left bow if the tee is inflated
    if( MP.boxinflate ) {
      for( int aa = 1; aa < pisteps/2; aa++ )
        glVertex3( MP.boxsize[2]-MP.boxsize[1]-MP.boxinflate*cos(pifrac*aa), \
                   MP.boxsize[1]+MP.boxinflate*sin(pifrac*aa), 0. );
    }
    glEnd();

    #if ( DIM == 3 )

    //display a line at 45 degrees
    glBegin(GL_LINE_LOOP);
    //the base cylinder in the lower back
    glVertex3(  MP.boxsize[0], -inflRad/sqrttwo, -inflRad/sqrttwo );
    glVertex3( -MP.boxsize[0], -inflRad/sqrttwo, -inflRad/sqrttwo );
    //the left cap
    for( int aa = 1; aa < pisteps; aa++ )
      glVertex3( -MP.boxsize[0]-inflRad*sin(pifrac*aa), \
                 -inflRad*cos(pifrac*aa)/sqrttwo, \
                 -inflRad*cos(pifrac*aa)/sqrttwo );
    //the base cylinder in the upper front, left
    glVertex3( -MP.boxsize[0],                       inflRad/sqrttwo, inflRad/sqrttwo );
    glVertex3(  MP.boxsize[2]-MP.boxsize[1]/sqrttwo, inflRad/sqrttwo, inflRad/sqrttwo );
    //the left bow in front if the tee is inflated
    if( MP.boxinflate ) {
      for( int gg = 1; gg <= pisteps/2; gg++ )
        glVertex3( MP.boxsize[2]-(MP.boxsize[1]+MP.boxinflate*sin(pifrac*gg))/sqrttwo, \
                   (MP.boxsize[1]+MP.boxinflate*cos(pifrac*gg))/sqrttwo, \
                   inflRad/sqrttwo );
    }
    //the arm front, left
//    glVertex3(  MP.boxsize[2]-inflRad/sqrttwo, MP.boxsize[1]/sqrttwo,       inflRad/sqrttwo );
    glVertex3(  MP.boxsize[2]-inflRad/sqrttwo, MP.boxsize[1]+MP.boxsize[3], inflRad/sqrttwo );
    //the cap of the arm
    for( int aa = 1; aa < pisteps; aa++ )
      glVertex3( MP.boxsize[2]-inflRad*cos(pifrac*aa)/sqrttwo, \
                 MP.boxsize[1]+MP.boxsize[3]+inflRad*sin(pifrac*aa), \
                 inflRad*cos(pifrac*aa)/sqrttwo );
    //the arm back, right
    glVertex3(  MP.boxsize[2]+inflRad/sqrttwo, MP.boxsize[1]+MP.boxsize[3], -inflRad/sqrttwo );
    glVertex3(  MP.boxsize[2]+inflRad/sqrttwo, MP.boxsize[1]/sqrttwo,       -inflRad/sqrttwo );
    //the right bow in the back if the tee is inflated
    if( MP.boxinflate ) {
      for( int gg = 1; gg <= pisteps/2; gg++ )
        glVertex3( MP.boxsize[2]+(MP.boxsize[1]+MP.boxinflate*cos(pifrac*gg))/sqrttwo, \
                   (MP.boxsize[1]+MP.boxinflate*sin(pifrac*gg))/sqrttwo, \
                   -inflRad/sqrttwo );
    }
    //the base cylinder in the upper back, right
//    glVertex3( MP.boxsize[2]+MP.boxsize[1]/sqrttwo, inflRad/sqrttwo, -inflRad/sqrttwo );
    glVertex3( MP.boxsize[0],                       inflRad/sqrttwo, -inflRad/sqrttwo );
    //the right cap
    for( int aa = 1; aa < pisteps; aa++ )
      glVertex3(  MP.boxsize[0]+inflRad*sin(pifrac*aa), \
                  inflRad*cos(pifrac*aa)/sqrttwo, \
                 -inflRad*cos(pifrac*aa)/sqrttwo );
    //the base cylinder in the lower front
    glVertex3(  MP.boxsize[0], -inflRad/sqrttwo, inflRad/sqrttwo );
    glVertex3( -MP.boxsize[0], -inflRad/sqrttwo, inflRad/sqrttwo );
    //the left cap
    for( int aa = 1; aa < pisteps; aa++ )
      glVertex3( -MP.boxsize[0]-inflRad*sin(pifrac*aa), \
                 -inflRad*cos(pifrac*aa)/sqrttwo, \
                  inflRad*cos(pifrac*aa)/sqrttwo );
    //the base cylinder in the upper back, left
    glVertex3( -MP.boxsize[0],                       inflRad/sqrttwo, -inflRad/sqrttwo );
    glVertex3(  MP.boxsize[2]-MP.boxsize[1]/sqrttwo, inflRad/sqrttwo, -inflRad/sqrttwo );
    //the left bow in the back if the tee is inflated
    if( MP.boxinflate ) {
      for( int gg = 1; gg <= pisteps/2; gg++ )
        glVertex3( MP.boxsize[2]-(MP.boxsize[1]+MP.boxinflate*sin(pifrac*gg))/sqrttwo, \
                   (MP.boxsize[1]+MP.boxinflate*cos(pifrac*gg))/sqrttwo, \
                   -inflRad/sqrttwo );
    }
    //the arm back, left
//    glVertex3(  MP.boxsize[2]-inflRad/sqrttwo, MP.boxsize[1]/sqrttwo,       -inflRad/sqrttwo );
    glVertex3(  MP.boxsize[2]-inflRad/sqrttwo, MP.boxsize[1]+MP.boxsize[3], -inflRad/sqrttwo );
    //the cap of the arm
    for( int aa = 1; aa < pisteps; aa++ )
      glVertex3(  MP.boxsize[2]-inflRad*cos(pifrac*aa)/sqrttwo, \
                  MP.boxsize[1]+MP.boxsize[3]+inflRad*sin(pifrac*aa), \
                 -inflRad*cos(pifrac*aa)/sqrttwo );
    //the arm front, right
    glVertex3(  MP.boxsize[2]+inflRad/sqrttwo, MP.boxsize[1]+MP.boxsize[3], inflRad/sqrttwo );
    glVertex3(  MP.boxsize[2]+inflRad/sqrttwo, MP.boxsize[1]/sqrttwo,       inflRad/sqrttwo );
    //the right bow in the front if the tee is inflated
    if( MP.boxinflate ) {
      for( int gg = 1; gg <= pisteps/2; gg++ )
        glVertex3( MP.boxsize[2]+(MP.boxsize[1]+MP.boxinflate*cos(pifrac*gg))/sqrttwo, \
                   (MP.boxsize[1]+MP.boxinflate*sin(pifrac*gg))/sqrttwo, \
                   inflRad/sqrttwo );
    }
    //the base cylinder in the upper front, right
//    glVertex3( MP.boxsize[2]+MP.boxsize[1]/sqrttwo, inflRad/sqrttwo, inflRad/sqrttwo );
    glVertex3( MP.boxsize[0],                       inflRad/sqrttwo, inflRad/sqrttwo );
    //the right cap
    for( int aa = 1; aa < pisteps; aa++ )
      glVertex3(  MP.boxsize[0]+inflRad*sin(pifrac*aa), \
                  inflRad*cos(pifrac*aa)/sqrttwo, \
                  inflRad*cos(pifrac*aa)/sqrttwo );
    glEnd();

    //display the borderline of the base cylinder
    glBegin(GL_LINE_LOOP);
    //the base cylinder in the back (-z)
    glVertex3(  MP.boxsize[0], 0., -inflRad );
    glVertex3( -MP.boxsize[0], 0., -inflRad );
    //the left cap
    for( int aa = 1; aa < pisteps; aa++ )
      glVertex3( -MP.boxsize[0]-inflRad*sin(pifrac*aa), 0., -inflRad*cos(pifrac*aa) );
    //the base cylinder in the front (+z)
    glVertex3( -MP.boxsize[0], 0., inflRad );
    glVertex3(  MP.boxsize[0], 0., inflRad );
    //the right cap
    for( int aa = 1; aa < pisteps; aa++ )
      glVertex3(  MP.boxsize[0]+inflRad*sin(pifrac*aa), 0.,  inflRad*cos(pifrac*aa) );
    glEnd();

    //display the borderline of the arm
    glBegin(GL_LINE_LOOP);
    //the cylinder in the back (-z)
    glVertex3( MP.boxsize[2], 0                          , -inflRad );
    glVertex3( MP.boxsize[2], MP.boxsize[1]+MP.boxsize[3], -inflRad );
    //the upper cap
    for( int aa = 1; aa < pisteps; aa++ )
      glVertex3( MP.boxsize[2], MP.boxsize[1]+MP.boxsize[3]+inflRad*sin(pifrac*aa), -inflRad*cos(pifrac*aa) );
    //the cylinder in the front (+z)
    glVertex3( MP.boxsize[2], MP.boxsize[1]+MP.boxsize[3], inflRad );
    glVertex3( MP.boxsize[2], 0                          , inflRad );
    //the lower cap
    for( int aa = 1; aa < pisteps; aa++ )
      glVertex3( MP.boxsize[2], -inflRad*sin(pifrac*aa), inflRad*cos(pifrac*aa) );
    glEnd();

    //display the intersection area
    glBegin(GL_LINE_LOOP);
    for( int aa = 0; aa < 2*pisteps; aa++ )
      glVertex3( MP.boxsize[2]-(MP.boxsize[1]+MP.boxinflate/sqrttwo)*cos(pifrac*aa), \
                 (MP.boxsize[1]+MP.boxinflate/sqrttwo)*fabs(cos(pifrac*aa)), \
                 inflRad*sin(pifrac*aa) );
    glEnd();
    //the edges of the curvature if the tee is inflated
    if( MP.boxinflate ) {
      //the left edge
      glBegin(GL_LINE_LOOP);
      for( int aa = 0; aa < 2*pisteps; aa++ )
        glVertex3( MP.boxsize[2]-MP.boxsize[1]*cos(pifrac*aa), \
                   inflRad*fabs(cos(pifrac*aa)), \
                   inflRad*sin(pifrac*aa) );
      glEnd();
      //the right edge
      glBegin(GL_LINE_LOOP);
      for( int aa = 0; aa < 2*pisteps; aa++ )
        glVertex3( MP.boxsize[2]-inflRad*cos(pifrac*aa), \
                   MP.boxsize[1]*fabs(cos(pifrac*aa)), \
                   inflRad*sin(pifrac*aa) );
      glEnd();
  }

    //display rings ~ every 1.2 um
    //the distance between the rings is different for the cylinders left and
    //right of the arm and for the arm itself. Maybe we should make the
    //distances equal.
    real ringPos;
    int nbSegmentsL  = (int)floor((MP.boxsize[2]-MP.boxsize[1] + MP.boxsize[0]) / 1.2 );
    int nbSegmentsR  = (int)floor((MP.boxsize[0] - (MP.boxsize[2]+MP.boxsize[1])) / 1.2 );
    int nbSegmentsUp = (int)floor(MP.boxsize[3] / 1.2 );
//    real dSegL       = (MP.boxsize[2]-MP.boxsize[1] + MP.boxsize[0]) / (real)nbSegmentsL;
//    real dSegR       = (MP.boxsize[0] - (MP.boxsize[2]+MP.boxsize[1])) / (real)nbSegmentsR;
//    real dSegUp      = MP.boxsize[3] / (real)nbSegmentsUp;
//    printf( "nbSegL: %d, nbSegR: %d, nbSegUp: %d\n", nbSegmentsL, nbSegmentsR, nbSegmentsUp );
//    printf( "dSegL:  %f, dSegR:  %f, dSegUp:  %f\n", dSegL, dSegR, dSegUp );

    //rings left of the arm
    for( int ii = 0; ii <= nbSegmentsL; ii++ ) {
      if( nbSegmentsL > 0 )
        ringPos = (-MP.boxsize[0] + (MP.boxsize[2]-MP.boxsize[1] + MP.boxsize[0])/nbSegmentsL*ii);
      else
        ringPos = MP.boxsize[2]-MP.boxsize[1];
      glBegin(GL_LINE_LOOP);
      for( int aa = 0; aa < 2*pisteps; aa++ ) {
        glVertex3( ringPos, inflRad*cos(pifrac*aa), inflRad*sin(pifrac*aa) );
      }
      glEnd();
    }
    //rings right of the arm
    for( int ii = 0; ii <= nbSegmentsR; ii++ ) {
      if( nbSegmentsR > 0 )
        ringPos = (MP.boxsize[2]+MP.boxsize[1] + (MP.boxsize[0] - (MP.boxsize[2]+MP.boxsize[1]))/nbSegmentsR*ii);
      else
        ringPos = MP.boxsize[2]+MP.boxsize[1];
      glBegin(GL_LINE_LOOP);
      for( int aa = 0; aa < 2*pisteps; aa++ ) {
        glVertex3( ringPos, inflRad*cos(pifrac*aa), inflRad*sin(pifrac*aa) );
      }
      glEnd();
    }
    //rings on the arm
    for( int ii = 0; ii <= nbSegmentsUp; ii++ ) {
      if( nbSegmentsUp > 0 )
        ringPos = (MP.boxsize[1] + MP.boxsize[3]/nbSegmentsUp*ii);
      else
        ringPos = MP.boxsize[1];
      glBegin(GL_LINE_LOOP);
      for( int aa = 0; aa < 2*pisteps; aa++ ) {
        glVertex3( MP.boxsize[2]+inflRad*cos(pifrac*aa), ringPos, inflRad*sin(pifrac*aa) );
      }
      glEnd();
    }
    #endif

    break;
  }

  case SHAPE_ROUND_SQUARE:
    if ( DIM == 3 ) {
      displayBoxSection1( space, 0, MP.boxsize[3]-MP.boxsize[0], 0.1 );
      displayBoxSection1( space, 0, MP.boxsize[0]-MP.boxsize[3], 0.1 );
      displayBoxSection1( space, 1, MP.boxsize[3]-MP.boxsize[1], 0.1 );
      displayBoxSection1( space, 1, MP.boxsize[1]-MP.boxsize[3], 0.1 );
      displayBoxSection1( space, 2, MP.boxsize[2]-MP.boxsize[3], 0.1 );
      displayBoxSection1( space, 2, MP.boxsize[3]-MP.boxsize[2], 0.1 );
    } else {
      displayBoxSection1( space, 2, 0, 0.1 );
    }
    break;

  default:
    displayBoxSection1( space, 2, 0, 0.1 );
    if ( DIM == 3 ) {
      displayBoxSection1( space, 0, 0, 0.1 );
      displayBoxSection1( space, 1, 0, 0.1 );
    }
  break;
  }
}

//---------------------------------------------------------------------------
void Player::drawFrame1()
{
  int ii;
  real ab;
  Vecteur xx, yy;
  //printf("disp %i\n", sim.frameInBuffer() );

  static GLUquadricObj* qobj;        // a pointer to quadric objects
  static bool initialised = 0;
  if( !initialised ) {
    qobj = gluNewQuadric();
    gluQuadricDrawStyle(qobj, GLU_SILHOUETTE);
    gluQuadricNormals(qobj, GLU_NONE);
    initialised = 1;
  }

  //----------------- draw a frame around the box:
  if ( PP.boxwidth && Microtub::space )
    displayBox1(Microtub::space);

 //--------------Display gradient
  if ( PP.showGrad && MP.fieldTyp[0]){
  	displayGradient( PP.showGrad );
  }



   //================================================================
  //---------------------------- microtubules ----------------------
  ///\todo display arrows on microtubules

  Microtub* mt;
  glColor( PP.mtcolor[0] );

  //display lines only:
  if( PP.mtwidth && PP.mtlines && !PP.mtratchets && !PP.mtrainbow ) {
    glLineWidth( PP.mtwidth );
    for( mt = sim.firstMicrotub(); mt; mt = mt->next())
      if ( mt->mtColor >= 0) {
        glBegin(GL_LINE_STRIP);
        glColor( PP.mtcolor[ mt->mtColor ] );
        for( ii = 0; ii < mt->nbPoints(); ++ii )
          glVertexPseudo( mt->pointAddr( ii ));
        glEnd();
      }
  }

  //display segments which gradient of color indicating polarity
  if( PP.mtwidth && PP.mtratchets && !PP.mtrainbow ) {
    glLineWidth( PP.mtwidth );
    for( mt = sim.firstMicrotub(); mt; mt = mt->next())
      if ( mt->mtColor >= 0) {
        glBegin(GL_LINE_STRIP);
        glColor( PP.bgcolor );
        glVertexPseudo( mt->pointAddr( 0 ) );
        for( ii = 1; ii < mt->nbSegments(); ++ii ) {
          glColor( PP.mtcolor[ mt->mtColor ] );
          glVertexPseudo( mt->pointAddr( ii ) );
          glColor( PP.bgcolor );
          glVertexPseudo( mt->pointAddr( ii ) );
        }
        glColor( PP.mtcolor[ mt->mtColor ] );
        glVertexPseudo( mt->pointAddr( ii ) );
        glEnd();
      }
  }

  //display fix points of the tubes (every micro-meter):
  if( PP.mtspeckles ) {
    glPointSize( PP.mtsize );
    glBegin(GL_POINTS);
    for( mt = sim.firstMicrotub(); mt; mt = mt->next())
      if ( mt->mtColor >= 0) {
        glColor( PP.mtcolor[ mt->mtColor ] );
        ab = PP.mtspeckledist * ceil( mt->abscissa(MT_MINUS_END) / PP.mtspeckledist );
        while( ab <= mt->abscissa(MT_PLUS_END) ) {
          glVertexPseudo( mt->where(ab, MT_ORIGIN) );
          ab += PP.mtspeckledist;
        }
      }
    glEnd();
  }

  //display simulated points of the tubes:
  if( PP.mtpoints ) {
    glPointSize( PP.mtsize );
    glBegin(GL_POINTS);
    for( mt = sim.firstMicrotub(); mt; mt = mt->next())
      if ( mt->mtColor >= 0) {
        glColor( PP.mtcolor[ mt->mtColor ] );
        for( ii = 0; ii < mt->nbPoints(); ++ii )
          glVertexPseudo( mt->pointAddr( ii ) );
      }
    glEnd();
  }

  //display points at minus-ends:
  if( PP.mtends[MT_MINUS_END] ) {
    glPointSize( 3*PP.mtsize );
    glBegin(GL_POINTS);
    for( mt = sim.firstMicrotub(); mt; mt = mt->next())
      if ( mt->mtColor >= 0) {
        glColor( PP.mtcolor[ mt->getDynamicState(MT_MINUS_END) ] );
        glVertexPseudo( mt->pointAddr( 0 ) );
      }
    glEnd();
  }

  //display points at plus-ends:
  if( PP.mtends[MT_PLUS_END] ) {
    glPointSize( 3*PP.mtsize );
    glBegin(GL_POINTS);
    for( mt = sim.firstMicrotub(); mt; mt = mt->next())
      if ( mt->mtColor >= 0) {
        glColor( PP.mtcolor[ mt->getDynamicState(MT_PLUS_END) ] );
        glVertexPseudo( mt->pointAddr( mt->nbSegments() ) );
      }
    glEnd();
  }

  //display forces acting on the tubes:
  if( PP.mtwidth && PP.mtforces ) {
    glLineWidth( PP.mtwidth );
    glColor( PP.mtfcolor );
    glBegin(GL_LINES);
    for( mt = sim.firstMicrotub(); mt; mt = mt->next())
      if (( mt->mtColor >= 0 ) && ( mt->getForcesAddr() )) {
        for( ii = 0; ii < mt->nbPoints(); ++ii ) {
          glVertex( mt->pointAddr( ii ) );
          xx = mt->whereP(ii) + PP.mtforces * mt->forcesP(ii);
          glVertex( xx );
        }
      }
    glEnd();
  }

  //display projection of the points on the box:
  if ( PP.mtwidth && PP.mtproject ) {
    glColor( PP.mtfcolor );
    glLineWidth( PP.mtwidth );
    glPointSize( PP.mtsize );
    for( mt = sim.firstMicrotub(); mt; mt = mt->next())
      if ( mt->mtColor >= 0 )
        for( ii = 0; ii < mt->nbPoints(); ++ii ) {
          xx = mt->whereP( ii );
          if ( Microtub::space->isOutside( xx ) ) {
            Microtub::space->project( xx, yy );
            glBegin(GL_LINES);
            glVertex( xx );
            glVertex( yy );
            glEnd();
            glBegin(GL_POINTS);
            glVertex( yy );
            glEnd();
          }
        }
  }


  //display the internal tensions of the microtubules:
  if ( PP.width && PP.mtrainbow ) {
    glLineWidth( PP.mtwidth );
    for( mt = sim.firstMicrotub(); mt; mt = mt->next())
      if ( mt->mtColor >= 0 ) {
        glBegin(GL_LINES);
        xx = mt->whereP( 0 );
        for( ii = 0; ii < mt->nbSegments(); ++ii ) {
          glColorRainbow( 0.2 - mt->getLagrange(ii) * PP.mtrainbow );
          glVertexPseudo( xx );
          xx = mt->whereP( ii+1 );
          glVertexPseudo( xx );
        }
        glEnd();
      }
  }

  //================================================================
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

  //================================================================
  //---------------------------- nuclei ----------------------------

  //display the nuclear envelope
  if ( PP.nuenvelope ) {
    glLineWidth( PP.nuenvwidth );
    glColor(PP.nuenvcolor);
    for(Nucleus * nu=sim.firstNucleus(); nu ; nu=nu->next()) {
      glMatrixMode(GL_MODELVIEW);
      glPushMatrix();
      if(DIM == 2) {
        glTranslate(nu->coord(0), nu->coord(1), 0.0);
        gluDisk(qobj, 0.0, (GLdouble)nu->radius(), 100, 1);
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
  }

  //display the points on the nucleus (without the center-point)
  if ( PP.nupoints ) {
    glPointSize(PP.nuptsize);
    glBegin(GL_POINTS);
    glColor(PP.nuptcolor);
    for(Nucleus * nu=sim.firstNucleus(); nu ; nu=nu->next()) {
      for( ii = 1; ii < nu->nbPoints()-nu->nbrefpts; ++ii )
        glVertex( nu->pointAddr( ii ) );
    }
    glEnd();
  }

  //display the reference point
  if ( PP.nurefpts ) {
    glPointSize(PP.nuptsize);
    glBegin(GL_POINTS);
    glColor(PP.nurefcolor);
    for(Nucleus * nu=sim.firstNucleus(); nu ; nu=nu->next()) {
      for( ii = nu->nbPoints()-nu->nbrefpts; ii < nu->nbPoints(); ++ii )
      glVertex( nu->pointAddr( ii ) );
    }
    glEnd();
  }

  //display the mtoc-microtubule links
  if ( PP.nulinks ) {
    PointExact pte1, pte2;
    glLineWidth( PP.mtwidth );
    glBegin(GL_LINES);
    for(Nucleus * nu=sim.firstNucleus(); nu ; nu=nu->next())
      for( int nuindx = 0 ; nuindx < nu->maxNbMicrotub(); nuindx++ ) {
	if ( nu -> setPointClamp( &pte1, &pte2, nuindx )) {
	  glColor( PP.nuptcolor );
	  glVertex( pte1.where() );
	  glColor( PP.nulinkcolor );
	  glVertex( pte2.where() );
	}
      }
    glEnd();
  }

  //================================================================
  //--------------------------- grafted ----------------------------

  Grafted * gh;
  glPointSize(PP.ghsize);


  //display the attached position:
  if ( PP.ghhands && ( PP.ghselect & 1 )) {
    glBegin(GL_POINTS);
    for(gh=sim.firstFreeGrafted(); gh ; gh=gh->next()) {
      glColor( PP.hacolor[ gh->getType() ], 2 );
      glVertex( gh->whereGrafted() );
    }
    glEnd();
  }

  // display the grafted positions
  if ( PP.ghhands && ( PP.ghselect & 2 )) {
    glBegin(GL_POINTS);
    for(gh=sim.firstBoundGrafted(); gh ; gh=gh->next()) {
      if ( PP.ghrainbow )
        glColorRainbow( gh->getStretch().norm() * PP.ghrainbow );
      else
        glColor( PP.hacolor[ gh->getType() ] );
      glVertex( gh->whereHand() );
    }
    glEnd();
  }

  // display the grafted link or forces
  if (( PP.ghlinks || PP.ghrainbow ) && ( PP.ghselect & 2 )) {
    glLineWidth( PP.ghwidth );
    glBegin(GL_LINES);
    for(gh=sim.firstBoundGrafted(); gh ; gh=gh->next()) {
      Space::moduloNear( xx, gh->whereHand(), gh->whereGrafted() );

      if ( PP.ghrainbow ) {
        glColorRainbow( gh->getStretch().norm() * PP.ghrainbow );
        glVertex( gh->whereGrafted() );
      } else {
        glColor( PP.hacolor[ gh->getType() ], 2 );
        glVertex( gh->whereGrafted() );
        glColor( PP.hacolor[ gh->getType() ] );
      }
      glVertex( xx );
    }
    glEnd();
  }

  // display the grafted displacements
  if ( PP.ghdisp && ( PP.ghselect & 1 )) {
    glLineWidth( PP.ghwidth );
    glBegin(GL_LINES);
    for(gh=sim.firstFreeGrafted(); gh ; gh=gh->next()) {
      glColor( PP.hacolor[ gh->getType() ], 2 );
      glVertex( gh->getPosition() );
      glColor( PP.bgcolor );
      glVertex( gh->getPositionOld() );
      gh->setPositionOld();
    }
    glEnd();
  }


  if ( PP.ghdisp && ( PP.ghselect & 2 )) {
    glLineWidth( PP.ghwidth );
    glBegin(GL_LINES);
    for(gh=sim.firstBoundGrafted(); gh ; gh=gh->next()) {
      glColor( PP.hacolor[ gh->getType() ] );
      glVertex( gh->getPosition() );
      glColor( PP.bgcolor );
      glVertex( gh->getPositionOld() );
      gh->setPositionOld();
    }
    glEnd();
  }

  //display a flash (a point of double size) for a recent change of state
  if ( PP.ghflash && PP.live ) {
    glPointSize(2 * PP.ghsize);
    glBegin(GL_POINTS);
    for(gh=sim.firstFreeGrafted(); gh ; gh=gh->next())
      if ( gh->isStateNew() ) {
        glColor( PP.hacolor[ gh->getType() ], 1 );
        glVertex( gh->whereGrafted() );
      }
    glEnd();
  }

  /*
  //display a flash (a point of double size) for change of state
  if ( PP.ghflash && PP.live  ) {
    glPointSize(2 * PP.ghsize);
    glBegin(GL_POINTS);
    for(gh=sim.firstBoundGrafted(); gh ; gh=gh->next())
      if ( gh->isStateNew() ) {
        glColor( PP.hacolor[ gh->getType() ] );
        glVertex( gh->whereHand() );
      }
      glEnd();
  }
  */

  //================================================================
  //------------------------- complexes ---------------------------

  Complex * cx;
  glPointSize(PP.cxsize);


  // display free complexes
  if ( PP.cxhands && ( PP.cxselect & 1 )) {
    glBegin(GL_POINTS);
    for(cx=sim.firstFreeComplex(); cx ; cx=cx->next()) {
      //glColor( PP.hacolor[ cx->getType1() ], PP.hacolor[ cx->getType2() ], 2);   // Nk: 31 - 03 -2016
	glColor( PP.cxcolor[ cx->getType1() ], PP.cxcolor[ cx->getType2() ], 2);
      glVertex( cx->getPosition() );
    }
    glEnd();
  }


  // display bound complexes
  if ( PP.cxhands && ( PP.cxselect & 2 )) {
    glBegin(GL_POINTS);
    for(cx=sim.firstBoundComplex(); cx ; cx=cx->next()) {
      if ( cx->isAttached1() )
        //glColor( PP.hacolor[ cx->getType1() ], 1 );   // NK: 31 -03 - 2016
	glColor( PP.cxcolor[ cx->getType1() ], PP.cxcolor[ cx->getType2() ], 2);
      else
        //glColor( PP.hacolor[ cx->getType2() ], 1 );  // NK: 31 -03 - 2016
	glColor( PP.cxcolor[ cx->getType1() ], PP.cxcolor[ cx->getType2() ], 2);
      glVertex( cx->calculatePosition() );
    }
    glEnd();
  }


  // display bridging complexes
  if ( PP.cxhands && ( PP.cxselect & 4 )) {
    glBegin(GL_POINTS);
    for(cx=sim.firstBridgeComplex(); cx ; cx=cx->next()) {
      if ( PP.cxrainbow ) {
        glColorRainbow( cx->getStretch().norm() * PP.cxrainbow );
        glVertex( cx->where1() );
      } else {
        //glColor( PP.hacolor[ cx->getType1() ] );   // NK: 31 -03 -2016
	glColor( PP.cxcolor[ cx->getType1() ], PP.cxcolor[ cx->getType2() ], 2);
        glVertex( cx->where1() );
        //glColor( PP.hacolor[ cx->getType2() ] );   // NK: 31 -03-2016
	glColor( PP.cxcolor[ cx->getType1() ], PP.cxcolor[ cx->getType2() ], 2);
      }
      glVertex( cx->where2() );
    }
    glEnd();
  }


  // display the link for bridging complexes
  if ( PP.cxlinks && ( PP.cxselect & 4 )) {
    glLineWidth( PP.cxwidth );
    glBegin(GL_LINES);
    for(cx=sim.firstBridgeComplex(); cx ; cx=cx->next()) {
      xx = cx->where1();
      yy = cx->where2();
      Space::moduloNear( yy, xx );

      if ( PP.cxrainbow ) {
        glColorRainbow( cx->getStretch().norm() * PP.cxrainbow );
        glVertex( xx );
      } else {
        glColor( PP.cxcolor[ cx->getType1() ] );  // Neha Khetan, 10 July 2018 - Modifying the color of links of the complex so that it is same as the color of complex - Changed hacolor to cxcolor
        glVertex( xx );
        glColor( PP.cxcolor[ cx->getType2() ] );
// Neha Khetan, 10 July 2018 - Modifying the color of links of the complex so that it is same as the color of complex - Changed hacolor to cxcolor
      }
      glVertex( yy );
    }
    glEnd();
  }



  //display the displacement of motor complexes
  if ( PP.cxdisp ) {

    glLineWidth( PP.cxwidth );
    glBegin(GL_LINES);

    for(cx=sim.firstFreeComplex(); cx ; cx=cx->next()) {
      if ( PP.cxselect & 1 ) {
        glColor( PP.hacolor[ cx->getType1() ], PP.hacolor[ cx->getType2() ], 2);
        glVertex( cx->getPosition() );
        glColor( PP.bgcolor );
        glVertex( cx->getPositionOld() );
      }
      cx->setPositionOld();
    }

    for(cx=sim.firstBoundComplex(); cx ; cx=cx->next()) {
      cx->getPosition();
      if ( PP.cxselect & 2 ) {
        glColor( PP.hacolor[ cx->getType1() ], 1 );
        glVertex( cx->getPosition() );
        glColor( PP.bgcolor );
        glVertex( cx->getPositionOld() );
      }
      cx->setPositionOld();
    }

    for(cx=sim.firstBridgeComplex(); cx ; cx=cx->next()) {
      cx->getPosition();
      if ( PP.cxselect & 4 ) {
        glColor( PP.hacolor[ cx->getType1() ], PP.hacolor[ cx->getType2() ], 0);
        glVertex( cx->getPosition() );
        glColor( PP.bgcolor );
        glVertex( cx->getPositionOld() );
      }
      cx->setPositionOld();
    }
    glEnd();
  }


  //display a flash (a point of double size) for a recent change of state
  if ( PP.cxflash && PP.live ) {
    glPointSize(2 * PP.cxsize);
    glBegin(GL_POINTS);
    for(cx=sim.firstFreeComplex(); cx ; cx=cx->next())
      if ( cx->isStateNew() ) {
        glColor( PP.hacolor[ cx->getType1() ] );
        glVertex( cx->getPosition() );
      }
    ///\todo flash display for the complexes: only when they break, not when they detach
    glEnd();
  }

}

