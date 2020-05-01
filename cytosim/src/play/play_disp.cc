//RCS: $Id: play_disp.cc,v 2.27 2005/04/14 11:20:46 foethke Exp $
//----------------- calculate the color / visibility of objects:
void Player::processState(){
  //if this is a periodic space, we bring all in the same representation
  if ( Space::isPeriodic() )
    sim.moduloPosition();

  //--------------Display gradient
  if ( PP.showGrad && MP.fieldTyp[0]){
  	displayGradient( PP.showGrad );
  }


  //set the color parameters of microtubules: mtcolor = -1 makes them invisible
  for(Microtub * mt = sim.firstMicrotub(); mt; mt = mt->next()) {

    switch( PP.mtcolormode ) {

      case 0:
        mt -> mtColor = mt -> getType() % 10;
        break;

      case 1:
        //we use the memory location of the organizer modulo 10...
        //this gives unpredictable colors, but it works!
        mt -> mtColor = long( mt->getOrganizer() ) % 10;
        break;

      case 2:
        mt -> mtColor = mt -> getDynamicState( MT_PLUS_END );
        break;
    }

    // check is we are showing specific colors
    if ((( mt->mtColor == 0 ) && !( PP.mtselect & 1 ))
        || (( mt->mtColor == 1 ) && !( PP.mtselect & 2 ))
        || (( mt->mtColor == 2 ) && !( PP.mtselect & 4 )))
      mt->mtColor = -1;

  }
}


//---------------------------------------------------------------------------
void Player::drawFrame()
{
  switch( PP.style ) {
  case 3: drawFrame3(); break; //solid rendering
  case 2: drawFrame2(); break; //not used
  default:
  case 1: drawFrame1(); break; //wireframe
  }
  glPrintErrors("drawFrame()");
}


//---------------------------------------------------------------------------
void Player::tileFrame()
{
  glMatrixMode(GL_MODELVIEW);

  if ( MP.boxshape == SHAPE_PERIODIC ) {
    //display four times:
    for( int dx = 0; dx <= 1; ++dx )
      for( int dy = 0; dy <= 1; ++dy )
        if ( dx || dy ) {
          glTranslate( 2*dx * MP.boxsize[0] , 2 * dy * MP.boxsize[1], 0 );
          drawFrame();
          glTranslate( -2*dx * MP.boxsize[0] , -2 * dy * MP.boxsize[1], 0 );
        }
  }

  if ( MP.boxshape == SHAPE_STRIP ) {
    //display three times:
    for( int dx = -1; dx <= 1; dx+=2 ) {
      glTranslate( 2*dx * MP.boxsize[0] , 0, 0 );
      drawFrame();
      glTranslate( -2*dx * MP.boxsize[0] , 0, 0 );
    }
  }
}


//---------------------------------------------------------------------------
void Player::displaySubScaleBar(const int vertical, const real scale, const real width )
{
  glLineWidth(1);
  glBegin(GL_LINES);
  glColor3f(0, 0, 1.0);
  if ( vertical ) {
    for(int ii = -5; ii <= 5; ++ii) {
      glVertex2( -width, ii*scale );
      glVertex2(      0, ii*scale );
    }
  } else {
    for(int ii = -5; ii <= 5; ++ii) {
      glVertex2( ii*scale,     0 );
      glVertex2( ii*scale, width );
    }
  }
  glEnd();
}

void Player::displayScaleBar(const int vertical)
{
  glPushMatrix();
  glLoadIdentity();

  if ( vertical )
    glTranslate( visibleRegion[0]-1, 0, 0 );
  else
    glTranslate( 0, 1-visibleRegion[1], 0 );

  glScaled(PP.zoom, PP.zoom,  PP.zoom);

  glLineWidth(2);
  //draw a long box:
  glBegin(GL_LINE_LOOP);
  glColor3f(0, 0, 1.0);

  if ( vertical ) {
    glVertex2( -0.5, -5.0 );
    glVertex2(  0.0, -5.0 );
    glVertex2(  0.0,  5.0 );
    glVertex2( -0.5,  5.0 );
  } else {
    glVertex2( -5.0,  0.0 );
    glVertex2( -5.0,  0.5 );
    glVertex2(  5.0,  0.5 );
    glVertex2(  5.0,  0.0 );
  }
  glEnd();

  //draw lines every 1 um
  displaySubScaleBar(vertical, 1.0, 0.5);

  //draw small lines every 0.1 um
  displaySubScaleBar(vertical, 0.1, 0.25);

  //draw very tiny lines every 0.01 um
  displaySubScaleBar(vertical, 0.01, 0.10);

  //cleanup
  glPopMatrix();
}

//-------------------------------------------------------
//CHAITANYA: display gradient of catastrophe and rescue in color using openGL commands
void Player::displayGradient( int gradTyp ){
	Vecteur v;
	int maxC=sim.distMap.nbCells();
	real alphaR = PP.alphaDisp;
	//real alphaR = 0.5;//transparency of map-display
	//real x_cor = -MP.boxsize[0], y_cor = -MP.boxsize[1];

	//glEnable(GL_BLEND);
	//glBlendFunc(GL_ONE, GL_ONE );
	/*glBegin(GL_POLYGON);
		glVertex2f(0.0, 0.0);
		glVertex2f(5.0, 0.0);
		glVertex2f(5.0, 5.0);
		glVertex2f(0.0, 5.0);
	glEnd();*/

	//printf("num. cells: %i\n", sim.resMap.nbCells() );
	if(gradTyp == 1){
		real maxv= sim.resMap.getMaxval();
		real minv= sim.resMap.getMinval();
		real rDiff = maxv - minv;//range in values
		//printf("@disp, max: %f, min: %f, val[0]: %f\n", maxv, minv, sim.resMap[ 0 ]);
		for ( int n=0; n< maxC; n++ ){
			//glColor3f( 0.1, 0.5*sim.resMap[ n ]/MP.mtresgradient[1], 0.7);
			//glColor4f( sim.resMap[ n ]/MP.mtresgradient[1], sim.resMap[ n ]/MP.mtresgradient[1], sim.resMap[ n ]/MP.mtresgradient[1], 0.1);
			glColor4f(0.1*(sim.resMap[ n ]-minv)/rDiff, 1*(sim.resMap[ n ]-minv)/rDiff, 0.1*(sim.resMap[ n ]-minv)/rDiff, alphaR);
			v = sim.resMap.positionFromIndex( n, 0 );
			//printf("coord: %f, %f \n", v[0], v[1] );
			glRectf( v[0], v[1], (v[0] + 2*MP.ficellsize), (v[1] + 2*MP.ficellsize) ); //lowercorner_x, lowercorner_y, width, height
			//vertical
			if( DIM == 3){
				glBegin( GL_QUADS );
				glColor4f(0.1*(sim.resMap[ n ]-minv)/rDiff, 1*(sim.resMap[ n ]-minv)/rDiff, 0.1*(sim.resMap[ n ]-minv)/rDiff, alphaR);

				glVertex3f( v[0], 0.0, v[2]);
				glVertex3f( v[0], 0.0, v[2]+2*MP.ficellsize);
				glVertex3f( v[0]+2*MP.ficellsize, 0.0, v[2]+2*MP.ficellsize);
				glVertex3f( v[0]+2*MP.ficellsize, 0.0, v[2] );
				glEnd();
			}
		}
	}
	else if( gradTyp == 2 ){
		real maxv= sim.catMap.getMaxval();
		real minv= sim.catMap.getMinval();
		real rDiff = maxv - minv;//range in values
		//printf("@disp, max: %f, min: %f, val[0]: %f\n", maxv, minv, sim.catMap[ 0 ] );
		for ( int m=0; m< maxC; m++ ){
			//glColor3f( 0.5 * sim.catMap[ m ]/MP.mtcatagradient[0], 0.0, 0.3);
			//glColor4f( 0.5*sim.catMap[ m ]/MP.mtcatagradient[0], 0.5*sim.catMap[ m ]/MP.mtcatagradient[0], 0.5*sim.catMap[ m ]/MP.mtcatagradient[0], 0.1);
			glColor4f( 1.0*(sim.catMap[ m ]-minv)/rDiff, 0.1*(sim.catMap[ m ]-minv)/rDiff, 0.1*(sim.catMap[ m ]-minv)/rDiff, alphaR);
			v = sim.catMap.positionFromIndex( m, 0 );
			//printf("coord: %f, %f \n", v[0], v[1] );
			glRectf(v[0]  , v[1] , ( v[0] + 2*MP.ficellsize), ( v[1] + 2*MP.ficellsize ) ); //lowercorner_x, lowercorner_y, width, height
			//vertical
			if( DIM == 3){
				glBegin( GL_QUADS );
				glColor4f( 1.0*(sim.catMap[ m ]-minv)/rDiff, 0.1*(sim.catMap[ m ]-minv)/rDiff, 0.1*(sim.catMap[ m ]-minv)/rDiff, alphaR);

				glVertex3f( v[0], 0.0, v[2]);
				glVertex3f( v[0], 0.0, v[2]+2*MP.ficellsize);
				glVertex3f( v[0]+2*MP.ficellsize, 0.0, v[2]+2*MP.ficellsize);
				glVertex3f( v[0]+2*MP.ficellsize, 0.0, v[2] );
				glEnd();
			}
		}
	}
	else if( gradTyp == 3 ){
		real maxv= sim.distMap.getMaxval();
		real minv= sim.distMap.getMinval();
		real rDiff = maxv - minv;
		//printf("@disp, max: %f, min: %f, val[0]: %f\n", maxv, minv, sim.distMap[ 0 ] );
		for ( int m=0; m< maxC; m++ ){
			//glColor3f( 0.5 * sim.catMap[ m ]/MP.mtcatagradient[0], 0.0, 0.3);
			glColor4f( 0.5*(sim.distMap[ m ]-minv)/rDiff, 0.5*(sim.distMap[ m ]-minv)/rDiff, 0.5*(sim.distMap[ m ]-minv)/rDiff, alphaR);
			v = sim.distMap.positionFromIndex( m, 0 );
			//printf("coord: %f, %f \n", v[0], v[1] );
			glRectf(v[0]  , v[1] , ( v[0] + 2*MP.ficellsize), ( v[1] + 2*MP.ficellsize ) ); //lowercorner_x, lowercorner_y, width, height
			//vertical
			if( DIM == 3){
				glBegin( GL_QUADS );
				glColor4f( 0.5*(sim.distMap[ m ]-minv)/rDiff, 0.5*(sim.distMap[ m ]-minv)/rDiff, 0.5*(sim.distMap[ m ]-minv)/rDiff, alphaR);
				glVertex3f( v[0], 0.0, v[2]);
				glVertex3f( v[0], 0.0, v[2]+2*MP.ficellsize);
				glVertex3f( v[0]+2*MP.ficellsize, 0.0, v[2]+2*MP.ficellsize);
				glVertex3f( v[0]+2*MP.ficellsize, 0.0, v[2] );
				glEnd();
			}

		}
	}
    //glEnd();
}

//---------------------------------------------------------------------------
void Player::displayMessage()
{
  // build the string to be displayed:
  char text[STRING_SIZE] = "\0";
  unsigned long pos = 0;

  if( PP.showtag ) {
    if ( PP.tag[0] != '\0' )
      pos += snprintf(text+pos, sizeof(text)-pos, "%s", PP.tag );
    else
      pos += snprintf(text+pos, sizeof(text)-pos, "Cytosim ");
  }

  if ( PP.showframe && ( sim.frameInBuffer() >= 0 ))
    pos += snprintf(text+pos, sizeof(text)-pos, " %3i", sim.frameInBuffer());

  if ( PP.showtime )
    pos += snprintf(text+pos, sizeof(text)-pos, " %8.3fs", sim.simTime());

  //report the current status of user_control_key[]
  for( int ii = 0; ii < 5; ++ii )
    if ( user_control_keys[ii] )
      pos += snprintf(text+pos, sizeof(text)-pos, " ALT-%i", ii);

  if ( PP.live > 1 )
    pos += snprintf(text+pos, sizeof(text)-pos, " 2^%i steps", PP.live-1 );
  else
    if ( reader.lastFrameBeforeEof() <= PP.frame )
      pos += snprintf(text+pos, sizeof(text)-pos, " eof");

  //display frame info / comments
  if ( PP.showinfo )
    pos += snprintf(text+pos, sizeof(text)-pos, "\n%s", sim.frameInfo() );

  displayText( text, 0x000FFFFF, PP.width, PP.height, 0 );

  //we should check more carefully for string overflow
  if ( pos >= sizeof(text) )
    printf("play_disp.cc::displayMessage(): string overflow %lu > sizeof(text)\n", pos);
}

//---------------------------------------------------------------------------
void Player::displayParameters(bool displayValues)
{
  char text[4096];
  if ( displayValues ) {
    //this displays the values of the parameters used
    MP.printLine( text, sizeof(text) );
    displayText( text,  0x555555FF, PP.width, PP.height, 4 );
  } else {
    //this display the file from which the parameters where set,
    // even if this file has been modified since!
    FILE * paramFile = fopen( MP.lastFileParsed(), "r" );
    if ( paramFile ) {
      if ( ! ferror( paramFile )) {
        IO.readLine( text, sizeof(text), paramFile, false );
        displayText( text,  0x555555FF, PP.width, PP.height, 4 );
      }
      fclose( paramFile );
    }
  }
}


//-------------------------------  display  ---------------------------------
// reads the frame if necessary, clear the buffer, and add time /frame tag
//---------------------------------------------------------------------------
void Player::display()
{
  //MSG(9, "display() PP.frame %i frameInBuffer %i eofFrame %i\n", PP.frame, sim.frameInBuffer(), reader.firstEofFrame());
  //----------------------load the frame in sim. if necessary:
  // PP.frame is the one requested, while sim.frameInBuffer() is in the buffer
  if ( sim.frameInBuffer() != PP.frame ) {
    MSG(9, "play::display() loading frame %i from file\n", PP.frame);
    try {
      if ( NO_ERROR == reader.readFrameNoCatch( PP.frame )) {
        //read was successful!
        //could perform post-processing of the simulation state here:

      } else {
        if ( reader.getInputFile() )
          flashText("Error reading frame %i", PP.frame);
      }
    } catch( IOException e ) {
      flashText("Error reading frame %i:\n%s", PP.frame, e.getMessage());
    }
  }

  //-------------------autozoom, usually on the first occurence:
  if ( PP.autozoom && Microtub::space ) {
      Vecteur bounds = Microtub::space->getBoundingRect();
      //if the dimensions of the surrounding box are not zero:
      if (( bounds[0] > 0 ) && ( bounds[1] > 0 )) {
          //we adjust the zoom to fit the box:
          PP.zoom *= minT( visibleRegion[0] / bounds[0], visibleRegion[1] / bounds[1] );
          PP.autozoom = 0;
          setModelView();
      }
  }


  //----------------- clear pixels:

  glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );



  //----------------- calculate the color / visibility of objects:

  processState();

  //----------------- draw objects:

  drawFrame();


  //--------------Display gradient
    if ( PP.showGrad && MP.fieldTyp[0]){
      	displayGradient( PP.showGrad );
    }



  //----------------- periodic duplication of image:

  if ( PP.periodic ) tileFrame();

  //------------------write time/frame info:

  //display some info:
  displayMessage();

  if ( ! PP.offscreen ) {
    //show the description of key bindings
    if ( PP.showhelp )
      displayText( descriptionOfKeys(), 0x55551FFF, PP.width, PP.height, 3 );

    //display the flash text in a random color:
    if ( PP.showflashtext > 0 ) {
      --PP.showflashtext; //count-down timer: at zero the text is hidden
      displayText( PP.flashtext, 0x8888FFFF, PP.width, PP.height, 2 );
    }
  }

  //--------------- add a scale bar:
  if ( PP.scalebar )
    displayScaleBar( PP.width > PP.height );


  //--------------- display the parameters on the window
  if ( PP.showparameters )
    displayParameters();

  //cleanup OpenGL
  glFinish();
  if ( !PP.offscreen )
    glutSwapBuffers();
}
