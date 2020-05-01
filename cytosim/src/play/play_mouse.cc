//RCS: $Id: play_mouse.cc,v 2.7 2005/03/15 19:58:09 nedelec Exp $


Player::MOUSE_ACTION Player::mouseAction(const int button, const int modifier ) {
  switch( button ) {
  case GLUT_LEFT_BUTTON:
    switch( modifier ) {
    case GLUT_ACTIVE_CTRL: 
      return MOUSE_SPIN;
    case GLUT_ACTIVE_SHIFT:
      return (DIM<3) ? MOUSE_GRAB : MOUSE_TRANSLATE;
    case GLUT_ACTIVE_SHIFT + GLUT_ACTIVE_CTRL: 
      return MOUSE_TRANSLATE_AXIS;
    case GLUT_ACTIVE_ALT: 
      return MOUSE_NOTHING;
    case 0:
      return (DIM<3) ? MOUSE_TRANSLATE : MOUSE_ROTATE;
    }
  case GLUT_MIDDLE_BUTTON:
    if ( modifier & GLUT_ACTIVE_SHIFT ) 
      return (DIM<3)? MOUSE_SPIN : MOUSE_ROTATE;
    return MOUSE_ZOOM;
  }
  return MOUSE_NOTHING;
}



//--------------------------------------------------------------------------
//-----------------------------  mouse -------------------------------------
//--------------------------------------------------------------------------

void Player::processMouse(const int button, const int state, const int x, const int y)
{
  if ( state != GLUT_DOWN ) { mouse_action = MOUSE_NOTHING; return; }
  //printf("mouse button %i %i %i %i\n", button, state, x, y);
  
  mouse_action = mouseAction( button, glutGetModifiers() );
  if (( mouse_action == MOUSE_GRAB ) && ( !sim.isReady() ))
    mouse_action = MOUSE_NOTHING;
  if ( mouse_action == MOUSE_NOTHING ) return;
  
  glGetDoublev(GL_MODELVIEW_MATRIX, modelviewMat);
  
  GLdouble dx, dy, dz;
  gluUnProject(x, viewport[3]-y, 0, modelviewMat, 
	       projectionMat, viewport, &dx, &dy, &dz);
  
  switch( mouse_action ) {
    
  case MOUSE_TRANSLATE: {
    mouse_origin.set( dx, dy, dz );
    focus_save = PP.focus;
  } break;
  
  case MOUSE_TRANSLATE_AXIS: {
    mouse_origin.set( dx, dy, dz );
    GLdouble cx, cy, cz;
    gluUnProject(viewport[2]/2.0, viewport[3]/2.0, 0, modelviewMat,
                 projectionMat, viewport, &cx, &cy, &cz);
    mouse_normal = Vecteur3(cx, cy, cz) - PP.focus;
    mouse_normal.normalize();
    gluUnProject(viewport[2]/2.0, viewport[3], 0, modelviewMat,
		 projectionMat, viewport, &dx, &dy, &dz);
    mouse_axis = Vecteur3( dx-cx, dy-cy, dz-cz ).normalized();
    focus_save = PP.focus;
  } break;
  
  case MOUSE_ROTATE: {
    mouse_origin.set( dx, dy, dz );
    mouse_normal = mouse_origin - Vecteur3( PP.focus );
    mouse_normal *= mouse_amplification / mouse_normal.normSquare();
    quat_save = PP.quat;
  } break;
  
  case MOUSE_SPIN: {
    GLdouble cx, cy, cz;
    gluUnProject(viewport[2]/2.0, viewport[3]/2.0, 0, modelviewMat,
		 projectionMat, viewport, &cx, &cy, &cz);
    mouse_origin.set( cx, cy, cz );
    mouse_axis = mouse_origin - Vecteur3( PP.focus );
    mouse_axis.normalize();
    mouse_normal = Vecteur3( dx, dy, dz ) - mouse_origin;
    quat_save = PP.quat;
  } break;
  
  case MOUSE_ZOOM: {		
    real xx = x - viewport[2]/2.0;
    real yy = y - viewport[3]/2.0;
    mouse_zoom_scalar = sqrt( xx*xx + yy*yy );
    if ( mouse_zoom_scalar > 0 ) 
      mouse_zoom_scalar = 1.0 / mouse_zoom_scalar;
    zoom_save = PP.zoom;
  } break;
  
  case MOUSE_GRAB: {  //grab MTs with the mouse, a 2D only feature so far
    if ( playGrafted == 0 ) {
      //register a grafted hand with a dummy type MAGIC_TYPE:
      playGrafted = sim.link( new Grafted( MAGIC_TYPE, Vecteur(dx, dy, 0)));
      playGraftedList.pushFront(playGrafted);
      
      //set its properties to make it bind as fast as possible:
      MP.haattachdist[ MAGIC_TYPE ]    = MP.MHGL;
      MP.haattachrate_dt[ MAGIC_TYPE ] = 1;
    } else {
      //find the closest grafted that is below the threshold:
      real distance = playGraftedDistanceThreshold;
      for(int ii=0; ii < playGraftedList.size(); ++ii) {
        real d = (playGraftedList[ii] -> whereGrafted() - Vecteur(dx, dy, 0)).normSquare();
        if ( d < distance ) {
          playGrafted = playGraftedList[ii];
          distance = d;
        }
      }
      if ( playGrafted ) {
        if (( distance >= playGraftedDistanceThreshold ) && playGrafted->isAttached() ) 
          playGrafted -> detach(DETACH_SPONTAN);
        playGrafted -> setPosition( Vecteur(dx, dy, 0) );
      }
    }
  } break;
  
  }
}



void Player::processMotion(const int x, const int y)
{
  if ( mouse_action == MOUSE_NOTHING ) return;
  //printf("mouse motion %i %i %i\n", mouse_action, x, y);

  GLdouble dx, dy, dz;
  gluUnProject(x, viewport[3]-y, 0, modelviewMat, 
			   projectionMat, viewport, &dx, &dy, &dz);
  Vecteur3 mouse_drag = Vecteur3( dx, dy, dz ) - mouse_origin;

  switch( mouse_action ) {

  case MOUSE_ROTATE: {
    PP.quat.setFromAxisAngle( mouse_normal ^ mouse_drag );
    PP.quat.leftMult( quat_save );
    flashText("Quat %.2f %.2f %.2f %.2f", PP.quat[0], PP.quat[1], PP.quat[2], PP.quat[3]);
    setModelView(1);
  } break;
  
  case MOUSE_SPIN: {
    real cos = mouse_normal * mouse_drag;
    real sin = ( mouse_normal ^ mouse_drag ) * mouse_axis;
    PP.quat.setFromAxisAngle( mouse_axis, atan2( sin, cos ) );
    PP.quat.leftMult( quat_save );
    flashText("Quat %.2f %.2f %.2f %.2f", PP.quat[0], PP.quat[1], PP.quat[2], PP.quat[3]);
    setModelView(1);
  } break;
  
  case MOUSE_TRANSLATE: {
    PP.focus = focus_save - mouse_drag;
    flashText("Focus %.2f %.2f %.2f", PP.focus[0], PP.focus[1], PP.focus[2]);
    setModelView(1);
  } break;
  
  case MOUSE_TRANSLATE_AXIS: {
    real S = mouse_drag * mouse_axis;
    mouse_drag -= S * ( mouse_normal + mouse_axis );
    PP.focus = focus_save - mouse_drag;
    flashText("Focus %.2f %.2f %.2f", PP.focus[0], PP.focus[1], PP.focus[2]);
    setModelView(1);
  } break;
  
  case MOUSE_ZOOM: {
    real xx = x - viewport[2]/2.0;
    real yy = y - viewport[3]/2.0;
    real Z = mouse_zoom_scalar * sqrt( xx*xx + yy*yy );
    if ( Z > 0 ) {
      PP.zoom = zoom_save * Z;
      flashText("Zoom %.2f", PP.zoom);
      setModelView(1);
    }
  } break;
  
  case MOUSE_GRAB: if ( playGrafted ) {
    playGrafted->setPosition( Vecteur(dx, dy, dz) );
  } break;
  }
}
