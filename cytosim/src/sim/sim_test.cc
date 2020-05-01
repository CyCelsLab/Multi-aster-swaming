//RCS: $Id: sim_test.cc,v 2.7 2005/02/02 13:30:20 nedelec Exp $
 
//TODO: replace all this by a XML parser for initial configurations

void SIM::pinch1( Microtub * mt, const real ab, const Vecteur shift, const int type)
{
  link( new Grafted(type,  mt->where(ab, MT_MINUS_END) + shift) )
    -> attachTo( mt, ab, MT_MINUS_END );
}

void SIM::pinch2( Microtub * mt, const real ab, const Vecteur shift, const int type)
{
  link( new Grafted(type,  mt->where(ab, MT_PLUS_END) + shift) )
  -> attachTo( mt, ab, MT_PLUS_END );
}


//========================================================================
//==================== TEST 1     MP.initcode=1??  =======================
//========================================================================


//diverse arrangements of microtubules, motors and complexes
void SIM::populateTest1()
{
  Vecteur shift=Vecteur(0, MP.boxsize[4], 0);
  Microtub * mt=0, * mt1=0, * mt2=0;
  Complex  * co=0;
  Vecteur   dir, place;
  real      len;

  for( int ii=0; ii < MP.mtmax; ++ii )
    {
      dir   = Microtub::initDirection();
      len   = Microtub::initLength();
      place = Microtub::initPosition( Microtub::space, digit(MP.mtinit,1)) - dir*len/2;

      switch( digit(MP.initcode, 1) ) {
 
        case 0:                            //pinched in the middle
          mt = link( new Microtub(place, dir, len) );
          pinch1( mt, len/2.0, shift, 0); 
          break;
          
        case 1:                            //pinched at minus end
          mt = link( new Microtub(place, dir, len) );
          pinch1( mt, 0, shift, 0);
          break;
	
        case 2:                            //pinched at both ends (compression)
          mt = link( new Microtub(place, dir, len) );
          shift = dir * MP.boxsize[4];
          pinch1( mt, 0,   shift, 0);
          pinch1( mt, len,-shift, 1);
          break;

      case 3:                            //pinched at both ends (traction)
        mt = link( new Microtub(place, dir, len) );
        shift = dir * MP.boxsize[4];
        pinch1( mt, 0,  -shift, 0);
        pinch1( mt, len, shift, 1);
        break;
        
       
      case 4:                            //pinched both size and middle
        mt = link( new Microtub(place, dir, len) );
        
        pinch1( mt, 0,  Vecteur(MP.boxsize[4], MP.boxsize[5], 0), 0);
        pinch2( mt, 0,  Vecteur(-MP.boxsize[4], MP.boxsize[5], 0), 0);
        pinch1( mt, len/2.0, Vecteur(0, -MP.boxsize[5], 0), 0);
        pinch1( mt, len/2.0, Vecteur(0, -MP.boxsize[5], 0), 0);
        break;
	
      case 5: { //parallel & antiparallel pairs of microtubules
    
        real offset = MP.initsize[0];
        
        mt1 = link( new Microtub(Vecteur(-len/2, -offset, 0), VX, len) );
        mt2 = link( new Microtub(Vecteur(-len/2, +offset, 0), VX, len) );
     
        pinch1( mt1, 0.0, VZERO, MAGIC_TYPE);
        pinch1( mt1, 1.0, VZERO, MAGIC_TYPE);
        pinch1( mt2, 0.0, VZERO, MAGIC_TYPE);
        pinch1( mt2, 1.0, VZERO, MAGIC_TYPE);
        
        mt1 = link( new Microtub(Vecteur(-0.5-len/2, MP.boxsize[1]/2-offset, 0),  VX, len) );
        mt2 = link( new Microtub(Vecteur(+0.5+len/2, MP.boxsize[1]/2+offset, 0), -VX, len) );
        
        pinch1( mt1, 0.0, VZERO, MAGIC_TYPE);
        pinch1( mt1, 1.0, VZERO, MAGIC_TYPE);
        pinch1( mt2, 0.0, VZERO, MAGIC_TYPE);
        pinch1( mt2, 1.0, VZERO, MAGIC_TYPE);

        mt1 = link( new Microtub(Vecteur(-len/2, -MP.boxsize[1]/2-offset, 0),  VX, len) );
        mt2 = link( new Microtub(Vecteur(+len/2, -MP.boxsize[1]/2+offset, 0), -VX, len) );
        
        pinch1( mt1, 0.0, VZERO, MAGIC_TYPE);
        pinch1( mt1, 1.0, VZERO, MAGIC_TYPE);
        pinch2( mt2, 0.0, VZERO, MAGIC_TYPE);
        pinch2( mt2, 1.0, VZERO, MAGIC_TYPE);
        
      } break;

      case 6:
        break;
	
      case 7:                 //internal links between the ends
        mt  = link( new Microtub(place, dir, len) );
        
        co = link( new Complex(0, VZERO) );
        co->attachTo1( mt, 0, MT_MINUS_END);
        co->attachTo2( mt, 0, MT_PLUS_END);
        break;
        
      case 8:         
        break;
        
      case 9: 
        break;
      default:  
        break;
      }
    }

  //add motors of all types
  for (int type=0; type < MP.MAX; ++type)
    for(int ii=0; ii < MP.cxmax[type]; ++ii)
      link( new Complex(type) );  
}


//========================================================================
//==============     TEST 2     MP.initcode=2??        ===================
//========================================================================

void SIM::populateTest2()
{
}

//========================================================================
//==============     TEST 3     MP.initcode=3??        ===================
//========================================================================



void SIM::populateTest3()
  //provide some incorrect mt-configurations
{
  Vecteur d1, d2, dir, place, shift=Vecteur(0, MP.boxsize[4], 0);
  real len;
  Microtub * mt=0;
  Complex * co=0;
  int p, i;

  for( int m=0; m < MP.mtmax; ++m ) {
    dir   = Microtub::initDirection();
    len   = Microtub::initLength();
    place = Microtub::initPosition( Microtub::space, digit(MP.mtinit,1)) - dir*len/2;
    
    switch( digit(MP.mtinit, 1) ) {
      
      case 0:             //perturb the line perpendicularily
        mt = link( new Microtub( place, dir, len) );
        
        d1 = Vecteur(-dir.YY, dir.XX, 0 );
        d2 = -d1;
        for( p = 0; p < mt->nbPoints(); ++p)
          if ( p % 2 )
            mt->movePoint(p, d1);
          else
            mt->movePoint(p, d2);   
        break;


      case 1:             //perturb the line in the axis
        mt = link( new Microtub( place, dir, len) );
        
        d1 = dir * mt->rodLength() / 4.0 ;
        d2 = -d1;
        for( p= 0; p < mt->nbPoints(); ++p)
          if ( p % 2 )
            mt->movePoint(p, d1);
          else
            mt->movePoint(p, d2);   
        break;


      case 2:             //collapse the tubule to one point
        mt = link( new Microtub( place, dir, len) );
        for( p= 0; p < mt->nbPoints(); ++p) mt->setPoint( p, place );
        break;

      case 3:             //fold one point of the tubule
        mt = link( new Microtub( place, dir, len) );
        p = mt->nbSegments();
        if ( p > 2 )
          mt->setPoint( p-1, mt->whereP( p-2 ));
      break;

      case 4:
        mt  = link( new Microtub( place, dir, len) );
        co = link( new Complex(0, VZERO) );
        co->attachTo1( mt, 0, MT_MINUS_END);
        co->attachTo2( mt, 0, MT_PLUS_END);
        i = mt->nbPoints();
        for( p = i/2; p < i; ++p)
          mt->setPoint(p, mt->whereP( i - 1 - p ));
        break;

      case 5:     // a bunch of complexes bound twice to the same microtubule:
        mt  = link( new Microtub( place, dir, len) );
        co = link( new Complex(0, VZERO) );
        co->attachTo1( mt, len * RNG.preal(), MT_MINUS_END);
        co->attachTo2( mt, len * RNG.preal(), MT_MINUS_END);
        break;


      case 7:
        co = link( new Complex(1, VZERO) );

        mt = link( new Microtub(place, dir, len, MT_ORIGIN) );
        pinch1( mt, 0, shift, 0);
        co->attachTo1( mt, len/2, MT_MINUS_END);
        
        len=Microtub::initLength();
        dir=Microtub::initDirection();
        
        mt = link( new Microtub(place, dir, len, MT_ORIGIN) );
        pinch1( mt, 0, -shift, 0);
        co->attachTo2( mt, len/2, MT_MINUS_END);
        break;

      case 8:
        co = link( new Complex(1, VZERO) );
        
        mt = link( new Microtub(place, dir, len, MT_ORIGIN) );
        pinch1( mt, 0.0, shift, 0);
        pinch1( mt, MP.mtrodlength, shift, 0);
        co->attachTo1( mt, len/2, MT_MINUS_END);
        
        len=Microtub::initLength();
        dir=Microtub::initDirection();
        
        mt = link( new Microtub(place, dir, len, MT_ORIGIN) );
        pinch1( mt, 0.0, -shift, 0);
        pinch1( mt, MP.mtrodlength, -shift, 0);
        co->attachTo2( mt, len/2, MT_MINUS_END);
        break;

      case 9:
        co = link( new Complex(1, VZERO) );
        
        mt = link( new Microtub(place, dir, len) );
        pinch1( mt, 0.0, shift, 0);
        pinch1( mt, 1.0, shift, 0);
        co->attachTo1( mt, len/2, MT_MINUS_END);
        
        len=Microtub::initLength();
        dir=Microtub::initDirection();
        
        mt = link( new Microtub(place, dir, len, MT_ORIGIN) );
        pinch1( mt, 0.0, -shift, 0);
        pinch1( mt, 1.0, -shift, 0);
        co->attachTo2( mt, len/2, MT_MINUS_END);
        break;
	
      default:  
        break;
      }
    }
}





//========================================================================
//==============     TEST 4     MP.initcode=4??        ===================
//========================================================================



void SIM::populateTest4()
{
  Microtub ** mt = new Microtub*[2];
  Complex * co = 0;
  real len, ang;
  Vecteur dir, pos;
  int i;

  for( i=0; i < MP.mtmax; ++i ) {
    ang = PI*RNG.sreal();
    pos = Microtub::initPosition( Microtub::space, digit(MP.mtinit,1));
    len  = MP.mtrodlength;
    dir  = Vecteur( cos(ang), sin(ang), 0);
    
    mt[0] = link( new Microtub(pos,  dir, len, MT_ORIGIN) );
    mt[1] = link( new Microtub(pos, -dir, len, MT_ORIGIN) );
    
    co = link( new Complex(0, VZERO) );
    co->attachTo1( mt[0], len/2, MT_MINUS_END);
    co->attachTo2( mt[1], len/2, MT_MINUS_END);
    
  }

  for( i=0; i < MP.mtmax; ++i )
    {
      ang = PI*RNG.sreal();
      pos = Microtub::initPosition( Microtub::space, digit(MP.mtinit,1));
      len  = 2 * MP.mtrodlength;
      dir  = Vecteur( cos(ang), sin(ang), 0);

      mt[0] = link( new Microtub(pos,  dir, len, MT_ORIGIN) );
      mt[1] = link( new Microtub(pos, -dir, len, MT_ORIGIN) );

      co = link( new Complex(1, VZERO) );
      co->attachTo1( mt[0], len/2, MT_MINUS_END);
      co->attachTo2( mt[1], len/2, MT_MINUS_END);
    }

  for( i=0; i < 2*MP.mtmax; ++i )
    {
      ang = PI*RNG.sreal();
      pos = Microtub::initPosition( Microtub::space, digit(MP.mtinit,1));
      len  = MP.mtrodlength;
      dir  = Vecteur( cos(ang), sin(ang), 0);

      mt[0] = link( new Microtub(pos,  dir, len, MT_ORIGIN) );

      link( new Grafted(2, pos) ) -> attachTo( mt[0], len/2, MT_MINUS_END );
    }

  Vecteur shift=Vecteur((MP.mtinitlength[0] - MP.boxsize[4])/2.0, 0, 0);
  for( i=0; i < 2*MP.asmax; ++i )
    {
      pos  = Microtub::initPosition( Microtub::space, digit(MP.mtinit,1));
      len  = Microtub::initLength();
      dir  = VX;

      mt[0] = link( new Microtub(pos, dir, len, MT_ORIGIN) );

      link( new Grafted(0, pos+shift) ) -> attachTo( mt[0], 0.1, MT_PLUS_END );
      link( new Grafted(0, pos-shift) ) -> attachTo( mt[0], 0.1, MT_MINUS_END  );

    }

  delete mt;

}



//========================================================================
//==============     TEST 5     MP.initcode=5??        ===================
//========================================================================


// 2D array of grafted to visually testing the attachment
void SIM::populateTest5()
{
  Vecteur w;
  switch( digit(MP.mtinit, 1) ) {
    
    case 0:	

      for( int ii=0; ii < 200; ++ii)
	for( int jj=0; jj < 200; ++jj) {
	  w.set( 0.1*(ii-100), 0.1*(jj-100), 0);
	  link( new Grafted(0,  w) );
	}
    break;	
  }

  for( int ii=0; ii < MP.mtmax; ++ii )
    link( new Microtub(0) );
}
