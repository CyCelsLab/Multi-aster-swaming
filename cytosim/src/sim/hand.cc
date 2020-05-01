//RCS: $Id: hand.cc,v 2.21 2005/04/22 18:54:52 nedelec Exp $

//-----------------------------------------------------------------------

//================================hand.cc================================

//-----------------------------------------------------------------------
// Neha Khetan, May  2016     		: Gliding Assay - hamodel case 5 added for kinesin like stepping behaior: v = vo*( 1 - f/fstall )
// Neha Khetan, May 2016      		: Gliding Assay - to mimic the scharrel model 2014 for kinesin with the similar detachment model considering the active n pasive motors
// Neha Khetan, June -july 2016     : Gliding Assay - Semi- disctrete continuous Stepping behaviour for discrete model when the dthastep > dt time.
// Neha Khetan, Sept 2016     		: Gliding Assay - hamodel, case 9 - to incorporate discrete stepping behaviour when dthastep > dt time.
// Neha Khetan, Sept 2016     		: Gliding Assay - getActionMinusNK addedunder hamodel, case 9 - to incorporate discrete stepping behaviour
// Neha Khetan, Sept  2016    		: Gliding Assay - Detachment increases exponentially with the load (Kramer's theory) - variable alpha or scaling factor added in the exponent: CASE 6 I.E. Instead of hrad-wiring the scaling factor, introduced a parameter 
// 								                      hadetachkramerscale
// Neha Khetan, Dec.  2016          : Gliding Assay - Modified the motor model case 6, Continuous model -> wherein the FV is applicable only in the case of opposing force while at assisting force v = v unloaded
// Neha Khetan, Sept  2016     		: Detachment model wherein, case 9
// Neha Khetan, May 2016 for Gliding Assay- for kinesin like detachment model -to mimic scharrel model & in general kinesin 5 like - wherein, stall force is different from the characteristic detachment force 
// Typically the motor detached before stalling - and that mean detachment force is termed as teh characteristic detachment force



#include "hand.h"

#include "grafted.h"

#include "complex.h"

#include "sim.h"

#include "exceptions.h"

#include "iowrapper.h"

#include <vector>
#include <iostream>

//------------------------------------------------------------------------

void Hand::handConstructorBase(const HandBaseType btype, Object * base, const int type)

{

  assert(( type >= 0 ) && ( type < MP.MAX ));

  haType            = type;

  haBaseType        = btype;

  haBase            = base;

  haStateUpdateTime = sim.iterationCnt();

  haBindingKey      = 0;

  hadtLastStep      = 0;



}



//------------------------------------------------------------------------

Hand::Hand(HandBaseType btype, Object * base, int type)

{

  assert(( 0 <= type ) && ( type < MP.MAX ));

  assert(( btype == 0 ) || ( base != 0 ));



  //if(( 0 > type ) || ( type >= MP.MAX ))

  //MSG.error("Hand::Hand","wrong type %i", type);



  handConstructorBase(btype, base, type);

  haState = 0;

}



//------------------------------------------------------------------------

Hand::~Hand()

{

  MSG(13, "Hand::destructor %p\n", this );

  if ( isAttached() ) detach( DETACH_DESTRUCT );

}





//------------------------------------------------------------------------

/**returns the distal position where the other hand is linked. 

If the Hand is grafted, the grafted position is returned.

If the Hand is in a complex, the position of the other Hand in the complex is returned

*/

Vecteur Hand::whereOtherSide() const

{

  switch( haBaseType ) {

    case HA_BASE_NONE:

      return VZERO;

      

    case HA_BASE_CX_SIDE1:

      //assert( dynamic_cast<Complex*>(haBase) );

      return static_cast<Complex*>(haBase)->cxHand2.whereHand();

      

    case HA_BASE_CX_SIDE2:

      //assert( dynamic_cast<Complex*>(haBase) );

      return static_cast<Complex*>(haBase)->cxHand1.whereHand();

      

    case HA_BASE_GRAFTED:

      //assert( dynamic_cast<Grafted*>(haBase) );

      return static_cast<Grafted*>(haBase)->whereGrafted();

      

    default:

      MSG.error("Hand::","unknown basetype in Hand::OtherWhere()");

      return VZERO;

    }

}



//------------------------------------------------------------------------

Vecteur Hand::getStretch() const

{

  Vecteur s = whereOtherSide() - whereHand(); 

  Space::modulo(s); 

  return s; 

}





//------------------------------------------------------------------------

//------------------------------------------------------------------------

//------------------------------------------------------------------------





int Hand::otherBound() const

{

  switch( haBaseType ) {

    case HA_BASE_CX_SIDE1: return static_cast<Complex*>(haBase)->cxHand2.isAttached();

    case HA_BASE_CX_SIDE2: return static_cast<Complex*>(haBase)->cxHand1.isAttached();

    case HA_BASE_GRAFTED:  return 1;

    default:

    case HA_BASE_NONE:     return 0;

    }

}





//------------------------------------------------------------------------

const Hand * Hand::otherHand() const

{

  switch( haBaseType ) {

    case HA_BASE_CX_SIDE1: 

      return &( static_cast<Complex*>(haBase)->cxHand2 );

    case HA_BASE_CX_SIDE2:

      return &( static_cast<Complex*>(haBase)->cxHand1 );

    default: return 0; 

  }

}





//------------------------------------------------------------------------

const Microtub * Hand::otherMicrotub() const

{

  switch( haBaseType ) {

    case HA_BASE_CX_SIDE1: 

      return static_cast<Complex*>(haBase)->cxHand2.getMT();

    case HA_BASE_CX_SIDE2: 

      return static_cast<Complex*>(haBase)->cxHand1.getMT();

    default: return 0;

  }

}





//------------------------------------------------------------------------

void Hand::setState()

{

  haState           = getState();

  haStateUpdateTime = sim.iterationCnt();

}



//------------------------------------------------------------------------

void Hand::updateState(const ReasonToUpdate why ) 

{

  int newstate = getState();

 

  if ( newstate != haState )

    {

      //MSG(9,"updateState H%lx : %i -> %i\n", this, haState, newstate );



      haState           = newstate;

      haStateUpdateTime = sim.iterationCnt();



      switch( haBaseType ) 	{

        case HA_BASE_CX_SIDE1:

        case HA_BASE_CX_SIDE2:

          static_cast<Complex*>(haBase) -> updateState( why );

          break;

        case HA_BASE_GRAFTED:

          static_cast<Grafted*>(haBase) -> updateState( why );

          break;

        case HA_BASE_NONE: break;

        default:

          MSG.error("Hand::updateState", "haBaseType not set");

      } 

    }

}



///true if the state was updated at this time step

bool Hand::isStateNew() const

{

  return haStateUpdateTime == sim.iterationCnt();

}





//------------------------------------------------------------------------

//                       ATTACH --- DETACH

//------------------------------------------------------------------------



bool Hand::bindingIsForbidden( const PointMicrotub & site ) const

{

#ifdef USE_LATTICE

  //check if lattice site is free for binding

  if ( ! site.getMT() -> siteFree ( site.abscissa() ) ) {

    return true;

  }

#endif



  //if Hand is part of a complex, we prevent attachements if that produces a useless link:

  if (( haBaseType == HA_BASE_CX_SIDE1 ) || ( haBaseType == HA_BASE_CX_SIDE2 )) {

    const PointMicrotub * otherSite = otherHand() -> getMTRod();

    if ( otherSite ) {

      

      //test if binding is directly adjacent to the place where the other Hand

      //of the complex is already attached. That would create a link internal to 

      //a microtubule, on the same rod or directly adjacent rods.

      //we compare the address of the rods in memory. This works because

      //getMTRod() returns pointers to objects stored consecutively in an array

      //the address differ by one only if they belong to the same MT, and adjacent

      if ( abs( &site - otherSite ) <= 1 ) 

        return true;

      

      //test if binding would make a link inside an aster, near the center:

      //i.e. a link between two microtubules from the same aster, very close to center

      if ( site.getMT() -> getOrganizer() == otherSite->getMT() -> getOrganizer() )

        if (( site.abscissa() < MP.assize[0] ) && ( otherSite->abscissa() < MP.assize[0] ))

          return true;

    }

  }

  

  //immediate end-detachment prohibits attachments directly at end:

  if ( site.getEnd() && MP.haenddetach[ haType ] ) {

    return true;

  }



  //if MP.haendattachdist is set and >0, this is a model of EB1

  if ( 0 < MP.haendattachdist[ haType ] ) {

    //FIX EXPERIMENTAL: special for Gohta simulations mars 2005. F. Nedelec

    // prevents EB1 to bind to the K-fibers (of type MAGIC_TYPE)

    if (( MP.initcode == 40 ) && ( site.getMT()->getType() == MAGIC_TYPE ))

      return true;

    

    //EB1 only binds to growing MT, near the plus-end:

    if ( MP.haendattachdist[ haType ] < site.abscissa(MT_PLUS_END) )

      return true;

    else

      return  ( site.getMT()->getDynamicState(MT_PLUS_END) == MT_SHRINKING );

  }



  

  //by default, binding is allowed:

  return false;

}



//------------------------------------------------------------------------

//this attach() is called by all the others !

//attachement is to the given MT site, to be interpreted as a InterpolatedPoint



bool Hand::attach( PointMicrotub & site )

{

  assert( isFree() );

  assert( isLinked() == 0 );

  //MSG(9,"attach to M%lx, rod %i, ab %.2f\n", site->mt->getName(), site->rod, site->abs );

  

  //the argument 'site' is taken to be an interpolated point.

  //by calling setFromInterpolated(), we set haEnd and haAbs

  //from the interpolated position (mPS, mPoint1, mCoef).

  site.setFromInterpolated();

  assert( site.looksWrong() == NO_ERROR );



  //test if the binding is possible, if so, return without binding

  if ( bindingIsForbidden( site ) )

    return false;

  

  //copy the site to set the position of the hand:

  *((PointMicrotub*)this) = site;

  

  //register the hand in the list of attached hands of the corresponding MT

  getMT() -> linkHand( this );



  //update lattice on MT

  updateLattice();



  //update the state of the hand

  updateState( ATTACH );

  return true;

}





//----------------------------------------------------------------------------

bool Hand::attachTo( const Microtub * mt, const real ab, const MTEnd from )

{

  static PointMicrotub site;

  mt->setInterpolation( & site, ab, from );

  return attach( site );

}



//----------------------------------------------------------------------------



/**detach 'this', the parameter is the "reason" for detaching, 

   used for debugging and to build transitions statistics.

   The meaning of the codes is in the enum ReasonToUpdate

*/

void Hand::detach(const ReasonToUpdate why)



{

  //MSG(9,"detach %lx %i\n", this, why);

  assert( isAttached() );

  

#ifdef USE_LATTICE

  //Update lattice

  freeFromLattice();

#endif



  //we calculate the position of the complex before detaching:

  if (( haBaseType == HA_BASE_CX_SIDE1 ) || ( haBaseType == HA_BASE_CX_SIDE2 ))

    static_cast<Complex*>(haBase) -> setPosition( whereHand() );



  //we calculate the position of the diffusing grafted before detaching:

  if (( haBaseType == HA_BASE_GRAFTED ) && ( static_cast<Grafted*>(haBase) -> isDiffusing() ))

    static_cast<Grafted*>(haBase) -> setPosition( whereHand() );

  

  pop();

  clear();



  updateState( why );

}



//------------------------------------------------------------------------

//------------------------------------------------------------------------

//------------------------------------------------------------------------



//called in stepUnderForce()

inline void limitstep(real & a, const real b)

{

  if ( b >= 0 )

    {

      if ( a > 2 * b )  a = 2 * b; 

      else if ( a < 0 ) a = 0;

    }

  else

    {

      if ( a < 2 * b )  a = 2 * b; 

      else if ( a > 0 ) a = 0;

    }

}



//------------------------------------------------------------------------

enum HandAction Hand::getAction( Vecteur force )

{

  real Pfor = 0.9; // probability to go forward under no load

  real Pback = 0.1; // probability to go backward if load > f2

  real f1 = 3.0; // force before prob to advance starts to decline

  real f2 = 5.0; // force where Pfor is zero and prob to go back is Pback

  real a = log(100.0*(Pfor/Pback + 1.0))/(f2-f1); // exponent

  real P;



  // force parallel to the tubule is used for movement along the tubule

  real force_parr = - force * dirMT();



  // force perpendicular to the tubule is used for detachment - NOT CALCULATED HERE YET

  // real force_perp = force * dirMT();



  enum HandAction action = STAY;



  //printf("Force: %f\n",force);

  // calculate the probability to go forward

  if ( force_parr < f1 ) {

    P = Pfor;

  } else if ( force_parr < f2 ) {

    P = (Pfor + Pback)/exp(a*(force_parr-f1)) - Pback;

  } else {

    P = -Pback;

  }



  //printf("%f %f\n",force,P);

  

  // find action from probability

  if ( P < 0.0 ) { // if probability is negative, we may go back

    if ( RNG.test( -P ))

      action = MINUSSTEP;

  } else {

    if ( RNG.test( P ))

      action = PLUSSTEP;

  }

  return action;

}


// ===================================================================================================================
/*
The stepping model for yeast dynein
Yeast minimal construct used by Kunalika Jain for Gliding Assays
Construct - 331 kDa Minimal construct Peterson et. al. 2006
In Gennerich et. al. 2007, force induced stepping behaviour of dynein was studied ( Full length dynein )
FUNCTION TO assign steps as per the forward prob , September 2016
  Using an overly simplified stepping model.
  1. As pairwise distance from MT gliding all show peak around 8 nm with few 4 nm and very few higher
     Thus fair to assume, in experiments we do not capture very fine spatial resolution.
  2. 
*/

real Hand::getActionMinusNK( Vecteur force ){
    // Based on Gennerich et al. (2007) model of dynein discrete steps in single molecule assays ( Over Simplified )   
    // Force cutoff values for different fwd and backwd stepping changed based on dyn331 construct reports and some 'guesswork'
  
}
// ===================================================================================================================


// ===================================================================================================================
/*
The stepping model for yeast FULL LENGTH dynein
Yeast : FULL LENGTH Peterson et. al. 2006
In Gennerich et. al. 2007, force induced stepping behaviour of dynein was studied ( Full length dynein )
FUNCTION TO assign steps as per the forward prob , September 2016
  Using an overly simplified stepping model.
  1. As pairwise distance from MT gliding all show peak around 8 nm with few 4 nm and very few higher
     Thus fair to assume, in experiments we do not capture very fine spatial resolution.
  2. 


real Hand::getActionMinusNK_WT( Vecteur force ){
    // Based on Gennerich et al. (2007) model of dynein discrete steps in single molecule assays ( Over Simplified )   
    // Force cutoff values for different fwd and backwd stepping changed based on dyn331 construct reports and some 'guesswork'
  
}
*/
// ===================================================================================================================

// ===================================================================================================================
/*
 here.....
*/


//------------------------------------------------------------------------

enum HandAction Hand::getActionMinus( Vecteur force )

{
	
	real Pfor = 0.9; // probability to go forward under no load
	
	real Pback = 0.1; // probability to go backward if load > f2
	
	real f1 = 3.0; // force before prob to advance starts to decline
	
	real f2 = 5.0; // for where Pfor is zero and prob to go back is Pback
	
	real a = log(100.0*(Pfor/Pback + 1.0))/(f2-f1); // exponent
	
	real P;
	
	
	
	// force parallel to the tubule is used for movement along the tubule
	
	//real force_parr = - force * dirMT(); //CA: multiplication of scalars to get magnitude
	real force_parr = force * dirMT();

	
	
	// force perpendicular to the tubule is used for detachment - NOT CALCULATED HERE YET
	
	// real force_perp = force * dirMT();
	
	
	
	enum HandAction action = STAY;
	
	
	
	//printf("Force: %f\n",force);
	
	// calculate the probability to go forward
	
	if ( force_parr < f1 ) {
		
		P = Pfor;
		
	} else if ( force_parr < f2 ) {
		
		P = (Pfor + Pback)/exp(a*(force_parr-f1)) - Pback;
		
	} else {
		
		P = -Pback;
		
	}
	
	
	
	//printf("%f %f\n",force,P);
	
	
	
	// find action from probability
	
	if ( P < 0.0 ) { // if probability is negative, we may go back
		
		if ( RNG.test( -P ))
			
			action = PLUSSTEP;
		
	} else {
		
		if ( RNG.test( P ))
			
			action = MINUSSTEP;
		
	}
	
	return action;
	
}



//------------------------------------------------------------------------

void Hand::stepUnderForce( const Vecteur & stretch )

{

  assert( isAttached() );

  assert( looksWrong() == NO_ERROR );

  

  //the load is the projection of the stretch on the local direction of tubule

  real load = stretch * dirMT();



  //we calculate the detachement of hands, according to different models:

  switch( MP.hadetachmodel[haType] ) 
      {
    
      //detachment occurs immediately if the load is above MP.force
      case 1: 
        
        if ( stretch.normSquare() > MP.habreakstretch_sq[haType] )
          { detach( DETACH_MAXFORCE );  return; }
          //detachment occurs also spontaneously: it goes to case 0: from here


       //detachment is independent of load
      case 0: {
          
        // If case is 0: also from case 1: for spontaneous detachment
        //there are two different detach. rates, depending if the motor is at end of MT

        real base_rate = ( getEnd() ? MP.haenddetachrate_dt[haType] : MP.hadetachrate_dt[haType] );
          if ( RNG.test( base_rate ))  
            { detach( DETACH_SPONTAN ); return; }
        } break;

      
      //detachment occurs if slightly pulled forward        
      case 2: {
  			if ( load * MP.haspeed[haType] > 0 )

  				{ detach( DETACH_MAXFORCE ); return; }

        	}break;

      case 3:{
  		  
  		    /* detachment increases exponentially with strech 	(based on kramer's theory)  
  		     Leduc et al. (2010) PRL: detachment rate (w_off) for a given strech value (y) is	  
  		                               w_off(y) = w_0 * exp( k*|y|/f_c )
  		                              where:
                                    f_c = characteristic detachment force, 
                                    w_0 = detachment rate in the absence of load, 
  		                              k = motor stiffness, 
                                    y = linker extension (stretch)
  		    */
  	     
          real base_rate = ( getEnd() ? MP.haenddetachrate_dt[haType] : MP.hadetachrate_dt[haType] );
        
          /* Neha Khetan, modifying for kinesin model , April 2017 commening the Lines below
            REPLACING THE STALL FORCE WITH A DETACHMENT FORCE OR CHARACTERISTIC FORCE
            Introduced a new parameter "charforce" 
            Now enables us to distinguis between the stall force as used in stepping models with 
            the detachment or characteristic force

  		      base_rate *= exp( MP.km * stretch.norm() / MP.haforce[haType] );
  		      printf("%9.4f\t%9.4f\t%2.2f\n", base_rate, stretch.norm() , MP.haforce[haType]  );
  		    */
  		 
  		    base_rate *= exp( MP.km * stretch.norm() /  MP.charforce[haType] );
  		    if ( RNG.test( base_rate ))
        		  { detach( DETACH_SPONTAN ); return; }
  	  
  	      } break;

      case 4:{

        /* ================================================================================================================
        CAUTION: Neha Khetan, In this model the description suggests that the scaling factor in Kramer's model is 2.
        However, the code implements it as the "MP.hadetachmodel[haType]" which suggests the value is 4.
        This NOW explains why the data.in file has a haforce value of 7 pN, as 4*F/7 which results in an effective stall force
        value of 1.75 pN.
        While USING THIS CONSIDER ----!!!!
        ================================================================================================================
        */


        // Used by CA , for Phys. Bio. 2014 paper

        // detachment increases exponentially with the load (Kramer's theory)

        // the detachment rate is multiplied by exp( 2 * force / max_force )

        // the factor 2 thus gives a off-rate amplification by ~10, when force ~ max_force 
  	
        real base_rate = ( getEnd() ? MP.haenddetachrate_dt[haType] : MP.hadetachrate_dt[haType] );

        // ----- NEHA KHETAN, Understanding n Checking the model
        // printf("%i : %9.4f \n", getEnd(), base_rate);
        // printf("basal dettach rate:%2.4f \n", base_rate);
        // base_rate *= exp( stretch.norm() / MP.hamaxstretch[haType] ); 
  	    // base_rate *= exp( 4 * stretch.norm() / MP.hamaxstretch[haType] ); 
  	  
        base_rate *= exp( MP.hadetachmodel[haType] * stretch.norm() / MP.hamaxstretch[haType] ); 

        //printf("%9.4f\t%d\t%9.4f\t%9.4f\t%2.2f\n", base_rate, MP.hadetachmodel[haType] ,stretch.norm() , MP.hamaxstretch[haType]  , MP.km);
        //printf(" %9.4f\n", base_rate);

        if ( RNG.test( base_rate ))

          { detach( DETACH_SPONTAN ); return; }

        } break;

      case 5: {
  		    /* 20 Feburary 2017, Neha Khetan, 
  		    For Gliding Assay project,
  		    Implementing the Anisotropic detachment of dynein molecules
          HARD-WIRED MODEL:
          USE ONLY FOR TRUNCATED MINIMAL YEAST DYNEIN
  		    
          Based on the experimental measurements from Nicholas et. al. 2015,
          Figure 1G: was digitized to extract coordinates and from data fitting a functional form was obtained
          */

  			real base_rate = 0;
  			real loadforce = MP.km*load;

  			//printf("Detach:%0.6lf\n", loadforce) ;\

  			//THIS IS BACKWARD LOAD- detachment acc. to one single expressio that has saturation behaviour
  			if( loadforce > 0 )
  			{
  				//EQUATION OPTIMIZED TO DATA FROM NICHOLAS ET AL. 2015 
  				base_rate = ( 2.3159 + ( 0.04 - 2.3159 ) *exp( - 0.8414 * loadforce) ) * MP.dt;

  			//THIS IS FORWARD LOAD- detachment goes through 1 threshold below which it's saturation, above which linear
  			}else if ( loadforce < 0 ){ 
  					//THRESHOLD FORCE ABOVE WHICH (SIGN) it follows saturation kinetics			
  					if( std::abs(loadforce) <= 3.6225 ){
  		 				base_rate = ( 7.2768 + ( 0.04 - 7.2768 ) *exp( - 0.6426 * std::abs(loadforce) ) ) * MP.dt;
                     //ALTERNATIVELY: linear 
  					}else{
  						base_rate = ( ( 3.5814* std::abs(loadforce) ) -6.4016 ) *MP.dt;
  					}

  			} else{

  				base_rate = MP.hadetachrate_dt[haType];	
  			}

  		    if ( RNG.test( base_rate ))

  				{ 
  				detach( DETACH_SPONTAN ); 
  				//printf("Detach\t%.6f\t%.6f\n", loadforce , base_rate);		

  				return; }
  	  
  	   } break;
  		  

  	 case 6: 
  		  {
        /* Modified by Neha Khetan, 5 September 2016: 
           Instead of hard-wiring the scaling factor, 
           Introduced a parameter hadetachkramerscale

           Detachment increases exponentially with the load (Kramer's theory)
           Scaling factor is multiplied in the exponent to scale the detachment rate , 
           rd = rb*exp( scaling factor * force / max_force )
           BASED ON STALL FORCE
    
        */      
    		  // printf("Hi! I am in hand detach model 6\n");
    		  // printf("Detachmodel:case 6 \n");
    		  // Modified by Neha Khetan, 5 September 2016
    		  
    		  real base_rate = ( getEnd() ? MP.haenddetachrate_dt[haType] : MP.hadetachrate_dt[haType] );
    		  base_rate *= exp( MP.hadetachkramerscale * stretch.norm() / MP.hamaxstretch[haType] ); 		  
    		  
    		  
    		  if ( RNG.test( base_rate ))
    		    { 
              detach( DETACH_SPONTAN ); 
              //printf("%2.3f \t %2.6f \t %2.6f\n", sim.simTime(), stretch.norm() , base_rate);
              return; }
          } break;



      case 7:{
  		  /* Neha Khetan, September- October 2016
           TRYING out different detachment model behaviors
           to be able to explain the Gliding assay behavior of dyneins

           Modified further based on Rose Loughlin Model
           Immediate detachment at the end,

        */


        if ( getEnd() )
           { detach( DETACH_SPONTAN ); return; }

        else 
  		    {
          real base_rate =  MP.hadetachrate_dt[haType];		 		  
  		    base_rate *= exp( MP.km * stretch.norm() / MP.charforce[haType] );        // Based on detachment force
  		    //base_rate *= exp( 2*MP.km * stretch.norm() / MP.haforce[haType] );      // Rose Loughlin model
            if ( RNG.test( base_rate ))
              { detach( DETACH_SPONTAN ); return; }
              //printf("%i : %9.4f \t %9.4f \n", getEnd(), base_rate , MP.charforce[haType] ) ;
  		    }
        }break;

      
      case 8: {
        /* Neha Khetan, 23 Dec. 2016,
           Attempt for Gliding Assay
           MP.charforce[haType] corresponds the the force threshold 
           Attempt to implement a sort of piece-wise model
           For a threshold the DETACH RATE IS A CONSTANT VALUE WHILE LESS THAN OR EQUAL UNDERGOES
           STOCHASTIC detachment based on the load free detachment value

           Attempt to optimize and scan few values
        */
  			real base_rate = 0;
  			real detachforce = MP.km*stretch.norm();
  			
  	
  			if ( detachforce > MP.charforce[haType]  ){

  				  base_rate =  2*MP.dt;

  			 }else{	
            base_rate = MP.hadetachrate_dt[haType];	
  			 }	
  		      
        if ( RNG.test( base_rate ))
            { detach( DETACH_SPONTAN );   return; }
        	  //printf("%0.9lf \t %0.6f \n",  base_rate , detachforce  ) ;

        }break;


  	 case 9: {
           /* Neha Khetan, 23 Dec. 2016,
              Attempt for Gliding Assay
              if detach force is less than the ' characteristic detach force' 
                    follow a kramer law as f( stall force )
              else 
                    follow kramer law as f( charforce )

              And, finally detachment is stochastic

              Motivated by Deepak Bhat and Gopalkrishnan Phys. Bio. 2012

           */

  			real detachforce =  MP.km*stretch.norm();
  			real base_rate   =  MP.hadetachrate_dt[haType];		 		  
  			
  			 if ( detachforce < MP.charforce[haType]  ){
              base_rate *= exp( MP.hadetachkramerscale * MP.km*stretch.norm() / MP.haforce[haType] ); 

  			 }else{	
              base_rate *= exp( MP.hadetachkramerscale * MP.km*stretch.norm() / MP.charforce[haType] ); 
  			 }	
  		    
         if ( RNG.test( base_rate ))
            { detach( DETACH_SPONTAN ); return; }
          
          //printf("%0.9lf \t %0.6f \n",  base_rate , detachforce  ) ;
        }break;
  	


    default:
        MSG.error("Hand::Step","unknowned hadetachmodel for type <%i>",haType);
        break;
    }

  
  // ====================================  STEPPING OF THE HAND =====================


  // we now calculate the dabslacement of the hand on the filament

  real dabs=0;


  enum HandAction action;
  // CA 2013-02-13
  // hamodel decides what kind of movement, 0, 1, 2, 3 are continuous models, 
  // 4 is a discrete model (depending on dthastep and mtdimersize)


  switch ( MP.hamodel[haType] ) 
    {

      case 0:                     //the motor cannot go backward

        dabs = MP.haspeed_dt[haType] + load * MP.havarspeed_dt[haType];

        limitstep( dabs, MP.haspeed_dt[haType] ); 

        break;



      case 1:                     //the motor can go backward under a reverse load

        dabs = MP.haspeed_dt[haType] + load * MP.havarspeed_dt[haType];

        break;



      case 2:                     //the motor speed is constant

        dabs = MP.haspeed_dt[haType];

        break;

      

      case 3:                     //the motor never keeps any stretch

        dabs = 0.5 * load;

        break;


      case 4:

	         //printf("Hi! I am in hand stepping model 4\n");
          hadtLastStep += MP.dt;
		      //printf("motor type: %i\n", haType);
           //Chaitanya- taking care of plus and minus-ended motion 2013-02-25
           while ( hadtLastStep > MP.dthastep ) 
            {
                  //if( haType == 0 ){ FN: Is very dangerous: you should check the sign of haspeed[], 2013/02/25
                  //if ( MP.haspeed[haType] >=0) {
                  if ( MP.haspeed[haType] >0) {
                    //plus-ended motor
                    action = getAction( stretch*MP.km ); // give force as argument
                    //}else if (haType == 1) { 2013/02/25
                    //}else if( MP.haspeed[haType]  < 0 ){ //indicative of negative velocity
                  }else{ 
                    //Minus-ended motor
                    action = getActionMinus( stretch*MP.km ); // give force as argument
                  }
                  //else {
                  //  printf("Hand type for this method not implemented");
                  //}


                  //real dabs = 0;  
                  switch ( action ) { 
                    case PLUSSTEP:
                        if ( getMT() -> siteFree ( haAbs + dabs + MP.mtdimersize ) )
        			             moveByCheckLattice(MP.mtdimersize);
                          //dabs += MP.mtdimersize;
                          break;

                    case MINUSSTEP:
                        if ( getMT() -> siteFree ( haAbs + dabs - MP.mtdimersize ) )
              			       moveByCheckLattice(-MP.mtdimersize);
                          //dabs -= MP.mtdimersize;
                         break;

                    case STAY:
                        break;

                    case DETACH: // is not used for now
           	            exit(1);
                        detach( DETACH_SPONTAN );
    		            return;
                    break;

                    default:
                        MSG.error("Hand::Step","unknowned HandAction");
                        break;

                    }
              hadtLastStep -= MP.dthastep;
            }     
	        return;

       
    case 5: 
            // Neha Khetan, 5 May 2016, case 5 an case 6 are inspiried from Scharrel et. al. 2015, 
            // multi-motor transport - active and inactive kinesins         
        		//printf("Hi! I am in hand stepping model 5\n");     
        		// Typical kinesin like model for stepping.		
            // SAME AS CASE 2

            // Neha Khetan, Modifying Model 5 with stepping rate as old 5 is same as case 2
            // Addind stepping rate
            // In attempt to understand the Gliding assays ( of Kunalika Jain)   
            // This is a typical Kinesin like model of stepping  - Corrected version of case 5
            // Stepping rate has been incorporated to account for the steeping
            // Last seen Jan 7 2017
            // WORKS FOR KINESIN GLIDING ASSAY ONLY
                
            // As Case 5, Incorporated a stepping rate criterion

            // Variable stepping rtae
           {
            
            }break;
            		


  	case 6:  
          {
  		     /* Neha Khetan, Modfying: Dec 2016
  		        Linear FV curve is applicable during the opposing loads , 
              under assiting loads - velocity is same as the unloaded velocity
           */

      		if ( MP.haspeed[haType] < 0) {
      		 		if ( load > 0.00 ){
      					dabs = ( MP.haspeed[haType] * ( 1 -  ( abs(MP.km * load ) /MP.haforce[haType] )) )* MP.dt ;
      				}else	        
      					dabs = MP.haspeed[haType] * MP.dt; 
      	
      			
      		}else{
      		 		if ( load < 0.00 ){
      					dabs = ( MP.haspeed[haType] * ( 1 -  ( abs(MP.km * load ) /MP.haforce[haType] )) )* MP.dt ;
      				}else	        
      					dabs = MP.haspeed[haType] * MP.dt; 
      			}

          	}break ;



	
  
 case 7:
          {
		      dabs = 0.0 ;
          //printf("%2.6f\n" , dabs);
          }break;

	
 case 8:
    // Applicable Only for Truncated Yeast dynein
 		// Neha Khetan, July 2016 - Stepping behaviour for discrete model with a STEPPING RATE
		// If dthaspeed > dt: 	
		// TO DO: Implement in case of '+' end motor
    // GLIDING ASSAY
		{

			real steprate      = 0.00;
			real Steploadforce = 0.00;
			Steploadforce      = ( MP.km * load ) ;	
			

			 // If + end directed motor
			 if ( MP.haspeed[haType] >0) {
			     printf("ERROR: Not implemented !!!!!!!!!!!!! \n");

			 // If - end directed motor
			 }else if ( MP.haspeed[haType] < 0) {

					// Check for stepping rate criterion
					      
					steprate  =  8* MP.dt ;      //( std::abs(MP.haspeed[haType])/(std::abs(dabs)) ) * MP.dt;  		
		
					real step_prob = 0.00;
					step_prob      = 1 - exp( - steprate ) ;			
				
					if (  RNG.preal() > step_prob  ){    // i.e. of RN <= steprate; then the hand steps else
						dabs = 0.00;
					 
					}else{
						dabs = getActionMinusNK( stretch*MP.km  );
						//printf("Step\t%2.6lf\t%2.6lf\n", Steploadforce, dabs);
						}
			// No motion
			}else{
				dabs  = 0.00;
		 	};
		//printf("Step:%2.6lf\n", Steploadforce);
		//printf("Time: %f \t F:%f \t dabs: %2.6lf\n",  haStateUpdateTime , force_parr, dabs);
		}break;

	

	case 9:
 		// Neha Khetan, July 2016 - Stepping behaviour for discrete model when the dthastep > dt time.
		// If dthaspeed > dt:
	
 	
		// TO DO: Implement in case of '+' end motor
		{
			
			 // If + end directed motor
			 if ( MP.haspeed[haType] >0) {

				printf("Not implemented\n");

			 // If - end directed motor
			 }else if ( MP.haspeed[haType] < 0) {
		
				dabs = getActionMinusNK( stretch*MP.km  );

			// No motion
			}else{
				dabs  = 0.00;
		 	};

			
			// Check for stepping rate criterion
			real steprate = 0.00;      
			steprate  =  ( std::abs(MP.haspeed[haType])/(std::abs(dabs)) ) * MP.dt;  		
		
			real step_prob = 0.00;
			step_prob = 1 - exp( - steprate ) ;			
				
			if (  RNG.preal() > step_prob  ){    // i.e. of RN <= steprate; then the hand steps else
				dabs = 0.00;
			 }
		//printf("!step%f\n", dabs);
		//printf("Time: %f \t F:%f \t dabs: %f\n",  haStateUpdateTime , force_parr, dabs);
		}break;

	

	case 10:
 		// Neha Khetan, July 2016 - Stepping behaviour for discrete model when the dthastep > dt time.
		// If dthaspeed > dt:
		// Stepping rate independent of the detachment module incorporated
		// Last seen 7 Jan 2017
	
 	
		// TO DO: Implement in case of '+' end motor
		{
			
			 // If + end directed motor
			 if ( MP.haspeed[haType] >0) {

				printf("Not implemented\n");

			 // If - end directed motor
			 }else if ( MP.haspeed[haType] < 0) {
		
				dabs = getActionMinusNK( stretch*MP.km  );

			// No motion
			}else{
				dabs  = 0.00;
		 	};

			
			// Check for stepping rate criterion
			real steprate = 0.00;      
			steprate  =  ( std::abs(MP.haspeed[haType])/(std::abs(dabs)) ) * MP.dt;  		
		
			real step_prob = 0.00;
			step_prob = 1 - exp( - steprate ) ;			
				
			if (  RNG.preal() > step_prob  ){    // i.e. of RN <= steprate; then the hand steps else
				dabs = 0.00;
			 }

		//printf("!step%f\n", dabs);
		//printf("Time: %f \t F:%f \t dabs: %f\n",  haStateUpdateTime , force_parr, dabs);
		}break;



    default:

      MSG.error("Hand::Step","unknowned hamodel for type <%i>", haType);

      break;

    }



  //if ( fabs(dabs) > 0.1 )

  //MSG(3, "Hand::step( load %8.3f ) : large step %8.3f um\n", load, dabs );

#ifdef USE_LATTICE

  moveByCheckLattice( dabs );

#else

  moveBy( dabs );

#endif

}





//------------------------------------------------------------------------

//STEP in the abscence of load:

void Hand::step()

{

  assert( looksWrong() == NO_ERROR );

  assert( isAttached() );

  //assert( otherBound() == 0 );



  //there are two different detach. rates, depending if the motor is at end of MT

  real base_rate = getEnd() ? MP.haenddetachrate_dt[haType] : MP.hadetachrate_dt[haType];

  if ( RNG.test( base_rate )) {

    detach( DETACH_SPONTAN ); 

    return; 

  }

  enum HandAction action;

  if ( MP.hamodel[haType] == 4 ) {

    hadtLastStep += MP.dt;

    while ( hadtLastStep > MP.dthastep ) {
		
		action = getAction( VZERO );

      //if ( getAction( VZERO ) == PLUSSTEP ) {
	  if ( action == PLUSSTEP ) {

        moveByCheckLattice( MP.mtdimersize );

      }

      //if ( getAction( VZERO ) == MINUSSTEP ) {
	  if ( action == MINUSSTEP ) {

        moveByCheckLattice( -MP.mtdimersize );

      }

      hadtLastStep -= MP.dthastep;

    }

	// Neha, 28 September 2016

	}else if ( MP.hamodel[haType] == 8 ){

		real P;
		P = RNG.preal();		
		real stepsize = 0.008;          // 8 nm
	    real Pfor = 0.80 ;              // Assuming at load free and assisting load condition the stepping behavior is similar to no load condition // assuming based on Gennerich 2007
		real dabs = 0.00 ;
			if ( P <= Pfor )
				dabs = -stepsize;      // ensures forward stepping implies movement towards the minus end
			else
				dabs = stepsize;
		
  } else {

    //motion is smooth, at the maximum speed:

#ifdef USE_LATTICE

    moveByCheckLattice ( MP.haspeed_dt[haType] );

#else

    moveBy ( MP.haspeed_dt[haType] );

#endif

  }

}



//------------------------------------------------------------------------

void Hand::moveByCheckLattice( real deltaAbs )

{

  if ( getMT() -> siteFree ( haAbs + deltaAbs ) ) {

    moveBy( deltaAbs );

  }

}



//------------------------------------------------------------------------

//-----------------------------------FILE---------------------------------

//------------------------------------------------------------------------



//------------------------------------------------------------------------

int Hand::looksWrong() const

{  

  if (( haType < 0 ) && ( haType >= MP.MAX )) return 1;



  if ( isFree() ) return NO_ERROR;



  return PointMicrotub::looksWrong();

}



//------------------------------------------------------------------------

/** hands are printed with the name of the microtubules they bind, 

an abscisse on this microtubule, and an indication of the 'end' state

*/

void Hand::write()

{

  IO.writeUInt8( haType );

  if ( isAttached() ) {

    checkBoundaries();

    assert( looksWrong() == NO_ERROR );

    IO.writeUInt16( getMT()->getName() );

    IO.writeReal32( haAbs );

    IO.writeUInt8( haEnd );

  }

  else

    IO.writeUInt16( 0 );

}







//------------------------------------------------------------------------

void Hand::read()

{

  //do not reset since motors might have already set the hand states...

  //that would screw up the lists



  haType = IO.readUInt8();

  

  if (( haType < 0 ) || ( haType >= MP.MAX ))

    throw IOException("wrong haType in Hand::read");

  

  Microtub * haMT = sim.readMicrotubName();



  MSG(80, "Hand::read()\n");

  

  if ( haMT == 0 ) {

  

    //the Hand is now detach, we need to deregister from any previously attached MT

    if ( isLinked() ) pop();

    clear();

    

  } else {

    

    //redo the link if the microtubule has changed:

    if ( haMT != getMT() ) {

      if ( isLinked() ) pop();

      mPS = haMT;

      haMT->linkHand( this );

    }



    haAbs = IO.readReal32();

    haEnd = (MTEnd) IO.readUInt8();

    

    //check the validity of haAbs:

    checkBoundaries();  

    updateInterpolated();

  }

  updateState( DETACH_FILEREAD );

}

