//RCS: $Id: scalar_field.cc,v 2.3 2005/05/04 15:34:34 clausen Exp $
//-------------------------------lattice.cc----------------------------------

#include "scalar_field.h"
#include "exceptions.h"
#include "iowrapper.h"
#include "iomessages.h"
#include "sim_param.h"

//----------------------------------------------------------------------------
ScalarField::ScalarField()
{
  if( MP.boxshape == 1) {        //square
    //sizeI = sizeJ = sizeK = 2 * int( round( (real)MP.boxsize[0] / (real)MP.ficellsize ) );

    sizeI = int( round( (real)MP.boxsize[0] / (real)MP.ficellsize ) );
    sizeJ = sizeI;
    sizeK = sizeI;

  }
  else
  if( MP.boxshape == 2) {  //periodic
    sizeI = 2 * int( round( (real)MP.boxsize[0] / (real)MP.ficellsize ) );
    sizeJ = sizeK = 2 * int( round( (real)MP.boxsize[1] / (real)MP.ficellsize ) );
  } else{
  printf("Error, not implemented for any other shapes other than square");
  }

  if( DIM == 1 ) {
    field = new real[ sizeI ];
  } else if( DIM == 2 ) {
    field = new real[ sizeI ][ sizeJ ];
  } else {
    field = new real[ sizeI ][ sizeJ][ sizeK ];
  }

  if( !field ) {
    MSG.error("scalar_field.cc:ScalarField", "Unable to allocate memory for field.");
    exit(1);
  }

  //initialize the field values
  //setFieldValues(); //CHAITANYA
  //CHAITANYA:
  //setFieldValues() is now set in sim_initial::populateRandom()

}

//----------------------------------------------------------------------------
//Here values are set to the field
void ScalarField::setFieldValues(int initCode){
//init code 0= initialize to zero, 1= initialize to some standard gradient(under development)
  switch(initCode){
    case 0:
      if( DIM == 1 ) {
	      for( int i=0; i<sizeI; i++ ) {
	        field[ i ] = 0.0;
	      }
      } else if( DIM == 2 ) {
      	      for( int i=0; i<sizeI; i++ ) {
				  for ( int j=0; i<sizeJ; j++ ) {
							  field[ i ] = 0.0;
					  }
			  	  }
      } else if( DIM == 3 ) {
    	      for( int i=0; i<sizeI; i++ ) {
				  for ( int j=0; i<sizeJ; j++ ) {
					  for ( int k=0; i<sizeK; j++ ) {
						  field[ i ] = 0.0;
					  }
				  }
      	      }
      }
     break;
    //---STANDARD GRADIENT of step shape- halfway in space
    case 1:
       //init code 1= initialize to some standard gradient(half-half)
       //1D
       printf("@scalarfield standard step gradient \n");
       switch(DIM){
       case 1:
      	   for( int i=0; i<sizeI; i++ ) {
       			if(i<sizeI/2){
      			  field[ i ] = 0.0;
       			}
       			else if(i>=sizeI/2 && i<sizeI){
       		          field[ i ] = 1.0;
       	   		}
      	   }
       break;

       //2D
       case 2:
      	   for( int i=0; i<sizeI; i++ ) {
			   for( int j=0; j<sizeI; j++ ) {
				   if( i < (sizeI/2) ){
					   field[ i ][ j ] = 0.0;
				   }
				   else if( i >= (sizeI/2) ){
					   field[ i ] = 1.0;
				   }
			   }
		   }
       break;

       case 3:
       //3D
          for( int i=0; i<sizeI; i++ ) {
			  for ( int j=0; i<sizeJ; j++ ) {
					  for ( int k=0; k<sizeK; k++ ) {
						  if( i < (sizeI/2) ) {
							  field[ i ][ j ][ k ] = 0.0;
						  }
						  else if(i >= (sizeI/2) && i < sizeI && j < sizeJ && k < sizeK){
							  field[ i ] = 1.0;
						  }
					  }
				  }
       	       }
        break;
      } //SWITCH dimension 1,2,3
  break;
  ////////////


  ////////////


  } //SWITCH initcode
} //setValue


//------------Alternative setField with coordinates of object that acts as a source of the field

void ScalarField::setFieldValues(int initCode, Vecteur sourcePos ) {

	real sDist; //dist of field point from source
	Vecteur fieldV; //vector of a given position in the grid

	  switch( initCode ){


		  //---gradient dependent on position of some point object- eg. bead
		  case 2:

	      //init code 1= initialize to some standard gradient(half-half)
	      //1D
	      printf("@scalarfield exponential gradient\n");
	               switch(DIM){
			       case 1:
			      		   for( int i=0; i<sizeI; i++ ) {
							   fieldV[0] = (i / sizeI) ;
							   sDist = euclDist( sourcePos, fieldV );
			      		   }
			       	break;

			       //2D
			       case 2:
			      		   for( int i=0; i<sizeI; i++ ) {
							   for( int j=0; j<sizeJ; j++ ) {
								   printf("@scalar field setField vals \n");

							   }
						   }
			       	break;

			       	//3D
			       	case 3:

			        	  for( int i=0; i<sizeI; i++ ) {
							  for ( int j=0; i<sizeJ; j++ ) {
									  for ( int k=0; k<sizeK; k++ ) {

									  }
								  }
			       		       }
			        	break;
			      	} //SWITCH dimension 1,2,3

		  default:
		     MSG.error("Space::Position","unknowned init code");
    	 }//switch initCode


}//method
////////////

real euclDist(Vecteur a, Vecteur b){
	if( DIM == 1){
		real eDist = pow( ( pow( (a[0]-b[0]), 2) + pow( (a[1]-b[1]), 2) + pow( (a[2]-b[2]), 2) ), 0.5 );
	}
	else if( DIM == 2){
		real eDist = pow( ( pow( (a[0]-b[0]), 2) + pow( (a[1]-b[1]), 2) + pow( (a[2]-b[2]), 2) ), 0.5 );
	}
	else if( DIM == 3){
		real eDist = pow( ( pow( (a[0]-b[0]), 2) + pow( (a[1]-b[1]), 2) + pow( (a[2]-b[2]), 2) ), 0.5 );
	}
	return eDist;
}

//----------------------------------------------------------------------------
//---CHAITANYA: this converts the x,y,z spatial coordinates (including negative vals) into field-box units
/*
Google.sci.math
let's say you want the linear address of (a, b, c) where a goes
from 0 to m-1, b goes from 0 to n-1, and c goes from 0 to p-1. This
gives m, n, and p values of the coordinates, respectively, begining
with 0. Then the formula for the linear address of (a, b, c) is:

L = a*np + b*p + c.

You have mnp positions numbered from 0 to mnp -1. You can check the
formula because the max values are m-1, n-1, and p-1 for a,b, and c
respectively giving:

Lmax = (m-1)np + (n-1)p + (p-1)

= mnp - np +np - p + p -1 = mnp -1.

--Lynn
*/
real ScalarField::value(Vecteur position)
{
  int cellI; //location scalar field value to be accessed
  int cellJ; //location of scalar field
  int cellK; //location of scalar field
  real v; //measure of change of MT polymerization rate. without scale for now.

  switch(DIM){
  case 1:
  	cellI =  int( round( position[0] * sizeI / (real)MP.boxsize[0] ) );
  	//cellI =  int( round( ( position[0] + (real)MP.boxsize[0] ) / (real)MP.ficellsize ) );
  	if( cellI < 0 || cellI >= sizeI ) {
  	    printf("Error: outside scalar field boundaries\n");
	    v = 0.0;
	}
	else {
	    v = field[ cellI ];
    	}
    	break;
  case 2:
 	if(MP.ficellsize >0){
  		cellI =  int( round( position[0] * sizeI / (real)MP.boxsize[0] ) );
  		cellJ =  int( round( position[1] * sizeJ / (real)MP.boxsize[1] ) );

    	}else {
       		cellI = 0;
       		cellJ = 0;
       		cellK = 0;
       	}
       	//TEST reset to 0
       	//	cellI = 0;
	   	//  cellJ = 0;
       	//	cellK = 0;

    	if( cellI < 0 || cellI >= sizeI || cellJ < 0 || cellJ >= sizeJ ||cellJ < 0 ) {
	  	    printf("Error: outside scalar field boundaries\n");
		    v = 0.0;
	}
	else {
	    v = field[ cellI][ cellJ  ];//////////////////////
	}
    	break;
   case 3:
     	//cellK =  int( round( ( position[2] + (real)MP.boxsize[2] ) / (real)MP.ficellsize ) );
     	if(MP.ficellsize > 0){
     	cellI =  int( round( position[0] * sizeI / (real)MP.boxsize[0] ) );
       	cellJ =  int( round( position[1] * sizeJ / (real)MP.boxsize[1] ) );
       	cellK =  int( round( position[2] * sizeI / (real)MP.boxsize[2] ) );
       	}
       	else {
       		cellI = 0;
       		cellJ = 0;
       		cellK = 0;
       	}
       	printf("Error: not implemented yet!\n");

    	if( cellI < 0 || cellI >= sizeI || cellK >= sizeK ) {
   	  	    printf("Error: outside scalar field boundaries\n");
   		    v = 0.0;
   		    }
   		else {
			v = field[ cellI][ cellJ ][ cellK ];
			}
        //CHAITANYA: verify and print to make sure (16/02/2006)
    	break;
    }

  return v;
}


//-------------------------------------------------------------------------
//elementary steps for all the microtubules in 'this' box
void ScalarField::stepField()
{

	printf("The boxshape: \t%d\t size: \t%d position_0: \t%d\n", MP.boxshape, MP.boxsize);

}



