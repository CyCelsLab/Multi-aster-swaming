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
    //sizeX = sizeY = sizeZ = 2 * int( round( (real)MP.boxsize[0] / (real)MP.ficellsize ) );

    sizeX = int( round( (real)MP.boxsize[0] / (real)MP.ficellsize ) );
    sizeY = sizeX;
    sizeZ = sizeX;
  }
  else
  if( MP.boxshape == 2) {  //periodic
    sizeX = 2 * int( round( (real)MP.boxsize[0] / (real)MP.ficellsize ) );
    sizeY = sizeZ = 2 * int( round( (real)MP.boxsize[1] / (real)MP.ficellsize ) );
  } else{
  printf("Error, not implemented for any other shapes other than square");
  }

  if( DIM == 1 ) {
    field = new real[ sizeX ];
    //real fmatrix[sizeX];

  } else if( DIM == 2 ) {
    field = new real[ sizeX * sizeY ];
    //real fmatrix[sizeX][sizeY];

  } else {
    field = new real[ sizeX * sizeY * sizeZ ];
    //real fmatrix[sizeX][sizeY][sizeZ];
  }

  if( !field || !fmatrix ) {
    MSG.error("scalar_field.cc:ScalarField", "Unable to allocate memory for field.");
    exit(1);
  }


  //real matrix1[sizeX][sizeY][sizeZ];

  //initialize the field values
  setFieldValues();
  //setFmatrix();

}

//----------------------------------------------------------------------------
//Here values are set to the field
void ScalarField::setFieldValues(int initCode){
//init code 0= initialize to zero, 1= initialize to some standard gradient(under development)
  switch(initCode){
    case 0:
      if( DIM == 1 ) {
	      for( int i=0; i<sizeX; i++ ) {
	        field[ i ] = 0.0;
	      }
      } else if( DIM == 2 ) {
      	      for( int i=0; i<sizeX*sizeY; i++ ) {
      		  field[ i ] = 0.0;
      	      }
      } else if( DIM == 3 ) {
    	      for( int i=0; i<sizeX*sizeY*sizeZ; i++ ) {
        	field[ i ] = 0.0;
      	      }
      }
     break;
    case 1:
  //init code 1= initialize to some standard gradient(half-half)
    //1D
    switch(DIM){
      case 1:
      	   for( int i=0; i<sizeX; i++ ) {
       			if(i<sizeX/2){
      			  field[ i ] = 0.0;
       			}
       			else if(i>=sizeX/2 && i<sizeX){
       		          field[ i ] = 1.0;
       	   		}
      	   }
      break;

      //2D
      case 2:
      	   for( int i=0; i<sizeX*sizeY; i++ ) {
           		if(i<(sizeX*sizeY)/2){
      			field[ i ] = 0.0;
      	   	}
      	   	else if(i>=(sizeX*sizeY)/2 && i <sizeX*sizeY){
        		field[ i ] = 1.0;
       	   	}
       	   }
      break;

      case 3:
      //3D
          for( int i=0; i<sizeX*sizeY*sizeZ; i++ ) {
          	if(i<(sizeX*sizeY*sizeZ)/2){
   			field[ i ] = 0.0;
	  	}
	  	else if(i>=(sizeX*sizeY*sizeZ)/2 && i <sizeX*sizeY*sizeZ){
	        	field[ i ] = 1.0;
   		}
       	  }
       break;
      } //SWITCH dimension 1,2,3
  break;
  } //SWITCH initcode
} //setValue



//////---------
/*void ScalarField::setFmatrix(int initCode){
	//init code 0= initialize to zero, 1= initialize to some standard gradient(under development)
	  switch(initCode){
	    case 0:
	    	if( DIM == 1 ) {
			    for( int i=0; i<sizeX; i++ ) {
			      fmatrix[ i ] = 0.00;
			    }
	    	} else if( DIM == 2 ) {
	    	      for( int i=0; i<sizeX; i++ ) {
					  for( int j=0; j<sizeY; j++ ) {
						  fmatrix[ i ][ j ] = 0.00;
					  }
				  }
	    	} else if( DIM == 3 ) {
	    	      for( int i=0; i<sizeX; i++ ) {
					  for( int j=0; j<sizeY; j++ ) {
						  for( int k=0; k<sizeZ; k++ ) {
							  fmatrix[ i ][ j ][ k ] = 0.00;
						  }
					  }
	    		  }
			 }
	     	   break;
	    case 1:
	  //init code 1= initialize to some standard gradient(half-half)
			if( DIM == 1 ) {
			    for( int i=0; i<sizeX; i++ ) {
					if( i < sizeX/2)
					   fmatrix[ i ] = 0.00;
					else
					   fmatrix[ i ]= 1.00;
			    }
	    	} else if( DIM == 2 ) {
	    	      for( int i=0; i<sizeX; i++ ) {
					  for( int j=0; j<sizeY; j++ ) {
						  if( i < sizeX/2)
						     fmatrix[ i ][ j ] = 0.00;
						  else
						     fmatrix[ i ][ j ] = 1.00;
					  }
				  }
	    	 } else if( DIM == 3 ) {
	    	      for( int i=0; i<sizeX; i++ ) {
					  for( int j=0; j<sizeY; j++ ) {
						  for( int k=0; k<sizeZ; k++ ) {
							  if( i < sizeX/2)
							  	fmatrix[ i ][ j ][ k ] = 0.00;
							  else
							    fmatrix[ i ][ j ][ k ] = 1.00;
						  }
					  }
	    		  }
			  }

	  break;
  } //SWITCH initcode



}
*/
/////-----------



//----------------------------------------------------------------------------
real ScalarField::value(Vecteur position){
  int cellX; //location scalar field value to be accessed
  int cellY; //location of scalar field
  int cellZ; //location of scalar field
  real v; //measure of change of MT polymerization rate. without scale for now.

  switch(DIM){
  case 1:
  	cellX =  int( round( position[0] * sizeX / (real)MP.boxsize[0] ) );
  	//cellX =  int( round( ( position[0] + (real)MP.boxsize[0] ) / (real)MP.ficellsize ) );
  	if( cellX < 0 || cellX >= sizeX ) {
  	    printf("Error: outside scalar field boundaries\n");
	    v = 0.0;
	}
	else {
	    v = field[ cellX ];
    	}
    	break;
  case 2:
 	if(MP.ficellsize >0){
  		cellX =  int( round( position[0] * sizeX / (real)MP.boxsize[0] ) );
  		cellY =  int( round( position[1] * sizeY / (real)MP.boxsize[1] ) );

    	}else {
       		cellX = 0;
       		cellY = 0;
       		cellZ = 0;
       	}
       	//TEST reset to 0
       	cellX = 0;
	   	cellY = 0;
       	cellZ = 0;

    	if( cellX < 0 || cellX >= sizeX || cellY < 0 || cellY >= sizeY ||cellY < 0 ) {
	  	    printf("Error: outside scalar field boundaries\n");
		    v = 0.0;
		}
		else {
		    v = field[ ( cellY*sizeX ) + cellX  ];//////////////////////
		}
    	break;
   case 3:
     	//cellZ =  int( round( ( position[2] + (real)MP.boxsize[2] ) / (real)MP.ficellsize ) );
     	if(MP.ficellsize > 0){
     	cellX =  int( round( position[0] * sizeX / (real)MP.boxsize[0] ) );
       	cellY =  int( round( position[1] * sizeY / (real)MP.boxsize[1] ) );
       	cellZ =  int( round( position[2] * sizeX / (real)MP.boxsize[2] ) );
       	}
       	else {
       		cellX = 0;
       		cellY = 0;
       		cellZ = 0;
       	}
       	printf("Error: not implemented yet!\n");
       	/*
       	if( cellX < 0 || cellX >= sizeX || cellZ >= sizeZ ) {
   	  	    printf("Error: outside scalar field boundaries\n");
   		    v = 0.0;
   	}
   	//else {
   	//    v = field[ cellX ];
   	//   	}
    	*/

    	if( cellX < 0 || cellX >= sizeX || cellY < 0 || cellY >= sizeY ||cellY < 0 ) {
		    printf("Error: outside scalar field boundaries\n");
		    v = 0.0;
		}
		else {
		    v = field[ ( cellY * sizeY * sizeZ ) + ( cellY * sizeZ ) + cellZ ];//////////////////////
		}
    	break;
    }

  return v;
}

////////////////////////
//-------------------------------------------------------------------------
//elementary steps for all the microtubules in 'this' box
void ScalarField::stepField()
{

	printf("The boxshape: \t%d\t size: \t%d position_0: \t%d\n", MP.boxshape, MP.boxsize);

}



