//RCS: $Id: map.h,v 2.10 2005/02/07 13:49:54 nedelec Exp $
// Map is a templated class to easily grid a finite rectangular region of space
// Francois Nedelec for Cytosim. This version Jan 2004. nedelec@embl.de

#ifndef MAP_H
#define MAP_H

#include <cstdio>
#include "types.h"
#include "assert_macro.h"
#include "iomessages.h"
#include "random.h"
#include "vecteur.h"
#include "iowrapper.h"

typedef unsigned char uchar;

///Map is a class to easily grid a finite region of a 1D, 2D or 3D space
/**
Map< int ORD, class CELL > associates one CELL to each grid division of
a rectangular piece in a space of dimensionality ORD.

The grid is regular, made of rectangular cells having all the same size
Access functions maps the real space onto a one-dimensional array,
holding type CELL which is specified as the second template argument
The cells are ordered sucessively, the x index varying the fastest
i.e. cell[ii+1] will in most cases be just right of cell[ii] along
axis X, but can also be one up on Y axis, if cell[ii] is already on the
right edge on the X axis.

Cells can be accessed in three ways:
 Position    a set of reals defining a point in space:  (x,y,z)
 Coordinates a set of ints defining a cell on the grid: (i,j,k)
 Index       is a int, defining a cell on the grid:      index

Conversions functions are given between these three types:
Assuming the grid size are (nbcells_x, nbcells_y, nbcells_z),
the size of the cells are calculated in delta:

  delta_x = ( right_x - left_x ) / nbcells_x;
  delta_y =  ... same with y
then
  i = int(  ( x - left_x ) / delta_x );
  j = ... same with y
  k = ... same with z
and
  index = i + nbcells_x * ( j + nbcells_y * k )

eg. a grid is numbered consecutively,
from left to right and bottom to top, then in the Z:
for a 4x4 2D grid, the index are like this:

12  13  14  15
8    9  10  11
4    5   6   7
0    1   2   3

Valid index are [0...nbCells()], where nbCells() is calculated by
setDimensions(). A point which is outside the rectangular region
defined by the bounds left and right, gives the closest cell in
indexFromPosition().
allocateEdgeTypeArray() calculate an int telling us what kind of border
each cell has, for the example above:

 7   6   6   8
 1   0   0   2
 1   0   0   2
 4   3   3   5

i.e. 0 is in the center, 1 is for left edge and 2 is for right edge
the total value is (value for X) + 3 * (value for Y)

The edge type of a cell is used to calculate the cells adjacent in space:
This is done by calling  allocateRegionArray();

After that, you can use for example:
void calculateRegionFromIndex(int list[], short & nbAdj, const int indx);
to fill in a list of cells adjacent to any cell in the grid:

for index = 0 it would return { 0 1 4 5 } in list[] and 4 in nbAdj
for index = 5 it would return { 5 0 1 2 4 6 8 9 10 } in list[] and nbAdj = 9

The map can be set to periodic boundary conditions with
setPeriodic();
for index = 0, it would now return { 0 15 12 13 3 2 7 4 5 } and  nbAdj = 9

cells can be accessed by:
    index        using operator [] with an argument int
    position     using the operator () with an argument  real[]
    coordinates  using the function cell( i, j, k );

EXAMPLE OF USAGE: a grid in 3D

Map< 3, someListClass > myGrid;

real left[]  = { -100, -100, -100 };   //lower range in each dimension
real right[] = { 100, 100, 100 };      //upper range in each dimension
int  size[]  = { 200,  200, 200 };     //number of cells in each dimension

myGrid.setDimensions(left, right, size);
myGrid.allocateRegionArray();

int nbAdjCells;    //number of adjacent cells (27 in 3D)
int * adjCells;    //list of adjacent (allocated by myGrid)


real aPosition[3] = { 1.3, 3.1, 5.0 }; //that would typically be a variable

int cellIndex = indexFromPosition( aPosition );
myGrid.setRegionFromIndex( adjCells, nbAdjCells, cellIndex );

for( int c = 0; c < nbAdjCells; ++c ) {
  someListClass  & local_list = myGrid[ adjCells[c] ];

  ... do something with this local_list;
}



TODO: if we would adjust the sizes to power of two, we might have a
faster conversions in indexToCoordinates.

The class is templated, and all the code need to be in the definition.
Sorry, this is a limitation of g++ in dealing with templates.
*/
template <int ORD, class CELL>
class Map
{
public:

    /// The possible mode for the wrapping of borders
    enum EdgeMode    { CONFINED = 0, PERIODIC = 1 };

    /// The type associated with cells at the edge
    enum BorderType  { CENTER = 0, LEFT_SIDE = 1, RIGHT_SIDE = 2 };



private:

    /// The number of cells in the map, i.e. allocation of val
    int maxIndx;

    /// The array of pointers to cells
    CELL * cells;

    /// The number of cells in each dimension of the real space
    int size[ ORD ];

    /// The real position of the left edge (min) in each dimension
    real left[ ORD ];

    /// The right position of the right edge (max) in each dimension
    real right[ ORD ];

    /// delta(d) = ( right(d)-left(d) ) / size( d )
    real delta[ ORD ];

    /// idelta(d) = 1.0 / delta(d)
    real idelta[ ORD ];

    //---------------------- EdgeType:
    /// tell if a cell is on a border / corner or in center

    /// The mode for the wrapping of borders
    EdgeMode edge_mode;

    /// Array holding the edge type of each cell, same size as cells
    uchar * edgeTypeForIndex;

    /// The maximum number of edge type = 3 to the power ORD
    uchar maxEdgeType;

    //---------------------- Region:
    //only adjacent cells are considered

    //because the numbering of cells is consitent, we can for
    //any cell find the adjacent cells by adding the same 'index offset'
    //to the initial cell index.
    //not all offsets are valid, however, and that depends on wheather
    //the cell is on the border of not. Below are list of Adjacent
    //offsets and their number as a function of the BorderType.

    /// max number of first neighbors
    short   maxNbAdj;

    /// number of neighboors as a function of edge type
    int  * nbAdjForEdgeType;

    /// array of index offset to neibhors, as a function of edge type
    int  * adjOffsetForEdgeType;

    /// a static array used to return in functions vicinity()
    int  * staticRegion;

    ///used to cleanup objects after reads
    short  mpFlag;

public:


    void test(){
		printf("@map.h, calling map \n");

	}

	//CHAITANYA: cells were orginall programmed to be private, now public!!
	/// The array of pointers to cells
    //CELL * cells;

    //--------------------------------------------------------------
    /// set all member arrays to zero
    void reset()
    {
      edge_mode            = CONFINED;
      // basic members:
      maxIndx              = 0;
      cells                = 0;
      // members needed for vicinity
      edgeTypeForIndex     = 0;
      nbAdjForEdgeType     = 0;
      adjOffsetForEdgeType = 0;
      staticRegion         = 0;
    }

    //--------------------------------------------------------------
    /// de-allocate all arrays
    void deleteAll()
    {
      deleteRegionArray();
      deleteEdgeTypeArray();
      deleteCellArray();
      edge_mode = CONFINED;
      maxIndx = 0;
    }

    //--------------------------------------------------------------
    /// Default creator
    Map()           { reset(); }

    //--------------------------------------------------------------
    /// Destructor
    virtual ~Map()  { deleteAll(); }

    //--------------------------------------------------------------
    /// set the sizes of the real space dimensions
    void setDimensions(const real min[ORD], const real max[ORD],
                       const int nbcells[ORD])
    {
      maxIndx    = 1;

      for (int d = 0; d < ORD; ++d) {

        if ( nbcells[d] <= 0 )
          MSG.error("map.h:setDimensions", "nbcells <= 0");

        if (  max[d] <= min[d] )
          MSG.error("map.h:setDimensions", "max <= min");

        maxIndx  *= nbcells[d];

        size[d]   = nbcells[d];
        left[d]   = min[d];
        right[d]  = max[d];
        delta[d]  = ( right[d] - left[d] ) / real( size[d] );
        idelta[d] = 1. / delta[d];
      }
    }


    //--------------------------------------------------------------
    /// create a map in 1D
    void create( real min, real max, int nbcells )
    {
      assert( ORD == 1 );
      setDimensions( &min, &max, &nbcells );
      allocateCellArray();
    }


    //--------------------------------------------------------------
    /// set the EdgeMode for periodic boundary conditions
    void setPeriodic()
    {
      if ( edge_mode != PERIODIC ) {
        edge_mode = PERIODIC;
        if ( staticRegion ) {
          deleteRegionArray();
          allocateRegionArray();
        }
      }
    }

    //--------------------------------------------------------------
    /// set the EdgeMode for confined boundary conditions
    void setConfined()
    {
      if ( edge_mode != CONFINED ) {
        edge_mode = CONFINED;
        if ( staticRegion ) {
          deleteRegionArray();
          allocateRegionArray();
        }
      }
    }

    //--------------------------------------------------------------
    /// the number of cells in the map
    int               nbCells()     const { return maxIndx; }

    //--------------------------------------------------------------
    /// the number of cells in the map
    CELL &            lastCell()    const { return cells[ maxIndx-1 ]; }

    //--------------------------------------------------------------
    /// the spatial position of the left border
    const real *      leftSide()    const { return left; }

    //--------------------------------------------------------------
    /// the position of the right border in space
    const real *      rightSide()   const { return right; }

    //--------------------------------------------------------------
    /// the width of the mapped region in space
    const real *      cellWidth()   const { return delta; }

    //--------------------------------------------------------------
    /// the length of the diagonal of a cell
    real diagonalCellLength() const
    {
      real result = 0;
      for( int d = 0; d < ORD; ++d )
        result += delta[ d ] * delta[ d ];
      return sqrt( result );
    }

    //--------------------------------------------------------------
    /// the volume of a cell
    real cellVolume() const
    {
      real result = 1;
      for( int d = 0; d < ORD; ++d)
        result *= delta[ d ];
      return result;
    }


    //----------getValue: CHAITANYA---------ONLY FLOATS---------------
    float getValue(int ii){
		float cval;
		cval = cells[ ii ];
		return cval;
	}

	//----------setValue: CHAITANYA---------ONLY FLOATS----------------
	void setValue(int ii, float cval ){
		//the value and the index position in the map

		cells[ ii ] = cval ;

		//return cells[ii];

	}


private:

    //--------------------------------------------------------------
    /// conversion from coordinates to index without border considerations
    int indexFromCoordinates_raw( const int coord[ORD] ) const
    {
      int result = 0;
      for( int d = ORD-1; d >= 0; --d ) {
        result = size[ d ] * result + coord[d];
      }
      return result;
    }

    //--------------------------------------------------------------
    /// conversion from coordinates to index without border considerations
    int indexFromCoordinates_confined( const int coord[ORD] ) const
    {
      register int c, result = 0;
      for( int d = ORD-1; d >= 0; --d ) {
        c = coord[ d ];
        if ( c < 0 )          c = 0;
        if ( c >= size[ d ] ) c = size[d]-1;
        result = size[ d ] * result + c;
      }
      assert( result < maxIndx );
      return result;
    }


    //--------------------------------------------------------------
    /// conversion from coordinates to index with periodic boundary cond.
    int indexFromCoordinates_periodic( int coord[ORD] ) const
    {
      register int c, result = 0;
      for( int d = ORD-1; d >= 0; --d ) {
        c = coord[ d ];
        while( c < 0 )        c += size[d];
        while( c >= size[d] ) c -= size[d];
        result = size[ d ] * result + c;
      }
      assert( result < maxIndx );
      return result;
    }


    //--------------------------------------------------------------
    /// returns the index of the cell whose center is closest to the point w[]
    int indexFromPosition_confined( const real w[ORD], const real shift=0) const
    {
      register int c, result = 0;
      for( int d = ORD-1; d >= 0; --d ) {
        c = int( shift + ( w[d] - left[d] ) * idelta[d] );
        if ( c < 0 )
          result *= size[d];
        else if ( c >= size[ d ] )
          result = size[d] * ( result + 1 ) -1;
        else
          result = size[d] * result + c;
      }
      assert( result < maxIndx );
      return result;
    }

    //--------------------------------------------------------------
    /// returns the index of the cell whose center is closest to the point w[]
    int indexFromPosition_periodic( const real w[ORD], const real shift=0) const
    {
      register int c, result = 0;
      for( int d = ORD-1; d >= 0; --d ) {
        c = int( shift + ( w[d] - left[d] ) * idelta[d] );
        while( c < 0 )        c += size[d];
        while( c >= size[d] ) c -= size[d];
        result = size[ d ] * result + c;
      }
      assert( result < maxIndx );
      return result;
    }

public:
    //--------------------------------------------------------------
    /// direct access to a cell by index
    CELL & operator[](const int indx) const
    {
      assert(( indx >= 0 ) && ( indx < maxIndx ));
      return cells[ indx ];
    }


    /*
    /// direct setting of a cell by index
	CELL & operator[](int i, float val)
	{
	   assert(( i >= 0 ) && ( i < maxIndx ));
	   return cells[ i ] = val;
    }
    */

    //--------------------------------------------------------------
    /// conversion from coordinates to index with border
    int indexFromCoordinates( const int coord[ORD] ) const
    {
      if ( edge_mode == CONFINED )
        return indexFromCoordinates_confined( coord );
      else
        return indexFromCoordinates_periodic( coord );
    }

    //--------------------------------------------------------------
    /// conversion from coordinates to index with border
    int indexFromPosition( const real w[ORD] ) const
    {
      if ( edge_mode == CONFINED )
        return indexFromPosition_confined( w );
      else
        return indexFromPosition_periodic( w );
    }

    //--------------------------------------------------------------
    /// conversion from coordinates to index with border
    int indexFromPosition_shifted( const real w[ORD] ) const
    {
      if ( edge_mode == CONFINED )
        return indexFromPosition_confined( w, 0.5 );
      else
        return indexFromPosition_periodic( w, 0.5 );
    }

    //--------------------------------------------------------------
    /// fast cell index by coordinates in 1D
    int indexFromCoordinates( int x ) const
    {
      assert ( ORD == 1 );
      if ( edge_mode == CONFINED ) {
        if ( x < 0 )          x = 0;
        if ( x >= size[0] )   x = size[0]-1;
      } else {
        while( x < 0 )        x += size[0];
        while( x >= size[0] ) x -= size[0];
      }
      return x;
    }

    ///direct access to cell in 1D
    CELL & cell(int x) {
      return cells[ indexFromCoordinates( x ) ];
    }

    //--------------------------------------------------------------
    /// fast cell index by coordinates in 2D
    int indexFromCoordinates( int x, int y ) const
    {
      assert ( ORD == 2 );
      if ( edge_mode == CONFINED ) {
        if ( x < 0 )          x = 0;
        if ( x >= size[0] )   x = size[0]-1;
        if ( y < 0 )          y = 0;
        if ( y >= size[1] )   y = size[1]-1;
      } else {
        while( x < 0 )        x += size[0];
        while( x >= size[0] ) x -= size[0];
        while( y < 0 )        y += size[1];
        while( y >= size[1] ) y -= size[1];
      }
      return x + size[0] * y;
    }

    ///direct access to cell in 2D
    CELL & cell(int x, int y) {
      return cells[ indexFromCoordinates( x, y ) ];
    }
    //--------------------------------------------------------------
    /// fast cell index by coordinates in 3D
    int indexFromCoordinates( int x, int y, int z ) const
    {
      assert ( ORD == 3 );
      if ( edge_mode == CONFINED ) {
        if ( x < 0 )          x = 0;
        if ( x >= size[0] )   x = size[0]-1;
        if ( y < 0 )          y = 0;
        if ( y >= size[1] )   y = size[1]-1;
        if ( z < 0 )          z = 0;
        if ( z >= size[2] )   z = size[2]-1;
      } else {
        while( x < 0 )        x += size[0];
        while( x >= size[0] ) x -= size[0];
        while( y < 0 )        y += size[1];
        while( y >= size[1] ) y -= size[1];
        while( z < 0 )        z += size[2];
        while( z >= size[2] ) z -= size[2];
      }
      return x + size[0] * ( y + size[1] * z );
    }

    ///direct access to cell in 3D:
    CELL & cell(int x, int y, int z) const
    {
      return cells[ indexFromCoordinates( x, y, z ) ];
    }

    //--------------------------------------------------------------
    /// access to a cell by position
    CELL & operator()(const real w[ORD]) const
    {
        return cells[ indexFromPosition(w) ];
    }

    //--------------------------------------------------------------
    /// access to a cell by position in 1D
    CELL & operator()(const real x) const
    {
        assert( ORD == 1 );
        return cells[ indexFromPosition(&x) ];
    }

    //--------------------------------------------------------------
    /// access to a cell by position in 2D
    CELL & operator()(const real x, const real y) const
    {
      assert( ORD == 2 );
      real w[ORD] = { x, y };
      return cells[ indexFromPosition(w) ];
    }

    //--------------------------------------------------------------
    /// access to a cell by position in 3D
    CELL & operator()(const real x, const real y, const real z) const
    {
      assert( ORD == 3 );
      real w[ORD] = { x, y, z };
      return cells[ indexFromPosition(w) ];
    }

    //--------------------------------------------------------------
    /// conversion from index back to coordinates
    void setCoordinatesFromIndex(int coord[ORD], int indx) const
    {
      for( int d = 0; d < ORD; ++d ) {
        coord[d] = indx % size[d];
        indx   /= size[d];
      }
    }


    //--------------------------------------------------------------
    /// conversion from Position to coordinates, return true if inside
    int setCoordinatesFromPosition(int coord[ORD], const real w[ORD], const real shift=0) const
    {
      if ( edge_mode == CONFINED ) {
        int inside = 1;
        for ( int d = 0; d < ORD; ++d ) {
          int c = int( shift + ( w[d] - left[d] ) * idelta[d] );
          if (( c <= 0 ) || ( c >= size[d] )) inside = 0;
          coord[ d ] = c;
        }
        return inside;
      } else {
        for ( int d = 0; d < ORD; ++d ) {
          int c = int( shift + ( w[d] - left[d] ) * idelta[d] );
          while( c < 0 )        c += size[0];
          while( c >= size[0] ) c -= size[0];
          coord[ d ] = c;
        }
        return 1;
      }
    }

    //--------------------------------------------------------------
    /// conversion from Index to Position
    void setPositionFromIndex( real w[ORD], int index, int center ) const
    {
      for( int d = 0; d < ORD; ++d ) {
        int coord = index % size[d];
        if ( center )
          w[d] = left[d] + delta[d] * ( 0.5 + coord );
        else
          w[d] = left[d] + delta[d] * ( coord );
        index /= size[d];
      }
    }

    //--------------------------------------------------------------
	///CHAITANYA Returns the VECTOR position (x,y,z) for a given index
	Vecteur positionFromIndex( int index, int center ) const
	{
		Vecteur vv;
		for( int d = 0; d < ORD; ++d ) {
			int coord = index % size[d];
			if ( center )
			    vv[d] = left[d] + delta[d] * ( 0.5 + coord );
		    else
		        vv[d] = left[d] + delta[d] * ( coord );
		    index /= size[d];
      }

	return vv;
    }


    //--------------------------------------------------------------
    /// conversion from Coordinates to Position
    void setPositionFromCoordinates( real w[ORD], const int coord[ORD], int center ) const
    {
      for ( int d = 0; d < ORD; ++d ) {
        if ( center )
          w[d] = left[d] + delta[d] * ( 0.5 + coord[d] );
        else
          w[d] = left[d] + delta[d] * coord[d];
      }
    }


    //--------------------------------------------------------------
    /// allocate the array of cells
    void allocateCellArray()
    {
      if ( maxIndx == 0 )
        printf("maxIndx==0 in allocateCellArray() : call setDimensions() first\n");

      if ( cells )
        delete[] cells;

      cells = new CELL[ maxIndx ];

      if ( cells == 0 ) {
        fprintf(stderr, "Map<,>::allocateCellArray():: memory allocation failed\n");
        exit(1);
      }
    }

	//--------------------------------------------------------------
	/// get maximum value in a map- for a map of reals!
	/// CHAITANYA
	real getMaxval(){
		real maxVal=cells[ 0 ];
		for (int mm = 0; mm < maxIndx; mm++){
			if( cells[ mm ] > maxVal)
			   maxVal = cells[ mm ];
			}
		return maxVal;
	}
	//--------------------------------------------------------------
	/// get maximum value in a map- for a map of reals!
	/// CHAITANYA
	real getMinval(){
		real minVal = cells[ 0 ];//initialize with first value
		for (int mm = 0; mm < maxIndx; mm++){
			if( cells[ mm ] <  minVal || minVal == 0)
			   minVal = cells[ mm ];
			}
		return minVal;
	}



    //--------------------------------------------------------------
    /// returns true if the array cells[] has been allocated
    bool cellsAllocated() const
    {
      return ( cells != 0 );
    }


    //--------------------------------------------------------------
    /// erase the array 'cells'
    void deleteCellArray()
    {
      if ( cells )
        delete[] cells;
      cells = 0;
    }


    //===================================================================
    //========================EDGE TYPE==================================

    /// the edge type of a particular cell by its index
    uchar edgeType( const int indx ) const
    {
      assert( edgeTypeForIndex );
      assert((0 <= indx ) && ( indx < maxIndx ));
      return edgeTypeForIndex[ indx ];
    }

    //--------------------------------------------------------------
    /// caculate the edge type from the coordinates
    uchar edgeTypeFromCoordinates( const int coord[ORD] ) const
    {
      uchar edgeType = CENTER;
      for (int d = ORD-1; d >= 0; --d) {
        edgeType *= 3;
        if ( coord[d] <= 0 )
          edgeType += LEFT_SIDE;
        else if ( coord[d] >= size[d]-1 )
          edgeType += RIGHT_SIDE;
      }
      return edgeType;
    }

    //--------------------------------------------------------------
    /// set-up the edge array
    void allocateEdgeTypeArray()
    {
      if ( maxIndx == 0 )
        printf("maxIndx==0 in allocateEdgeTypeArray() : call setDimensions() first\n");

      maxEdgeType = 1;
      for( int d = 0; d < ORD; ++d )
        maxEdgeType *= 3;
      edgeTypeForIndex = new uchar[ maxIndx ];

      if ( edgeTypeForIndex == 0 ) {
        fprintf(stderr, "Map<,>::allocateEdgeTypeArray():: memory allocation failed\n");
        exit(1);
      }

      int coord[ ORD ];
      for (int index = 0; index < maxIndx; index++) {
        setCoordinatesFromIndex( coord, index );
        edgeTypeForIndex[ index ] = edgeTypeFromCoordinates( coord );
        assert( edgeTypeForIndex[ index ] <= maxEdgeType );
      }
    }

    //--------------------------------------------------------------
    /// erase the edge array
    void deleteEdgeTypeArray()
    {
      if ( edgeTypeForIndex )
        delete[] edgeTypeForIndex;
      edgeTypeForIndex = 0;
    }


    //===================================================================
    //=================== OFFSET TO ADJACENT CELLS ======================

    /// calculate the index offset to neareset neighbors
    void setAdjacentOffsetFromEdgeType(int * array, int & nb, uchar edgeType) const
    {
      assert( maxEdgeType > 0 );

      //extract the edge type for each dimension:
      short edge[ ORD ];
      short b = edgeType;
      for (int d = 0; d < ORD; ++d) {
        edge[ d ] = b % 3;
        b /= 3;
      }


      nb = 0;
      int dx[ORD];
      for ( int power = 0; power < maxEdgeType; ++power ) {

        int reject = 0;
        int jj = power;

        for (int d = 0; d < ORD; ++d) {
          dx[d] = jj % 3 - 1;
          jj /= 3;

          if ( ( size[d] < 2 ) && ( dx != 0 ) )           reject=1;
          if ( ( size[d] < 3 ) && ( edge[d] == CENTER ) ) reject=1;

          if ( edge_mode == PERIODIC ) {
            if ( dx[d] == -1 ) {
              if ( edge[d] == LEFT_SIDE ) dx[d] = size[d]-1;
            } else if ( dx[d] ==  1 ) {
              if ( size[d] < 3 )
                reject=1;
              else
                if ( edge[d] == RIGHT_SIDE ) dx[d] = 1-size[d];
            }
          } else {
            if ( ( dx[d] == -1 ) && ( edge[d] == LEFT_SIDE  ) ) reject=1;
            if ( ( dx[d] ==  1 ) && ( edge[d] == RIGHT_SIDE ) ) reject=1;
          }
        }

        //calculate the offset, store it in the array
        if (reject == 0)
          array[ nb++ ] = indexFromCoordinates_raw( dx );
        //printf("%i %i : %i %i : %2i %2i : %i %i\n", edgeType, nb-1, edge[0], edge[1], dx[0], dx[1], array[ nb-1], reject );
      }

      assert ( nb <= maxNbAdj );
      //put the value zero (center cell) on top of the array:
      for (int jj = 1; jj < nb; ++jj)
        if ( array[jj] == 0 ) {
          array[jj] = array[0];
          array[0] = 0;
          break;
        }

    }


    //--------------------------------------------------------------
    /// calculate the index offset for all Edge type
    void allocateRegionArray()
    {
      if ( edgeTypeForIndex == 0 )
        allocateEdgeTypeArray();

      maxNbAdj = 1;
      for (int d = 0; d < ORD; ++d) {
        if (size[d]==2) maxNbAdj *= 2;
        if (size[d]>2)  maxNbAdj *= 3;
      }

      staticRegion         = new int[ maxNbAdj ];
      nbAdjForEdgeType     = new int[ maxEdgeType+1 ];
      adjOffsetForEdgeType = new int[ maxNbAdj * (maxEdgeType+1) ];

      if (( staticRegion == 0 ) || ( nbAdjForEdgeType == 0 ) || ( adjOffsetForEdgeType == 0 )) {
        fprintf(stderr, "Map<,>::allocateRegionArray():: memory allocation failed\n");
        exit(1);
      }

      for (uchar et = 0; et <= maxEdgeType; ++et)
        setAdjacentOffsetFromEdgeType( &adjOffsetForEdgeType[ maxNbAdj * et ], nbAdjForEdgeType[ et ], et );
    }


    //--------------------------------------------------------------
    /// erase the array of index offset
    void deleteRegionArray()
    {
      if ( staticRegion ) {
        delete[] staticRegion;
        delete[] nbAdjForEdgeType;
        delete[] adjOffsetForEdgeType;
      }
      staticRegion         = 0;
      nbAdjForEdgeType     = 0;
      adjOffsetForEdgeType = 0;
    }


    //===================================================================
    //==================== LIST OF ADJACENT CELLS =======================


    /// returns the list of offsets to the nearest neighbors
    void getRegionOffsetFromIndex(int *& list, short & nbAdj, const int indx, const int mix=1)
    {
      assert( 0 <= indx ); assert( indx < maxIndx );
      assert( edgeTypeForIndex );
      assert( adjOffsetForEdgeType );
      assert( list );

      int edge = edgeTypeForIndex[indx];
      list  = adjOffsetForEdgeType + edge * maxNbAdj;
      nbAdj = nbAdjForEdgeType[ edge ];
      if ( mix ) RNG.mixWell( list, nbAdj );
    }

    /// calculate a list of nearest neighbors
    /** the user must make sure that list[] is correctly allocated */
    void calculateRegionFromIndex(int list[], short & nbAdj, const int indx, const int mix=1)
    {
      assert( 0 <= indx ); assert( indx < maxIndx );
      assert( edgeTypeForIndex );
      assert( adjOffsetForEdgeType );
      assert( list );

      int edge = edgeTypeForIndex[indx];
      int * adjOffset = adjOffsetForEdgeType + edge * maxNbAdj;
      nbAdj = nbAdjForEdgeType[ edge ];

      if ( mix ) RNG.mixWell( adjOffset, nbAdj );
      for ( int ii = 0; ii < nbAdj; ++ii )
        list[ ii ] = indx + adjOffset[ ii ];
    }

    //--------------------------------------------------------------
    /// get a list of nearest neighbors
    /** use this if the order in which cells are considered is meaningless,
        for example to calculate a sum of all forces with all objects */
    void setRegionFromIndex(int *& list, short & nbAdj, const int indx, const int mix=1)
    {
      assert( 0 <= indx ); assert( indx < maxIndx );

      calculateRegionFromIndex( staticRegion, nbAdj, indx, mix );
      list = staticRegion;
    }

    //===================================================================
    //======================== EMPTY  REGION ============================

    //--------------------------------------------------------------
    /// set all values of cells[] to zero
    virtual void clear()
    {
      for (int ii = 0; ii < maxIndx; ++ii)
        cells[ ii ] = 0;
    }


    //--------------------------------------------------------------
    /// the sum of all cell, works only if CELL is arithmetic
    CELL sum() const
    {
      CELL result = 0;
      for (int ii = 0; ii < maxIndx; ++ii)
        result += cells[ii];
      return result;
    }


    //--------------------------------------------------------------
    /// write to a file, works only if CELL is float or double
    void write( FILE * file, int center ) const
    {
      real w[ ORD ];
      for (int ii = 0; ii < maxIndx; ++ii ) {
        setPositionFromIndex(w, ii, center);
        fprintf(file, "%i ", ii);
        for (int d=0; d < ORD; ++d)
          fprintf(file, " %7.2f", w[ d ]);
          fprintf(file,"  %f\n", cells[ ii ]);
      }
    }


    /// write to a file, works only if CELL is float or double
   void read( FILE * file, int center ) const
    {
       //need to include methods which could alternative to result.out be used to read in map


    }


	    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
	    //------------------------------------------------------------------------
		//-----------------------------------FILE---------------------------------
		//------------------------------------------------------------------------

//---------------------I/O START

//Chaitanya Athale 23-oct-2008
/// write using IOwrapper methods- for writing into result file

  void write() const
	   {
		   //
		   //write 'this' on file "s" as a list of rods
		   //IO.writeRecordTag( "mp" );

		   //printf("@map.h write() maxInd: %i\n",maxIndx );

		   IO.writeUInt32( maxIndx );


		   int center = 0;

		   for (int ii = 0; ii < maxIndx; ++ii ) {
				//setPositionFromIndex(w, ii, center);
				//for (int d=0; d < ORD; ++d){
				//	fprintf(file, " %7.2f", w[ d ]);
				//}
				IO.writeUInt32( ii );
				IO.writeReal32( cells[ ii ] );
				//printf("@map.h::write() cellV[%i]= %5.5f \n",ii,cells[ ii ]  );
	        }



	   	 }


	//read the map index and cell-values
	void read()
		{
           //not implemented

	   	 }







//---------------------I/O END

//--------------------------------------------------------------
//printf() function for checking/debug

	void dump( int dump_offset = 0 )
    {
      int jj;
      int * v = 0;
      short nv;

      real w[ORD];
      int  coord[ORD];
      printf("dump of Map %p with ORD=%i isPeriodic = %i\n", (void*)this, ORD, edge_mode);

      if ( dump_offset )
        for ( uchar et = 0; et < maxEdgeType; ++et ) {
          printf("edge type %i : offsets ", et);
          for( jj=0; jj < nbAdjForEdgeType[ et ]; ++jj )
            printf("%3i", adjOffsetForEdgeType[ maxNbAdj * et + jj] );
          printf("\n");
        }

      for (int index = 0; index < maxIndx; index++) {

        printf("cell %2d", index);

        if ( edgeTypeForIndex )
          printf(" : bor %2d", edgeTypeForIndex[index]);

        printf(" : ijk");
        setCoordinatesFromIndex( coord, index );
        for( jj=0; jj<ORD; ++jj) printf(" %2d", coord[jj]);

        printf(" : xyz");
        setPositionFromIndex( w, index, 1 );
        for( jj=0; jj<ORD; ++jj) printf(" %2.2f", w[jj]);

        if ( staticRegion ) {
          printf(" : neib");
          setRegionFromIndex(v, nv, index, 0);
          for( jj=0; jj<nv; ++jj) printf(" %2d", v[jj]);
        }

        printf("\n");
      }
    }

    //---------------------------------------------------------
	//                  Flag functions
    //---------------------------------------------------------

	///Set the Flag
	void   setFlag(short value)     { mpFlag = value; }

	///Return the Flag
    int    getFlag()  const { return mpFlag; }


};


#endif
