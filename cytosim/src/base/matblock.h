//RCS: $Id: matblock.h,v 2.4 2005/01/18 11:19:48 nedelec Exp $
//===========================================================================

#ifndef MATBLOCK_H
#define MATBLOCK_H

#include "matrix.h"
#include <cstdio>
#include "assert_macro.h"

///a non-symmetric real Matrix, diagonal by blocks.
class MatrixBlock : public Matrix
{
  
  static const int BCHUNCK = 2;

private:

  /// number of blocks
  int block_nb;
  
  /// array of pointers to the blocks
  real * * block_;

  /// array containing the size of each block
  int * block_size; 
  
  /// array specifying the allocated size of each block
  int * block_alc;
  
  /// array specifying if block should be considered equal to one
  bool * block_identity;

public:

  /// allocation
  void allocate();

  /// the deallocation
  void deallocate();

  /// allocate the matrix to be able to hold nb_block (arg 1) blocks
  void allocate( const int nb_block );
  
  /// allocate block ii (arg 1) to be capable of holding (size*size) (arg 2)
  void setBlockSize( const int ii, const int size, const bool treatBlockLikeIdentity=false );
  
  /// default constructor
  MatrixBlock() { allocate(); }

  /// default destructor
  virtual ~MatrixBlock() { deallocate(); }

  
  /// return the address of first element in block ii
  real * getBlock( const int ii ) {
    assert( block_ );
    assert(( ii >= 0 ) && ( ii < block_nb ));
    assert( block_[ ii ] );
    assert( block_identity[ii] == false );
    assert( block_size[ii] <= block_alc[ii] );
    return block_[ ii ];
  }

  /// returns the number of blocks
  int getNbBlocks() const {
    return block_nb;
  }
  
  /// return the total size of the matrix
  int computeSize();
  
  /// true if block ii will be considered equal to Identity by multiply()
  bool treatBlockLikeIdentity( const int ii ) const {
    return block_identity[ ii ];
  }
  
  /// returns the size of block ii
  int getBlockSize( const int ii ) const {
    assert( block_ );
    assert(( ii >= 0 ) && ( ii < block_nb ));
    return block_size[ ii ];
  }	

  /// the address holding element (ii, jj)
  real * addr( int ii, int jj ) const;

  /// returns the address of element at (x, y), allocating if necessary
  virtual real & operator()( int x, int y);

  /// reset all the values in block ii
  void setBlockToZero( const int ii );
  
  /// reset the entire matrix
  void setToZero();
 
  /// scale block ii
  void scaleBlock( const int ii, const real a );
  
  /// scale the entire matrix
  void scale( const real a );
  
  /// total number of elements allocated
  int  nbNonZeroElements() const;

  /// size of the biggest block
  int  maxBlockSize() const;
  
  /// vector multiplication: Y <- M * X
  void multiplyVectorAdd( const real * X, real * Y) const;

  /// vector multiplication: Y <- M * X
  void multiplyVector( const real * X, real * Y) const;

  /// vector multiplication: Y <- M * X
  void multiplyVectorAddIsotropic( const real * X, real * Y) const;

  /// maximum of the absolute value of all elements
  real maxNorm() const;
  
  /// printf debug function in sparse mode: i, j : value
  void printSparse(FILE * f=stdout) const;

};



#endif
