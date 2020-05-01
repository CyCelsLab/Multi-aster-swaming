//RCS: $Id: matblock.cc,v 2.4 2005/02/07 13:50:00 nedelec Exp $

#include "matblock.h"
#include "cblas.h"
#include <cstdlib>

//----------------------------------------------------------------------
void MatrixBlock::allocate() { 
  allocated  = 0;
  block_     = 0;
  block_alc  = 0;
  block_nb   = 0;
  block_size = 0;
  block_identity   = 0;
}

//----------------------------------------------------------------------
void MatrixBlock::deallocate() {

  if ( block_ == 0 ) return;

  for ( int ii=0; ii < block_nb; ++ii ) {
    if ( block_[ ii ] ) 
      delete[] block_[ ii ];
  }
  
  delete[] block_;
  delete[] block_size;
  delete[] block_alc;
  delete[] block_identity;
  allocated = 0;
}

//----------------------------------------------------------------------
void MatrixBlock::allocate( const int nbblock ) 
{
  assert( nbblock > 0 );
  block_nb = nbblock;
  if ( nbblock > allocated ) {
    
    int allocated_new = nbblock + BCHUNCK;
    //printf("new block-matrix sz %i\n", nbblock );
    real ** block_new           = new real*[ allocated_new ];
    int  *  block_alc_new       = new   int[ allocated_new ];
    int  *  block_size_new      = new   int[ allocated_new ];
    bool *  block_identity_new  = new  bool[ allocated_new ];
    
    if (( block_new == 0 ) || ( block_alc_new == 0 ) || ( block_size_new == 0 ) || (block_identity_new == 0)) {
      fprintf(stderr, "MatrixBlock::allocate(int):: memory allocation failed\n");
      exit(1);
    }      
    
    
    int ii = 0;
    
    if ( block_ ) {
      for( ii = 0; ii < allocated; ++ii) {
        block_new[ ii ]            = block_[ ii ];
        block_alc_new[ ii ]        = block_alc[ ii ];
        block_size_new[ ii ]       = block_size[ ii ];
        block_identity_new[ ii ]   = block_identity[ ii ];
      }
      delete[] block_;
      delete[] block_alc;
      delete[] block_size;
      delete[] block_identity;
    }

    for( ; ii < allocated_new; ++ii) {
      block_new[ ii ]            = 0;
      block_alc_new[ ii ]        = 0;
      block_size_new[ ii ]       = 0;
      block_identity_new[ ii ]   = false;
    }
    
    block_         = block_new;
    block_alc      = block_alc_new;
    block_size     = block_size_new;
    block_identity = block_identity_new;
    allocated      = allocated_new;
  }
}

//----------------------------------------------------------------------
void MatrixBlock::setBlockSize( const int ii, const int size, const bool treatBlockLikeIdentity ) 
{
  assert( block_ != 0 );
  assert( ( ii >= 0 ) && ( ii < block_nb ) );
  
  block_size[ ii ] = size;
  computeSize();

  if ( treatBlockLikeIdentity ) {
    
    block_identity[ ii ] = true;
    
  } else {
	
    block_identity[ ii ] = false;
    if ( block_alc[ ii ] < size ) {
      //printf("MatrixBlock::new block %i size %i\n", ii, size );
      block_alc[ ii ] = size + BCHUNCK;
      if ( block_[ ii ] )
        delete[] block_[ ii ];
      block_[ ii ] = new real[ block_alc[ii]*block_alc[ii] ];
      
      if ( block_[ ii ] == 0 ) {
        fprintf(stderr, "MatrixBlock::setBlockSize(int,int, bool):: memory allocation failed\n");
        exit(1);
      }      
      
    }
  
  }
}

//----------------------------------------------------------------------
int MatrixBlock::computeSize()
{
  size = 0;
  for(int ii = 0; ii < block_nb; ++ii ) 
	size += block_size[ii];
  return size;
}
	


//----------------------------------------------------------------------
real * MatrixBlock::addr(int ii, int jj) const {
  
  int bx = 0;  //block index on x
  int by = 0;  //block index on y
  
  while( ii >= block_size[bx] ) ii -= block_size[bx++];
  while( jj >= block_size[by] ) jj -= block_size[by++];
  
  if ( bx != by ) // the element is not on a diagonal block
    return 0;
  
  if ( block_identity[ bx ] ) //the block is identity
    return 0;
       
  if ( block_[ bx ] == 0 )  //the block in not allocated
    return 0;
  
  assert( ii < block_size[bx] );
  assert( jj < block_size[bx] );
  assert( block_size[bx] <= block_alc[bx] );
  
  return & block_[ bx ][ ii + block_size[bx] * jj ];
}

//----------------------------------------------------------------------
real & MatrixBlock::operator()( int x, int y)
{
  return * addr( x, y );
}

//----------------------------------------------------------------------
int MatrixBlock::nbNonZeroElements() const
{  
  int result=0;
  for( int ii = 0; ii < block_nb; ++ii )
	result += block_size[ ii ] * block_size[ ii ];
  return result;
}

//----------------------------------------------------------------------
int MatrixBlock::maxBlockSize()  const 
{
  int result = 0;
  for( int ii=0; ii < block_nb; ++ii )
	if ( result < block_size[ ii ] ) result = block_size[ ii ];
  return result;
}
  
//----------------------------------------------------------------------
void MatrixBlock::setBlockToZero( const int ii )
{
  real * BS; 
  if ( ( BS = block_[ ii ] ) )
	for(int kk = 0; kk < block_size[ ii ] * block_size[ ii ]; ++kk )
	  BS[ kk ] = 0;
}

//----------------------------------------------------------------------
void MatrixBlock::setToZero() {
  for( int ii = 0; ii < block_nb ; ++ii )
    setBlockToZero( ii );
}

//----------------------------------------------------------------------
void MatrixBlock::scaleBlock( const int ii, const real a )
{
  real * BS; 
  if ( ( BS = block_[ ii ] ) )
	for(int kk = 0; kk < block_size[ ii ] * block_size[ ii ]; ++kk )
	  BS[ kk ] = a * BS[ kk ];
}

//----------------------------------------------------------------------
void MatrixBlock::scale(const real a) {
  for( int ii = 0; ii < block_nb ; ++ii )
    scaleBlock( ii, a );
}

//----------------------------------------------------------------------
void MatrixBlock::multiplyVectorAdd( const real * X, real * Y )  const
{
  for(int xx = 0, ii = 0; ii < block_nb; xx += block_size[ ii++ ] )
	if ( block_identity[ ii ] )
	  blas_xcopy( block_size[ ii ], X+xx, 1, Y+xx, 1);
	else {
	  assert ( block_[ ii ] );
	  blas_xgemv('N', block_size[ ii ], block_size[ ii ], 1.0, 
                 block_[ ii ], block_size[ ii ], X+xx, 1, 1.0, Y+xx, 1);
	}
}

//----------------------------------------------------------------------
 void MatrixBlock::multiplyVector( const real * X, real * Y )  const
{
   for(int xx = 0, ii = 0; ii < block_nb; xx += block_size[ ii++ ] )
     if ( block_identity[ ii ] )
       blas_xcopy( block_size[ ii ], X+xx, 1, Y+xx, 1);
     else {
       assert( block_[ii] );
       blas_xgemv('N', block_size[ ii ], block_size[ ii ], 1.0, 
                  block_[ ii ], block_size[ ii ], X+xx, 1, 0.0, Y+xx, 1);
     }
 }

//----------------------------------------------------------------------
void MatrixBlock::multiplyVectorAddIsotropic( const real * X, real * Y )  const
{
  printf("MatrixBlock::multiplyVectorAddIsotropic() not implemented");
  exit(EXIT_FAILURE );
}

//----------------------------------------------------------------------
real MatrixBlock::maxNorm() const 
{
  real mx = 0;
  for(int b=0; b < block_nb; ++b)  {
    real * X = block_[ b ];
    int m = block_size[ b ];
    for(int ii = 0; ii < m*m; ++ii)
      if ( X[ ii ] > mx ) mx = X[ ii ];
  }
  return mx;
}


//----------------------------------------------------------------------
/// printf debug function in sparse mode: i, j : value
void MatrixBlock::printSparse(FILE * f) const
{
  int offset=0;
  for(int bx=0; bx < block_nb; ++bx) {
    const int bsize = block_size[bx];
    if ( block_identity[bx] ) {
      for(int ii=0; ii < bsize; ++ii)
        fprintf(f, "%i %i 1.0\n", ii+offset, ii+offset);
    } else {
      for(int ii=0; ii < bsize; ++ii)
        for(int jj=0; jj < bsize; ++jj)
          fprintf(f, "%i %i %16.8e\n", ii+offset, jj+offset, block_[bx][ ii + bsize * jj ]);
    }
    offset += bsize;
  }
}

