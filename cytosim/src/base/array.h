//RCS: $Id: array.h,v 2.11 2005/04/22 12:31:09 nedelec Exp $

#ifndef ARRAY_H
#define ARRAY_H

#include "assert_macro.h"
#include "random.h"

///templated class Array<OBJ> : linear array of OBJ with resizing capability
/** 
    function array.allocate( size ) ensure that array[size-1] 
    can be accessed, Just like in normal C-array.
    array[ii] is the object stored at position ii.
    The mixWell() function randomly permutes the array

    If new memory is needed, a little more than requested (CHUNCK)
    is allocated. This makes this costly operation less frequent

    The class does not handle destruction of stored object
    Note: This is a templated class, due to limitation in g++,
    we need to keep all code inline (and in the header file).
*/
template<class OBJ>
class Array
{
  /// chuncks in which new memory is allocated
  static const int ACHUNCK = 4;
  
 private:

  /// array holding the OBJ
  OBJ * val;
  
  /// number of objects currently in the array
  int nbo;
  
  /// size of memory allocated for array val[]
  int allocated; 

 public:

  /// Default creator without allocation
  Array() 
  {
    allocated = 0;
    val = 0;
    nbo = 0;
  }
  
  /// Creator with allocation
  Array(int max) 
  {
    allocated = max; 
    val = new OBJ[max];
    if ( val == 0 ) {
      fprintf(stderr, "Array<T>::Array(int):: memory allocation failed\n");
      exit(1);
    }
    nbo = 0;
  }

  /// Destructor
  virtual ~Array() 
  { 
    erase(); 
  }

  /// returns the size() or number of objects
  int  size()     const   { return nbo; }
  
  /// the currently allocated size
  int  maxSize()  const   { return allocated; }

  ///the number of zero values in the array
  int  nbZeroValues() const 
  {
    int result = 0;
    for( int ii = 0; ii < nbo; ++ii )
      if ( val[ii] == 0 ) ++result;
    return result;
  }
  
  ///the number of non-zero objects in the list
  int  nbNonZeroValues() const   { return nbo - nbZeroValues(); }
    
  /// access operator, overloading the standard C-style
  OBJ & operator[](const int ii) const 
  {
    //verify that the given index is valid:
    assert( ii >=  0 );
    assert( ii < nbo );
    
    return val[ ii ];
  }
      
  /// allocation function, after, val[size-1] is valid
  int allocate(const int size) 
  {
    if ( val && ( size <= allocated ))
      return 0;
      
    //MSG("Array::allocate(%i) for list %p\n", size, this);
    // when new memory is needed, we more than needed
    // this saves some calls to this function:
    allocated = size + ACHUNCK;
    
    //create a new array:
    OBJ * val_new = new OBJ[ allocated ];
    if ( val_new == 0 ) {
      fprintf(stderr, "Array<T>::allocate(int):: memory allocation failed\n");
      exit(1);
    }    
    
    //copy the current values to the new array
    int ii = 0;
    if ( val ) {
      for(ii = 0; ii < nbo; ++ii) 
        val_new[ii] = val[ii];
      delete [] val;
    }
    
    //set the rest of the new array to zero
    //for(; ii < allocated; ++ii) 
    //  val_new[ii] = 0;
    
    //the new array becomes the current one
    val = val_new;
    
    //return the size that was allocated
    return allocated;
  }

  /// force the size of the array to N, filling with value def
  /** while allocate() does not change nbo, this does! */
  void setSize(const int size) 
  { 
    if ( size < nbo ) 
      nbo = size;
    else if ( size > nbo ) {
      allocate( size );
      nbo = size;
    }
  }
  
  /// force the size of the array to N, filling with value def
  /** while allocate() does not change nbo, this does! */
  void setSizeAndReset(const int size) 
  { 
    if ( size < nbo ) 
      nbo = size;
    else if ( size > nbo ) {
      allocate( size );
      do {
        val[ nbo++ ] = 0;
      } while( nbo < size );
    }
  }
    
  /// make a new entry at the end of the array, returning the index
  int pushFront(const OBJ np) 
  {
    if ( nbo >= allocated ) 
      allocate( nbo+1 );
    val[ nbo ] = np;
    return nbo++;
  }
  
  /// make a new entry at the end of the array, returning the new object
  OBJ & pushFront() 
  {
    if ( nbo >= allocated ) 
      allocate( nbo+1 );
    return val[nbo++];
  }
  
  /// make a new entry at the beginning of the array, returning the index
  int pushBack(const OBJ np)
  {
    for( int ii = 0; ii < nbo; ++ii ) {
      if ( val[ii] == 0 ) {
        val[ ii ] = nbo;
        return ii;
      }
    }
    if ( nbo >= allocated ) 
      allocate( nbo+1 );
    val[ nbo ] = np;
    return nbo++;
  }
  
  ///find an object in the list (slow scan search), return -1 if not found
  int find(const OBJ obj) const
  {
    for( int ii = 0; ii < nbo; ++ii )
      if ( val[ii] == obj ) return ii;
    return -1;
  }
  
  ///pack the array to ensure that no null value is stored
  void pack()
  {
    int stop = nbo;
    int nbo  = 0;
    int b    = 0;
    
    while( b < stop ) {
      //find next empty spot:
      while( val[nbo] != 0 )
        if ( ++nbo >= stop ) return;
      //from there, find the next filled spot:
      b = nbo+1;
      while( val[b] == 0 )
        if ( ++b >= stop ) return;
      
      val[nbo] = val[b];
      val[b] = 0;
      ++nbo;
      ++b;
    }          
  }
  
  
  /// clear all allocated values of the array
  void clearAll() 
  {
    for(int ii=0; ii < allocated; ++ii) 
      val[ii] = 0;
    nbo = 0; 
  }
  
  /// clear the values of the array until the current size
  void clear() 
  {
    for(int ii=0; ii < nbo; ++ii) 
      val[ii] = 0; 
    nbo = 0;
  }
  
  /// clear by just reseting the number of objects
  void clearFast() 
  { 
    nbo = 0;
  }
  
  /// clear, but also deallocate, i.e. deletes val[]
  void erase() 
  {
    if ( val ) {
      delete [] val;
      val = 0;
    }
    allocated = 0;
    nbo = 0;
  }


  /// put the last object on top, pushing everybody down one slot
  void turn() 
  {
    if ( nbo < 2 ) return;
    
    OBJ * tmp = val[ 0 ];
    for( int ii = 0; ii < nbo-1; ++ii ) 
      val[ ii ] = val[ ii+1 ];
    val[ nbo-1 ] = tmp;
  }

  /// permute two random values in the array
  void permute() 
  {
    int ii = RNG.pint_exc( nbo );
    int jj = RNG.pint_exc( nbo );
    if ( ii != jj ) {
      OBJ   tmp = val[ ii ];
      val[ ii ] = val[ jj ];
      val[ jj ] = tmp;
    }
  }
  
  
  /// mixWell() randomly permutes the objects in the array
  /** this produces uniform shuffling in linear time.
   see Knuth's The Art of Programming, Vol 2 chp. 3.4.2 
  */
  void mixWell() 
  {
    int jj = nbo, kk;
    while( jj > 1 ) {
      kk = RNG.pint_exc( jj );  //between 0 and j-1 
      --jj;
      OBJ  tmp = val[ jj ];
      val[ jj ] = val[ kk ];
      val[ kk ] = tmp;
    }
  }
};




#endif
