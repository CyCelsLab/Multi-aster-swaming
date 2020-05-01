//RCS: $Id: name_list.cc,v 2.3 2005/02/07 13:50:27 nedelec Exp $

#include "name_list.h"
#include <cstdlib>

//------------------------------------------------------------
//beware: allocate(size) ensures that we can accesss byNames[ size ]
//this is not the usual C convention, which allocate to byNames[ size-1 ]
void NameList::allocate( Name size )
{
  if ( size > allocated ) {
    
    Name    allocated_new = size + NACHUNCK;
    Node **   byNames_new = new Node*[ 1 + allocated_new ];
    
    if ( byNames_new == 0 ) {
      fprintf(stderr, "NameList::allocate(int):: memory allocation failed\n");
      exit(1);
    }
    
    Name x = 0;
    for( ; x <= allocated;     ++x ) byNames_new[ x ] = byNames[ x ];
    delete[] byNames;
        
    for( ; x <= allocated_new; ++x ) byNames_new[ x ] = 0;
    byNames   = byNames_new;
    allocated = allocated_new;
    
  }
}

//------------------------------------------------------------
NameList::NameList() 
{
  next      = 1;
  allocated = 1;
  byNames   = new Node*[ 1 + allocated ];
    
  if ( byNames == 0 ) {
    fprintf(stderr, "NameList::NameList():: memory allocation failed\n");
    exit(1);
  }
  
  for( Name ii = 0; ii <= allocated; ++ii )
    byNames[ ii ] = 0;
}

//------------------------------------------------------------
NameList::~NameList()
{
  delete[] byNames;
}

//------------------------------------------------------------
//a name is normally given only once 
//(unless renameAll() or forgetAll() is called):
//names of deleted nodes are not recycled, to simplify 
//the analysis of simulations where dynamic objects are
//deleted/created continuously
Name NameList::nameIt( Node * obj )
{
  allocate( next );
  assert( byNames[ next ] == 0 );
  byNames[ next ] = obj;
  //MSG("named %lx %li\n", obj, next);
  obj->name = next++;
  return obj->name;
}

//------------------------------------------------------------
Name NameList::firstUnusedName()
{
  Name result = 1;
  while( result <= allocated ) {
    if( byNames[ result ] == 0 )
      return result;
    else
      ++result;
  }
  allocate( result );
  return result;
}

//------------------------------------------------------------
Name NameList::firstUsedName() const
{
  Name result = 1;
  while( result <= allocated ) {
    if( byNames[ result ] )
      return result;
    else
      ++result;
  }
  return 0;
}

//------------------------------------------------------------
Name NameList::lastUsedName() const
{
  Name result = allocated;
  while( result > 0 ) {
    if( byNames[ result ] )
      return result;
    else
      --result;
  }
  return 0;
}


//------------------------------------------------------------
void NameList::remember( Node * obj )
{
  if ( 0 == obj->name ) {
    nameIt( obj );
  } else {
    allocate( obj->name );
    if ( next <= obj->name )
      next = obj->name + 1;
    if (( byNames[ obj->name ] ) && ( byNames[ obj->name ] != obj ))
	    MSG.error("NameList::remember()", "conflict: name %i is already registered", obj->name );
    byNames[ obj->name ] = obj;
  }
}

//------------------------------------------------------------
void NameList::forget( Node * obj )
{
    assert( obj->name );
    assert( allocated >= obj->name );
    assert( byNames[ obj->name ] == obj );
    
    byNames[ obj->name ] = 0;
}

//------------------------------------------------------------
Node * NameList::nodeWithThisName( Name name ) const
{
  if (( name >= 0 ) && ( name <= allocated )) {
    assert(( byNames[name] == 0 ) || ( byNames[ name ] -> name == name ));
    return byNames[ name ];
  } else
    return 0;
}

//------------------------------------------------------------
void NameList::forgetAll()
{
  for( Name x = 0; x <= allocated; ++x )
    byNames[ x ] = 0;
  next = 1;
}

//------------------------------------------------------------
void NameList::renameAll()
//rename all objects in a continuous maner, preserving the current order.
//i.e. names are repacked from 1 to (nb of objects)
{
  next = 1;
  Name nn = 1;
    
  while( nn <= allocated ) {
    while (( nn <= allocated ) && ( byNames[ nn ] == 0 )) ++nn;
    if ( nn > allocated ) return;
    if ( next < nn ) {
      byNames[ next ] = byNames[ nn ];
      byNames[ next ]->name = next;
      byNames[ nn ] = 0;
    }
    ++next;
    ++nn;
  }
}
