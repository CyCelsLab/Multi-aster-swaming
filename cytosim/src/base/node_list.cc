//RCS: $Id: node_list.cc,v 2.3 2005/02/02 19:27:42 nedelec Exp $
//---------------------------------node.cc----------------------------------
// doubly linked list, STL style, with acces by iterators,
// some additions to manipulate the list: sorting, unsorting, etc.

#include "node_list.h"
#include "random.h"
#include "iomessages.h"
#include <cstring>

//-------------------------------------------------------------------------
Node * NodeList::pushFront(Node * np)
{
  assert( np );
  assert( np->mylist == 0 );
  //MSG("NodeList: push_front   %lx in   %lx\n", np, this);
  
  np->son = first;
  np->dad = 0;
  if ( first ) {
    first->dad = np;
    //TODO: test here if the new node type is compatible with the list
  } else { 
    last = np;
  }
  first = np;

  ++nbo;
  np->mylist = this;
  if ( namelist )
    namelist -> remember( np );
  return np;
}

//-------------------------------------------------------------------------
Node * NodeList::pushBack(Node * np)
{    
  //MSG("NodeList: push_back   %lx in   %lx\n", np, this);
  assert( np->mylist == 0 );
    
  np->son = 0;
  np->dad = last;
  if (last) {
    last->son = np;
    //TODO: test here if the new node type is compatible with the list
  } else {  
    first = np;
  }
  last = np;

  ++nbo;
  np->mylist = this;
  if ( namelist )
    namelist -> remember( np );
  return np;
}

//-------------------------------------------------------------------------
Node * NodeList::pop(Node * op)
{
  //MSG("NodeList: pop node %lx name %lx\n", op, op->name);
  assert( op->mylist == this );
  
  if (op->dad)
    (op->dad)->son = op->son;
  else
    first = op->son;
  
  if (op->son)
    (op->son)->dad = op->dad;
  else
    last = op->dad;
      
  op->dad    = 0;
  op->son    = 0;
  op->mylist = 0;

  --nbo;
  if ( namelist ) 
    namelist -> forget( op );
  //we do not set name = 0, since that might be used later
  return op;
}

//-------------------------------------------------------------------------
void NodeList::popAndPush(Node * np)
  //this assumes that np is linked in a list sharing the same
  //namelist as *this, hence we skip the namelist update parts
{
  assert( np->mylist );
  assert( namelist == np->mylist->namelist );

  //--------------pop() from the current list:
  if (np->dad)
    (np->dad)->son    = np->son;
  else
    np->mylist->first = np->son;
  
  if (np->son)
    (np->son)->dad   = np->dad;
  else
    np->mylist->last = np->dad;
      
  --(np->mylist->nbo);

  //-------------push_getFront() from this list:
  np->son = first;
  np->dad = 0;
  if (first)
    first->dad = np;
  else  
    last = np;
  first = np;

  ++nbo;
  np->mylist = this;
}

//-------------------------------------------------------------------------
void NodeList::clear() 
{
  while( first )
    pop( first );
  assert( nbo == 0 );
}

//-------------------------------------------------------------------------
void NodeList::erase() 
{
  //it is possible that deletion of objects lead to creation of new ones
  //because new objects can be allocated in a destructor
  //we try to detect this 'resurections' to avoid an infinite loop here
  for(int nboi = nbo; nboi > 0; --nboi) 
    delete( pop(first) );
  if ( nbo > 0 )
    MSG.error("NodeList::erase()", "Resurection detected while deleting objects");
  //do not call namelist.forgetAll(), since namelist might have objects for other nodelist;
}

//-------------------------------------------------------------------------
void NodeList::deleteIfFlagDifferentThan( int flagToKeep )
{
  Node * i = first, * j;
  
  while( i ) {
    j = i -> son;
    if ( flagToKeep != i->ndFlag )
      delete( pop( i ));
    i = j;
  }
}

//-------------------------------------------------------------------------
Node * NodeList::nodeWithThisName( Name n ) const
{
  if ( namelist ) {
    Node * result = namelist->nodeWithThisName( n );
    if ( result && ( result->mylist == this ))
      return result;
    else
      return 0;
  }
  return 0;
}

//-------------------------------------------------------------------------
Node * NodeList::nodeWithThisName_scan( Name n ) const
  //search without using the namelist array, just slow scanning:
{
  Node * i = first;
  while( i ) {
    if ( i->name == n ) return i;
    i = i->son;
  }
  return 0;
}

//-------------------------------------------------------------------------
void NodeList::turn()
  
  //put the last on top of the list
{
  if (last)
    {
      last->son = first;
      first->dad=last;
      first = first->son;
      first->dad = 0;
      last = last->son;
      last->son = 0;
    }
}

//-------------------------------------------------------------------------
void NodeList::exchange(Node * pivot)
{
  if (pivot && (pivot->son)) {
      last->son  = first;
      first->dad = last;
      first      = pivot->son;
      last       = pivot;
      last->son  = 0;
      first->dad = 0;
    }
}

//-------------------------------------------------------------------------
void NodeList::shuffle(Node * p1, Node * p2)
{
  //p2 has to be downstream of p1, otherwise the data is lost!
  Node * tmp;
  if (p2 && (p2-> son)) {
      tmp        = p2->son;
      p2->son    = first;
      first->dad = p2;
      first      = p1->son;
      first->dad = 0;
      p1->son    = tmp;
      tmp->dad   = p1;
    }
}

//-------------------------------------------------------------------------
void NodeList::mixWell()
{
  if ( nbo < 2 ) return;
  Node *p1, *p2 = 0;
  
  int n = RNG.pint_exc( nbo-1 );
  
  for (p1 = first; n > 0; --n )
    p1 = p1->son;
  
  assert( p1->son );
  if ( nbo > 3 ) {
    n = RNG.pint_inc( nbo/3 ) + 1;
    for (p2 = p1; p2 && ( n > 0 ); --n)
      p2 = p2->son;
  }
  
  if ( p2 && ( p2->son ))
    shuffle(p1, p2);
  else
    exchange(p1);
}

//-------------------------------------------------------------------------
int NodeList::looksWrong() const
{
  int cnt = 0;
  Node * p = first, * q;

  while( p ) {
    if ( p->mylist != this )
      return 1;
    q = p->son;
    if ( q ) {
      if ( q->dad != p ) return 2;
    } else {
      if ( p != last ) return 3;
    }
    p = q;
    ++cnt;
  }
  
  if (cnt == nbo) return NO_ERROR;

  return 4;
}










