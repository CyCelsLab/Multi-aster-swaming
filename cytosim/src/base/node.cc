//RCS: $Id: node.cc,v 2.0 2004/08/16 16:36:38 nedelec Exp $
//---------------------------------node.cc----------------------------------

#include "node.h"
#include "node_list.h"

//-------------------------------------------------------------------------
Node::Node()
{
  name   = 0; 
  mylist = 0; 
  son    = 0; 
  dad    = 0; 
  ndFlag = 0;
}

//-------------------------------------------------------------------------
Node::~Node()
{
  if( isLinked() )
    pop(); 
}

//-------------------------------------------------------------------------
void Node::push_before(Node * np)
{
  mylist->nbo++;
  np->mylist=mylist;
  np->son=this;
  np->dad=dad;
  if (dad == 0)
    mylist->first=np;
  else
    dad->son=np;
  dad=np;
}

//-------------------------------------------------------------------------
void Node::push_after(Node * np)
{
  mylist->nbo++;
  np->mylist=mylist;
  np->dad=this;
  np->son=son;
  if (son == 0)
    mylist->last=np;
  else
    son->dad=np;
  son=np;
}

//-------------------------------------------------------------------------
Node * Node::pop()
{
  assert( mylist != 0 );
  mylist->pop(this);
  return this;
}









