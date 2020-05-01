//RCS: $Id: node.h,v 2.5 2005/02/02 19:27:32 nedelec Exp $
//---------------------------------node.h----------------------------------
// Class Node is the base for most objects in the simulations.
// it is essentially contains Name, son and data, allowing the objects to be
// stored in doubly-linked list NodeList, with the following properties:
//1. Nodes can only be linked once, since the "dad" and "son" are members
//2. Nodes contain a pointer to the list they belong, which allows to 
//   unlink a node without giving a pointer to the list...
//3. Removal or addition of nodes in a NodeList is linear in time


#ifndef NODE_H
#define NODE_H

#include "types.h"                 //where type Name is defined
#include "assert_macro.h"

class NodeList;

///Node holds dad and son pointers, allowing derived objects to be linked in a NodeList
class Node
{
  ///NodeList and NameList are friend class, to access son / dad / name
  //we like friends, we can go to the movies with them on friday!
   friend class NodeList;
   friend class NameList;

protected:
     
    ///used to cleanup objects after reads
    short  ndFlag;
     
    ///used for finding objects in reads,stats...
    Name   name; 
    
    ///the list in which I am linked, or zero if I am not linked
    NodeList * mylist;
    
    ///the next object in the list, or zero if last in list
    Node * son;
    
    ///the previous object in the list, or zero if first in list
    Node * dad;
    
public:

    ///Creator
    Node();
    
    ///Destructor
    virtual ~Node();
    
    //---------------------------------------------------------
    //                  Flag functions
    //---------------------------------------------------------    
    
    ///Set the Flag
    void   setFlag(short value)     { ndFlag = value; }
    
    ///Return the Flag
    int    getFlag()          const { return ndFlag; }
    
    //---------------------------------------------------------
    //                  Name functions
    //---------------------------------------------------------

    ///Set the name
    void setName(const Name value)  { name = value; }
        
    ///Returns the name
    Name getName()            const { return name; }

    //---------------------------------------------------------
    //          Node (in List) related functions
    //---------------------------------------------------------

    ///return true if Node is linked in a list.
    bool isLinked()           const { return ( mylist != 0 ); }
    
    ///return my list in which I'm linked
    const NodeList* getList() const { return mylist; }
    
    ///return next object in the list
    Node * next()             const { return son; }
    
    ///return previous object in the list
    Node * prev()             const { return dad; }
    
    ///unlink the current object from its list (mylist)
    Node * pop();
    
    ///link another Node np before this one, in the same list
    void push_before(Node * np);
    
    ///link another Node np after this one, in the same list
    void push_after(Node * np);
};


#endif
