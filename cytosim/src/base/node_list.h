//RCS: $Id: node_list.h,v 2.1 2004/12/06 18:00:39 nedelec Exp $

#ifndef NODE_LIST_H
#define NODE_LIST_H

#include "node.h"
#include "name_list.h"

/// Class NodeList is a doubly linked (up and down) list to hold variable number of Node.
/**
Linking or Unlinking of Nodes is done in constant time. The List also
keeps track of how many objects are linked, and has some functions to 
randomize the order of the Nodes in the list.
 
Note: Nodes are named by the class NameList when they are linked into the list. 
The NameList keeps track of assigned Names, and can be used to find objects
by their names. Several NodeList can share the same NameList.
*/
class NodeList
{
  /// Node is a friend, to allow fast linking / unlinking
  friend class Node;
  
 protected:

  /// First Node of the list
  Node * first;
    
  /// Last Node of the list
  Node * last;
    
  ///number of objects in the list
  int nbo;

  ///NameList when the names of the nodes are being given/recorded 
  NameList * namelist;


 public:

  /// Constructor initializing all members to zero
  NodeList()                      { nbo=0; last=0; first=0; namelist=0;  }
    
  /// constructor giving a name list
  NodeList(NameList * ns)         { nbo=0; last=0; first=0; namelist=ns; }

  /// Destructor
  virtual ~NodeList()             { erase(); }


  /// First Node
  Node * getFront()  const        { return first; }
    
  /// Last Node
  Node * getBack()   const        { return last; }
    
  /// Number of objects in the list
  int size()      const           { return nbo; }

  /// put new Node at the begining of the list
  Node * pushFront(Node * np);
    
  /// put new Node at the end of the list
  Node * pushBack(Node * np);

  /// Remove Node op from the list
  Node * pop(Node * op);
  
  /// Remove First Node  from the list
  Node * popFront()               { return pop(first); }
  
  /// Remove First Node  from the list
  Node * popBack()                { return pop(last); }
  
  /// remove np from its list and the put it at at the end of *this
  void popAndPush(Node * np);

  /// Clear the list by resetting first and last
  void clearFast()               { first = 0; last = 0; nbo = 0; }
  
  /// clear the list by calling pop(first) until empty
  virtual void clear();
  
  /// clear the list as above, calling delete( ) for each node
  virtual void erase();

  /// delete objects in the list which flag is different than flagToKeep
  void deleteIfFlagDifferentThan( int flagToKeep );

  /// Give the a NameList for the NodeList
  void setNameList(NameList * ns) { namelist=ns; }
  
  /// Getting the NameList 
  NameList * getNameList()  const { return namelist; }
    
  /// Return the next name on the list
  Name nextName()        const    { return namelist->nextName(); }

  /// search in the list node with specific name
  Node * nodeWithThisName(Name)      const;
  
  /// search without using the namelist array, ie. linear scanning
  Node * nodeWithThisName_scan(Name) const;

  /// testing coherence in list (dad, son relation or belonging to same list)
  int  looksWrong() const;

  /// put the last on top of the list
  void turn();
  
  /// divide the list in two parts starting from pivot node and exchange them
  void exchange(Node * pivot);
  
  /// divide the list in three parts (first-p1-p2-last) and move p1-p2 in front
  void shuffle(Node * p1, Node * p2);
  
  /// randomly using exchange and shuffle function.  
  void mixWell();

};


#endif
