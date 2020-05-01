//RCS: $Id: object_list.h,v 1.3 2005/04/12 20:47:18 nedelec Exp $

#ifndef OBJECTLIST_H
#define OBJECTLIST_H

#include "object.h"
#include "node_list.h"

///A list of Object. Not very useful, as we define more specialized lists
class ObjectList : public NodeList
{
  
public:
  ///Creator
  ObjectList() {}
  
  ///Destructor
  virtual ~ObjectList() {}
  
  ///return the first object in the list
  Object * getFront() const { return static_cast<Object*>( first ); }
  
  ///return the last object in the list
  Object * getBack()  const { return static_cast<Object*>( last ); }

};

#endif
