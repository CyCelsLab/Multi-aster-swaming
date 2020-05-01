//RCS: $Id: name_list.h,v 2.4 2005/02/02 19:27:25 nedelec Exp $

/** NameList is a class to assign Names (unsigned long int) to objects, 
    keeping a pointer to the object in an array byNames[],
    so that objects can be recalled later, by name in linear time.
    
    \author Nedelec, August 2003. EMBL Heidelberg. nedelec@embl.de
 */
#ifndef NAME_LIST_H
#define NAME_LIST_H

#include "types.h"
#include "iomessages.h"
#include "node.h"

/// NameList attribute and remember names of object. It keeps an array byNames of Ojects.
class NameList
{
  private:
    
    /// size of the chuncks in which memory is allocated
    static const int NACHUNCK = 64;
    
    /// array of objects, by their name
    Node ** byNames;    ///< we should have ( byNames[ N ]->name == N ) for any N.

    
    /// size of memory allocated
    Name    allocated;
    
    /// next name to be given
    Name    next;
    
    ///allocation function
    /** beware: after allocation, we should be able to accesss byNames[ size ]
     this is not the usual C convention, which allocate to byNames[ size-1 ]
   */
    void allocate( Name size );
 
 public:
        
    /// Constructor
    NameList();
    
    /// Destructor
    virtual ~NameList();
    
    /// next name to be given
    Name nextName() const { return next; }
  
    /// the first (smallest) name that is not used in this list
    Name firstUnusedName();
    
    /// the first (smallest) name that is used
    Name firstUsedName() const;
    
    /// the last (biggest) name that is used
    Name lastUsedName() const;
    
    /// give a name to the object, and remember it
    Name nameIt( Node * obj );
  
    /// remember a name which was not given by nameIt()
    void remember( Node * obj );
  
    /// forget a name
    void forget( Node * obj );

    /// find the object with a certain name
    Node * nodeWithThisName( Name name ) const;
  
    /// clear the name list
    void forgetAll();

    /// redistribute the name, in a continuous way
    void renameAll();
};


#endif // NAME_LIST_H
