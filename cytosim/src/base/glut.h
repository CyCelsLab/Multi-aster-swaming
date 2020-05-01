//RCS: $Id: glut.h,v 2.0 2004/08/16 16:35:23 nedelec Exp $

// GLUT is a platform-independent windowing library built on OpenGL
// The place of the header 'glut.h' is platform-dependent, so we
// include this file "glut.h" rather than <glut.h>

#ifdef __APPLE__
    #include <GLUT/glut.h>
#else
    #include <GL/glut.h>
#endif


