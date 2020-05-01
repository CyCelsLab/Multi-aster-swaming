//
// C++ Interface: globjects3d.h
//
// Description: functions for drawing objects with style 3
//
//
// CVS: $Id: globjects3d.h,v 1.5 2005/04/18 18:37:28 foethke Exp $


#ifndef GLOBJECTS3D_H
#define GLOBJECTS3D_H


#include "glut.h"
#include "types.h"


namespace glObjects3d
{

// surface for lightning
enum surface { FRONT, BACK, FRONT_AND_BACK };


//-----------------------------------------------------------------------------
//these functions are needed to draw half spheres for the pombe caps
//-----------------------------------------------------------------------------

void drawTriangle(GLdouble* v1, GLdouble* v2, GLdouble* v3, real radius);
//draw a triangle


void subdivide(GLdouble* v1, GLdouble* v2, GLdouble* v3, real radius, int depth);
//subdivide one triangle into four smaller ones


void halfSphere3D(real radius, int depth);
//draw a halfSphere in 3D including normals for lightning


void oval3D(GLUquadricObj* qobj, real length, real radius, int depth, int stacks);
//draw an oval in 3D by drawing a cylinder and two half spheres as caps
//drawing includes normals for lightning


//-----------------------------------------------------------------------------
//The following code for drawing a cylinder was taken from the original
//implementation of GLU. It was expanded to support smooth shading along
//the long axis of the cylinder.
//The sources were imported from the rpm XFree86-4.3.0-2.90.55.src.rpm
//which is included in the RedHat-Linux distribution. After unpacking
//the XFree-4.3.0.tar.gz the GLU sources can be found in
//"xc/extras/ogl-sample/main/gfx/lib/glu"

//from include/gluos.h
#if !defined(_WIN32)
#define GLAPIENTRY
#endif

//from libutil/gluint.h
#define COS cos
#define SIN sin
#define SQRT sqrt

//from libutil/quad.c
#define CACHE_SIZE 240


struct GLUExtQuadric {
    GLint 	normals;
    GLboolean	textureCoords;
    GLint	orientation;
    GLint	drawStyle;
    void	(GLAPIENTRY *errorCallback)( GLint );
};


GLUExtQuadric* GLAPIENTRY gluNewExtQuadric(void);


void gluExtQuadricError(GLUExtQuadric* qobj, GLenum which);


void GLAPIENTRY gluExtQuadricNormals(GLUExtQuadric *extqobj, GLenum normals);


void GLAPIENTRY gluExtQuadricDrawStyle(GLUExtQuadric *extqobj, GLenum drawStyle);


//modified "gluCylinder()" from libutil/quad.c
void gluEnhancedCylinder( GLUExtQuadric* qobj, GLdouble baseRadius, GLdouble topRadius, 
                          GLdouble height, GLint slices, GLint stacks,
                          unsigned long rgba1, unsigned long rgba2, surface surf=FRONT );

                          
void gluTIntersec(GLUExtQuadric* extqobj, GLdouble radius, GLdouble inflate, GLint slices);

                                                    
//modified "gluCylinder()" from libutil/quad.c
void gluBentCylinder( GLUExtQuadric* extqobj, GLdouble baseRadius, GLdouble topRadius, 
                      GLdouble height, GLdouble arcLength, GLint slices, GLint stacks,
                      unsigned long rgba1, unsigned long rgba2, surface surf=FRONT );


void banana3D(GLUExtQuadric* extqobj, real length, real capRadius, real arcLength, int depth, int stacks);
//draw an oval in 3D by drawing a cylinder and two half spheres as caps
//drawing includes normals for lightning


void tee3D(GLUExtQuadric* extqobj, GLUquadricObj* qobj, real length, real capRadius, real tJucnction, real tArmLength, real inflate, int depth, int stacks);

}

#endif // GLOBJECTS3D_H
