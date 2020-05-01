//
// C++ Implementation: globjects3d.cc
//
// Description: functions for drawing objects with style 3
//
//
// CVS: $Id: globjects3d.cc,v 1.8 2005/04/19 22:26:47 foethke Exp $


#include "globjects3d.h"
#include "smath.h"
#include "iomessages.h"
#include "glextensions.h"
#include "play_param.h"


using namespace glExtensions;


//-----------------------------------------------------------------------------
//these functions are needed to draw half spheres for the pombe caps
//-----------------------------------------------------------------------------

void glObjects3d::drawTriangle(GLdouble* v1, GLdouble* v2, GLdouble* v3, real radius)
//draw a triangle
{
  glBegin(GL_TRIANGLES);
    for(int ii = 0; ii < 3; ii++) {
      v3[ii] /= radius;
    }
    glNormal3dv(v3);
    for(int ii = 0; ii < 3; ii++) {
      v3[ii] *= radius;
    }
    glVertex3dv(v3);
    for(int ii = 0; ii < 3; ii++) {
      v2[ii] /= radius;
    }
    glNormal3dv(v2);
    for(int ii = 0; ii < 3; ii++) {
      v2[ii] *= radius;
    }
    glVertex3dv(v2);
    for(int ii = 0; ii < 3; ii++) {
      v1[ii] /= radius;
    }
    glNormal3dv(v1);
    for(int ii = 0; ii < 3; ii++) {
      v1[ii] *= radius;
    }
    glVertex3dv(v1);
  glEnd();
}


void glObjects3d::subdivide(GLdouble* v1, GLdouble* v2, GLdouble* v3, real radius, int depth)
//subdivide one triangle into four smaller ones
{
  GLdouble v12[3], v23[3], v31[3];
  real norm12=0.;
  real norm23=0.;
  real norm31=0.;
  
  if(depth == 0) {
    drawTriangle(v1, v2, v3, radius);
    return;
  }
  
  //calculate the new vertices and normalize them
  for(int ii = 0; ii < 3; ii++) {
    v12[ii] = (v1[ii]+v2[ii])/2.;
    v23[ii] = (v2[ii]+v3[ii])/2.;
    v31[ii] = (v3[ii]+v1[ii])/2.;
    norm12 += v12[ii]*v12[ii];
    norm23 += v23[ii]*v23[ii];
    norm31 += v31[ii]*v31[ii];
  }
  norm12 = sqrt(norm12);
  norm23 = sqrt(norm23);
  norm31 = sqrt(norm31);
  for(int ii = 0; ii < 3; ii++) {
    v12[ii] *= radius/norm12;
    v23[ii] *= radius/norm23;
    v31[ii] *= radius/norm31;
  }
  
  subdivide(v1,  v12, v31, radius, depth-1);
  subdivide(v2,  v23, v12, radius, depth-1);
  subdivide(v3,  v31, v23, radius, depth-1);
  subdivide(v12, v23, v31, radius, depth-1);
}


void glObjects3d::halfSphere3D(real radius, int depth)
//draw a halfSphere in 3D including normals for lightning
{
  //we start drawing the half sphere as a tetraeder
  GLdouble v1[3] = {0, radius*1.0,0};
  GLdouble v2[3] = {0,-radius*0.5,  radius*sqrt(0.75)};
  GLdouble v3[3] = {0,-radius*0.5, -radius*sqrt(0.75)};
  GLdouble v4[3] = {-radius,0,0};
  
  // now we subdevide the triangles "depth" times
  subdivide(v1, v2, v4, radius, depth);
  subdivide(v3, v4, v2, radius, depth);
  subdivide(v3, v1, v4, radius, depth);
}


void glObjects3d::oval3D(GLUquadricObj* qobj, real length, real radius, int depth, int stacks)
//draw an oval in 3D by drawing a cylinder and two half spheres as caps
//drawing includes normals for lightning
{
  char cylDepth = 1<<depth;
  
  glPushMatrix();
  glTranslate(-length,0,0);
  halfSphere3D(radius, depth);
  glRotated(90, 0, 1., 0);
  //gluEnhancedCylinder(qobj, radius, radius, 2.*length, 3*cylDepth, stacks, PP.pombecolor, PP.nuenvcolor);
  gluCylinder(qobj, radius, radius, 2.*length, 3*cylDepth, stacks);
  glRotated(90, 0,1., 0);
  glTranslate(-2.*length,0,0);
  halfSphere3D(radius, depth);
  glPopMatrix();
}


//-----------------------------------------------------------------------------
//The following code for drawing a cylinder was taken from the original
//implementation of GLU. It was expanded to support smooth shading along
//the long axis of the cylinder.
//The sources were imported from the rpm XFree86-4.3.0-2.90.55.src.rpm
//which is included in the RedHat-Linux distribution. After unpacking
//the XFree-4.3.0.tar.gz the GLU sources can be found in
//"xc/extras/ogl-sample/main/gfx/lib/glu"


glObjects3d::GLUExtQuadric* GLAPIENTRY glObjects3d::gluNewExtQuadric(void)
{
    GLUExtQuadric *newstate;

    newstate = (GLUExtQuadric *) new GLUExtQuadric;
    if (newstate == 0) {
        /* Can't report an error at this point... */
        return 0;
    }
    newstate->normals = GLU_SMOOTH;
    newstate->textureCoords = GL_FALSE;
    newstate->orientation = GLU_OUTSIDE;
    newstate->drawStyle = GLU_FILL;
    newstate->errorCallback = 0;
    return newstate;
}


void glObjects3d::gluExtQuadricError(GLUExtQuadric* qobj, GLenum which)
{
  if (qobj->errorCallback) {
    qobj->errorCallback(which);
  }
}


void GLAPIENTRY glObjects3d::gluExtQuadricNormals(GLUExtQuadric *extqobj, GLenum normals)
{
    switch (normals) {
      case GLU_SMOOTH:
      case GLU_FLAT:
      case GLU_NONE:
        break;
      default:
        gluExtQuadricError(extqobj, GLU_INVALID_ENUM);
        return;
    }
    extqobj->normals = normals;
}


void GLAPIENTRY glObjects3d::gluExtQuadricDrawStyle(GLUExtQuadric *extqobj, GLenum drawStyle)
{
    switch(drawStyle) {
      case GLU_POINT:
      case GLU_LINE:
      case GLU_FILL:
      case GLU_SILHOUETTE:
        break;
      default:
        gluExtQuadricError(extqobj, GLU_INVALID_ENUM);
        return;
    }
    extqobj->drawStyle = drawStyle;
}


void glObjects3d::gluEnhancedCylinder( GLUExtQuadric* qobj, GLdouble baseRadius, GLdouble topRadius, 
                                       GLdouble height, GLint slices, GLint stacks,
                                       unsigned long rgba1, unsigned long rgba2, surface surf )
{
    GLint i, j;
    GLfloat sinCache[CACHE_SIZE];
    GLfloat cosCache[CACHE_SIZE];
    GLfloat sinCache2[CACHE_SIZE];
    GLfloat cosCache2[CACHE_SIZE];
    GLfloat sinCache3[CACHE_SIZE];
    GLfloat cosCache3[CACHE_SIZE];
    GLfloat angle;
    GLfloat zLow, zHigh;
    GLfloat sintemp, costemp;
    GLfloat length;
    GLfloat deltaRadius;
    GLfloat zNormal;
    GLfloat xyNormalRatio;
    GLfloat radiusLow, radiusHigh;
    int needCache2, needCache3;

    // allocate memory for the enhanced color cache
    static GLfloat* colorCache;
    static int allocated = 0;
    if( allocated < stacks + 1 ) {
      if ( colorCache ) delete[] colorCache;
      allocated  = stacks + 1;
      colorCache = new GLfloat[allocated*4];
    }

    if (slices >= CACHE_SIZE) slices = CACHE_SIZE-1;

    if (slices < 2 || stacks < 1 || baseRadius < 0.0 || topRadius < 0.0 ||
	    height < 0.0) {
	gluExtQuadricError(qobj, GLU_INVALID_VALUE);
	return;
    }

    /* Compute length (needed for normal calculations) */
    deltaRadius = baseRadius - topRadius;
    length = SQRT(deltaRadius*deltaRadius + height*height);
    if (length == 0.0) {
	gluExtQuadricError(qobj, GLU_INVALID_VALUE);
	return;
    }

    /* Cache is the vertex locations cache */
    /* Cache2 is the various normals at the vertices themselves */
    /* Cache3 is the various normals for the faces */
    needCache2 = needCache3 = 0;
    if (qobj->normals == GLU_SMOOTH) {
	needCache2 = 1;
    }

    if (qobj->normals == GLU_FLAT) {
	if (qobj->drawStyle != GLU_POINT) {
	    needCache3 = 1;
	}
	if (qobj->drawStyle == GLU_LINE) {
	    needCache2 = 1;
	}
    }

    zNormal = deltaRadius / length;
    xyNormalRatio = height / length;

    for (i = 0; i < slices; i++) {
	angle = 2 * PI * i / slices;
	if (needCache2) {
	    if (qobj->orientation == GLU_OUTSIDE) {
		sinCache2[i] = xyNormalRatio * SIN(angle);
		cosCache2[i] = xyNormalRatio * COS(angle);
	    } else {
		sinCache2[i] = -xyNormalRatio * SIN(angle);
		cosCache2[i] = -xyNormalRatio * COS(angle);
	    }
	} 
	sinCache[i] = SIN(angle);
	cosCache[i] = COS(angle);
    }

    if (needCache3) {
	for (i = 0; i < slices; i++) {
	    angle = 2 * PI * (i-0.5) / slices;
	    if (qobj->orientation == GLU_OUTSIDE) {
		sinCache3[i] = xyNormalRatio * SIN(angle);
		cosCache3[i] = xyNormalRatio * COS(angle);
	    } else {
		sinCache3[i] = -xyNormalRatio * SIN(angle);
		cosCache3[i] = -xyNormalRatio * COS(angle);
	    }
	}
    } 

    sinCache[slices] = sinCache[0];
    cosCache[slices] = cosCache[0];
    if (needCache2) {
	sinCache2[slices] = sinCache2[0];
	cosCache2[slices] = cosCache2[0];
    }
    if (needCache3) {
	sinCache3[slices] = sinCache3[0];
	cosCache3[slices] = cosCache3[0];
    }

    //calculate color values for all stacks
    unsigned char rgbaChar1[4];
    unsigned char rgbaChar2[4];
    rgbaChar1[0] = ( rgba1 >> 24 ) & 0xFF;
    rgbaChar1[1] = ( rgba1 >> 16 ) & 0xFF;
    rgbaChar1[2] = ( rgba1 >> 8  ) & 0xFF;
    rgbaChar1[3] =   rgba1         & 0xFF;
    rgbaChar2[0] = ( rgba2 >> 24 ) & 0xFF;
    rgbaChar2[1] = ( rgba2 >> 16 ) & 0xFF;
    rgbaChar2[2] = ( rgba2 >> 8  ) & 0xFF;
    rgbaChar2[3] =   rgba2         & 0xFF;
    GLfloat rstep = (GLfloat)(rgbaChar2[0] - rgbaChar1[0])/(255.*stacks);
    GLfloat gstep = (GLfloat)(rgbaChar2[1] - rgbaChar1[1])/(255.*stacks);
    GLfloat bstep = (GLfloat)(rgbaChar2[2] - rgbaChar1[2])/(255.*stacks);
    GLfloat astep = (GLfloat)(rgbaChar2[3] - rgbaChar1[3])/(255.*stacks);
    for (i = 0; i <= stacks; i++) {
	colorCache[i*4    ] = (GLfloat)rgbaChar1[0]/255. + i*rstep;
	colorCache[i*4 + 1] = (GLfloat)rgbaChar1[1]/255. + i*gstep;
	colorCache[i*4 + 2] = (GLfloat)rgbaChar1[2]/255. + i*bstep;
	colorCache[i*4 + 3] = (GLfloat)rgbaChar1[3]/255. + i*astep;
    }
    glEnable(GL_COLOR_MATERIAL);
    switch (surf) {
      case FRONT:
        glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);
        break;
      case BACK:
        glColorMaterial(GL_BACK, GL_AMBIENT_AND_DIFFUSE);
        break;
      case FRONT_AND_BACK:
        glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
        break;
      default:
        MSG.error("gluEnhancedCylinder", "Invalid color surface.");
        break;
    }

    switch (qobj->drawStyle) {
      case GLU_FILL:
	/* Note:
	** An argument could be made for using a TRIANGLE_FAN for the end
	** of the cylinder of either radii is 0.0 (a cone).  However, a 
	** TRIANGLE_FAN would not work in smooth shading mode (the common 
	** case) because the normal for the apex is different for every
	** triangle (and TRIANGLE_FAN doesn't let me respecify that normal).
	** Now, my choice is GL_TRIANGLES, or leave the GL_QUAD_STRIP and
	** just let the GL trivially reject one of the two triangles of the
	** QUAD.  GL_QUAD_STRIP is probably faster, so I will leave this code
	** alone.
	*/
	for (j = 0; j < stacks; j++) {
	    zLow = j * height / stacks;
	    zHigh = (j + 1) * height / stacks;
	    radiusLow = baseRadius - deltaRadius * ((float) j / stacks);
	    radiusHigh = baseRadius - deltaRadius * ((float) (j + 1) / stacks);

	    glBegin(GL_QUAD_STRIP);
	    for (i = 0; i <= slices; i++) {
		switch(qobj->normals) {
		  case GLU_FLAT:
		    glNormal3f(sinCache3[i], cosCache3[i], zNormal);
		    break;
		  case GLU_SMOOTH:
		    glNormal3f(sinCache2[i], cosCache2[i], zNormal);
		    break;
		  case GLU_NONE:
		  default:
		    break;
		}
		if (qobj->orientation == GLU_OUTSIDE) {
		    if (qobj->textureCoords) {
			glTexCoord2f(1 - (float) i / slices,
				(float) j / stacks);
		    }
		    glColor4fv(colorCache+(j*4));
		    glVertex3f(radiusLow * sinCache[i], 
			    radiusLow * cosCache[i], zLow);
		    if (qobj->textureCoords) {
			glTexCoord2f(1 - (float) i / slices,
				(float) (j+1) / stacks);
		    }
		    glColor4fv(colorCache+((j+1)*4));
		    glVertex3f(radiusHigh * sinCache[i], 
			    radiusHigh * cosCache[i], zHigh);
		} else {
		    if (qobj->textureCoords) {
			glTexCoord2f(1 - (float) i / slices,
				(float) (j+1) / stacks);
		    }
		    glVertex3f(radiusHigh * sinCache[i], 
			    radiusHigh * cosCache[i], zHigh);
		    if (qobj->textureCoords) {
			glTexCoord2f(1 - (float) i / slices,
				(float) j / stacks);
		    }
		    glVertex3f(radiusLow * sinCache[i], 
			    radiusLow * cosCache[i], zLow);
		}
	    }
	    glEnd();
	}
	break;
      case GLU_POINT:
	glBegin(GL_POINTS);
	for (i = 0; i < slices; i++) {
	    switch(qobj->normals) {
	      case GLU_FLAT:
	      case GLU_SMOOTH:
		glNormal3f(sinCache2[i], cosCache2[i], zNormal);
		break;
	      case GLU_NONE:
	      default:
		break;
	    }
	    sintemp = sinCache[i];
	    costemp = cosCache[i];
	    for (j = 0; j <= stacks; j++) {
		zLow = j * height / stacks;
		radiusLow = baseRadius - deltaRadius * ((float) j / stacks);

		if (qobj->textureCoords) {
		    glTexCoord2f(1 - (float) i / slices,
			    (float) j / stacks);
		}
		glVertex3f(radiusLow * sintemp, 
			radiusLow * costemp, zLow);
	    }
	}
	glEnd();
	break;
      case GLU_LINE:
	for (j = 1; j < stacks; j++) {
	    zLow = j * height / stacks;
	    radiusLow = baseRadius - deltaRadius * ((float) j / stacks);

	    glBegin(GL_LINE_STRIP);
	    for (i = 0; i <= slices; i++) {
		switch(qobj->normals) {
		  case GLU_FLAT:
		    glNormal3f(sinCache3[i], cosCache3[i], zNormal);
		    break;
		  case GLU_SMOOTH:
		    glNormal3f(sinCache2[i], cosCache2[i], zNormal);
		    break;
		  case GLU_NONE:
		  default:
		    break;
		}
		if (qobj->textureCoords) {
		    glTexCoord2f(1 - (float) i / slices,
			    (float) j / stacks);
		}
		glVertex3f(radiusLow * sinCache[i], 
			radiusLow * cosCache[i], zLow);
	    }
	    glEnd();
	}
	/* Intentionally fall through here... */
      case GLU_SILHOUETTE:
	for (j = 0; j <= stacks; j += stacks) {
	    zLow = j * height / stacks;
	    radiusLow = baseRadius - deltaRadius * ((float) j / stacks);

	    glBegin(GL_LINE_STRIP);
	    for (i = 0; i <= slices; i++) {
		switch(qobj->normals) {
		  case GLU_FLAT:
		    glNormal3f(sinCache3[i], cosCache3[i], zNormal);
		    break;
		  case GLU_SMOOTH:
		    glNormal3f(sinCache2[i], cosCache2[i], zNormal);
		    break;
		  case GLU_NONE:
		  default:
		    break;
		}
		if (qobj->textureCoords) {
		    glTexCoord2f(1 - (float) i / slices,
			    (float) j / stacks);
		}
		glVertex3f(radiusLow * sinCache[i], radiusLow * cosCache[i], 
			zLow);
	    }
	    glEnd();
	}
	for (i = 0; i < slices; i++) {
	    switch(qobj->normals) {
	      case GLU_FLAT:
	      case GLU_SMOOTH:
		glNormal3f(sinCache2[i], cosCache2[i], 0.0);
		break;
	      case GLU_NONE:
	      default:
		break;
	    }
	    sintemp = sinCache[i];
	    costemp = cosCache[i];
	    glBegin(GL_LINE_STRIP);
	    for (j = 0; j <= stacks; j++) {
		zLow = j * height / stacks;
		radiusLow = baseRadius - deltaRadius * ((float) j / stacks);

		if (qobj->textureCoords) {
		    glTexCoord2f(1 - (float) i / slices,
			    (float) j / stacks);
		}
		glVertex3f(radiusLow * sintemp, 
			radiusLow * costemp, zLow);
	    }
	    glEnd();
	}
	break;
      default:
	break;
    }
    glDisable(GL_COLOR_MATERIAL);
}


void glObjects3d::gluTIntersec(GLUExtQuadric* extqobj, GLdouble radius, GLdouble inflate, GLint slices)
{
  GLdouble sinCache[CACHE_SIZE];
  GLdouble cosCache[CACHE_SIZE];
  GLdouble yCache[2*CACHE_SIZE];
  GLdouble zCache[2*CACHE_SIZE];
  GLdouble xNormalCache[2*CACHE_SIZE];
  GLdouble yNormalCache[2*CACHE_SIZE];
  GLdouble zNormalCache[2*CACHE_SIZE];
  GLdouble angle;
  GLdouble inflRad;
  GLdouble stackstep;
  GLdouble xLow, xHigh;
  GLdouble sing, cosg;
  GLint    i, j;
  
  if (slices >= CACHE_SIZE) slices = CACHE_SIZE-1;
//  if (stacks >= CACHE_SIZE) stacks = CACHE_SIZE-1;
    
  if (slices < 2 || radius < 0.0 || inflate > 0.0) {
    gluExtQuadricError(extqobj, GLU_INVALID_VALUE);
    return;
  }
  
  /* Cache is the vertex locations cache */
  for( i = 0; i < slices; i++ ) {
    angle = 2 * PI * i / slices;
    sinCache[i] = SIN(angle);
    cosCache[i] = COS(angle);
  }

  sinCache[slices] = sinCache[0];
  cosCache[slices] = cosCache[0];

  /* we need these to be exactly zero */
  cosCache[slices/4]   = 0;
  cosCache[slices*3/4] = 0;
  sinCache[slices/2]   = 0;
  
  inflRad = radius + inflate;
  
  //glLightningColor( GL_FRONT_AND_BACK, 0x777777FF );
  
  //load first array
  int indx = 0;
  for( i = slices/4; i <= slices*3/4; i++ ) {
    xNormalCache[indx] = 0.;
    yNormalCache[indx] = cosCache[i];
    zNormalCache[indx] = sinCache[i];
    yCache[indx]       = inflRad * cosCache[i]; 
    zCache[indx]       = inflRad * sinCache[i];
    indx++;
  }
  
  //load next arrays alternating
  bool bit = 0;
  for( j = 1; j <= slices/4; j++ ) {
    GLdouble angleLow  = (j - 1)*PI/2/(slices/4);
    GLdouble angleHigh =  j     *PI/2/(slices/4);
    xLow               = (radius+inflate/sqrt(2.))*sin(angleLow);
    xHigh              = (radius+inflate/sqrt(2.))*sin(angleHigh);
    
    //alternate between the two buffers
    bit ? bit = 0 : bit = 1;
    
    indx = bit*CACHE_SIZE;
    for( i = slices/4-j; i <= slices*3/4+j; i++ ) {
      
      if( cosCache[i] > fabs(xHigh)/radius ) {
        //calculate cos(gamma)
	sing = -(radius - xHigh/cosCache[i])/inflate;
	//printf("sing: %f cosCache[i]: %f\n", sing, cosCache[i]);
	cosg = 1.0 - sing*sing;
	cosg > 0 ? cosg=sqrt(cosg) : cosg = 0;
// 	xNormalCache[indx] = 0.;
// 	yNormalCache[indx] = (radius + inflate*cosg) * fabs(cosCache[i]);
// 	zNormalCache[indx] = inflRad * sinCache[i];
	GLdouble normalNorm =   inflRad*inflRad*inflate*inflate*sinCache[i]*sinCache[i]*sinCache[i]*sinCache[i] \
	  + (radius*radius*inflate*inflate*(2*sing*cosg) + inflate*inflate*inflate*inflate)*sinCache[i]*sinCache[i]*cosCache[i]*cosCache[i];
	//normalNorm = sqrt(normalNorm);
	//normalNorm = 1;
  	xNormalCache[indx] = -inflRad*inflate*sing*cosCache[i]*fabs(cosCache[i]) / normalNorm;
  	yNormalCache[indx] = -inflRad*inflate*cosg*cosCache[i]*cosCache[i] / normalNorm;
  	zNormalCache[indx] = -((radius+inflate*sing)*inflate*sing*sinCache[i]*fabs(cosCache[i]) + \
	                      (radius+inflate*cosg)*inflate*cosg*fabs(sinCache[i])*cosCache[i]) / normalNorm;
	//	printf("x: %f, y: %f, z: %f\n", xNormalCache[indx], yNormalCache[indx], zNormalCache[indx]);
	yCache[indx]       = (radius + inflate*cosg) * fabs(cosCache[i]);
	zCache[indx]       = inflRad * sinCache[i];
      } else {
	xNormalCache[indx] = 0.;
	yNormalCache[indx] = cosCache[i];
        zNormalCache[indx] = sinCache[i];
        yCache[indx]       = inflRad * cosCache[i];
        zCache[indx]       = inflRad * sinCache[i];
      }
      indx++;
    }
    
    //    if( ! bit )
      //      glLightningColor( GL_FRONT_AND_BACK, 0xFF0000FF );
      //      glLightningColorOI( 0xFF0000FF, 0x00000000 );
      //    else
      //      glLightningColor( GL_FRONT_AND_BACK, 0x00FF00FF );
      //      glLightningColorOI( 0x00FF00FF, 0x00000000 );
    
    int kk;
    int lindx = ((bit == 0)*CACHE_SIZE);
    int hindx = bit*CACHE_SIZE;
    //    printf("Bit & 0: %d\n", (bit == 0));
    //the right side of the base cylinder
    glBegin(GL_TRIANGLE_STRIP);
    glNormal3d( xNormalCache[hindx], yNormalCache[hindx],  zNormalCache[hindx]);
    glVertex3d( xHigh,               yCache      [hindx],  zCache      [hindx]);
    for( kk = 0; kk <= slices/4+j-1; kk++ ) {
      glNormal3d( xNormalCache[lindx+kk],   yNormalCache[lindx+kk],    zNormalCache[lindx+kk]);
      glVertex3d( xLow,                     yCache      [lindx+kk],    zCache      [lindx+kk]);
      glNormal3d( xNormalCache[hindx+kk+1], yNormalCache[hindx+kk+1],  zNormalCache[hindx+kk+1]);
      glVertex3d( xHigh,                    yCache      [hindx+kk+1],  zCache      [hindx+kk+1]);
    }
    for( kk = slices/4+j-2; kk >= 0; kk-- ) {
      glNormal3d( xNormalCache[lindx+kk],   yNormalCache[lindx+kk],   -zNormalCache[lindx+kk]);
      glVertex3d( xLow,                     yCache      [lindx+kk],   -zCache      [lindx+kk]);
      glNormal3d( xNormalCache[hindx+kk+1], yNormalCache[hindx+kk+1], -zNormalCache[hindx+kk+1]);
      glVertex3d( xHigh,                    yCache      [hindx+kk+1], -zCache      [hindx+kk+1]);
    }
    glNormal3d( xNormalCache[hindx], yNormalCache[hindx], -zNormalCache[hindx]);
    glVertex3d( xHigh,               yCache      [hindx], -zCache      [hindx]);
    glEnd();

    //the left side of the base cylinder
    glBegin(GL_TRIANGLE_STRIP);
    glNormal3d( -xNormalCache[hindx], yNormalCache[hindx], -zNormalCache[hindx]);
    glVertex3d( -xHigh,               yCache      [hindx], -zCache      [hindx]);
    for( kk = 0; kk <= slices/4+j-1; kk++ ) {
      glNormal3d( -xNormalCache[lindx+kk],   yNormalCache[lindx+kk],   -zNormalCache[lindx+kk]);
      glVertex3d( -xLow,                     yCache      [lindx+kk],   -zCache      [lindx+kk]);
      glNormal3d( -xNormalCache[hindx+kk+1], yNormalCache[hindx+kk+1], -zNormalCache[hindx+kk+1]);
      glVertex3d( -xHigh,                    yCache      [hindx+kk+1], -zCache      [hindx+kk+1]);
    }
    for( kk = slices/4+j-2; kk >= 0; kk-- ) {
      glNormal3d( -xNormalCache[lindx+kk],   yNormalCache[lindx+kk],    zNormalCache[lindx+kk]);
      glVertex3d( -xLow,                     yCache      [lindx+kk],    zCache      [lindx+kk]);
      glNormal3d( -xNormalCache[hindx+kk+1], yNormalCache[hindx+kk+1],  zNormalCache[hindx+kk+1]);
      glVertex3d( -xHigh,                    yCache      [hindx+kk+1],  zCache      [hindx+kk+1]);
    }
    glNormal3d( -xNormalCache[hindx], yNormalCache[hindx], zNormalCache[hindx]);
    glVertex3d( -xHigh,               yCache      [hindx], zCache      [hindx]);
    glEnd();
    
    //the arm section is a mirror image of the base with the
    //symmetry axis being the line y=-x

    //the front of the arm (z>0)
    glBegin(GL_TRIANGLE_STRIP);
    glNormal3d( -yNormalCache[hindx], xNormalCache[hindx], zNormalCache[hindx]);
    glVertex3d( -yCache      [hindx], xHigh,               zCache      [hindx]);
    for( kk = 0; kk <= j-1; kk++ ) {
      glNormal3d( -yNormalCache[lindx+kk],   xNormalCache[lindx+kk],   zNormalCache[lindx+kk]);
      glVertex3d( -yCache      [lindx+kk],   xLow,                     zCache      [lindx+kk]);
      glNormal3d( -yNormalCache[hindx+kk+1], xNormalCache[hindx+kk+1], zNormalCache[hindx+kk+1]);
      glVertex3d( -yCache      [hindx+kk+1], xHigh,                    zCache      [hindx+kk+1]);
    }
    for( kk = j-2; kk >= 0; kk-- ) {
      glNormal3d(  yNormalCache[lindx+kk],   xNormalCache[lindx+kk],   zNormalCache[lindx+kk]);
      glVertex3d(  yCache      [lindx+kk],   xLow,                     zCache      [lindx+kk]);
      glNormal3d(  yNormalCache[hindx+kk+1], xNormalCache[hindx+kk+1], zNormalCache[hindx+kk+1]);
      glVertex3d(  yCache      [hindx+kk+1], xHigh,                    zCache      [hindx+kk+1]);
    }
    glNormal3d(  yNormalCache[hindx], xNormalCache[hindx], zNormalCache[hindx]);
    glVertex3d(  yCache      [hindx], xHigh,               zCache      [hindx]);
    glEnd();

    //the back of the arm (z<0)
    glBegin(GL_TRIANGLE_STRIP);
    glNormal3d(  yNormalCache[hindx], xNormalCache[hindx],    -zNormalCache[hindx]);
    glVertex3d(  yCache      [hindx], xHigh,                  -zCache      [hindx]);
    for( kk = 0; kk <= j-1; kk++ ) {
      glNormal3d(  yNormalCache[lindx+kk],   xNormalCache[lindx+kk],   -zNormalCache[lindx+kk]);
      glVertex3d(  yCache      [lindx+kk],   xLow,                     -zCache      [lindx+kk]);
      glNormal3d(  yNormalCache[hindx+kk+1], xNormalCache[hindx+kk+1], -zNormalCache[hindx+kk+1]);
      glVertex3d(  yCache      [hindx+kk+1], xHigh,                    -zCache      [hindx+kk+1]);
    }
    for( kk = j-2; kk >= 0; kk-- ) {
      glNormal3d( -yNormalCache[lindx+kk],   xNormalCache[lindx+kk],   -zNormalCache[lindx+kk]);
      glVertex3d( -yCache      [lindx+kk],   xLow,                     -zCache      [lindx+kk]);
      glNormal3d( -yNormalCache[hindx+kk+1], xNormalCache[hindx+kk+1], -zNormalCache[hindx+kk+1]);
      glVertex3d( -yCache      [hindx+kk+1], xHigh,                    -zCache      [hindx+kk+1]);
    }
    glNormal3d( -yNormalCache[hindx], xNormalCache[hindx], -zNormalCache[hindx]);
    glVertex3d( -yCache      [hindx], xHigh,               -zCache      [hindx]);
    glEnd();
  }
  
  //
  //===========================================================================
  //
  
  GLint stacks=5;

  //load next arrays alternating
  bit = 1;
  for( j = 0; j <= stacks; j++ ) {
    xLow  = radius+inflate/sqrt(2.) - (j-1)*inflate/sqrt(2.)/stacks;
    xHigh = radius+inflate/sqrt(2.) -  j   *inflate/sqrt(2.)/stacks;
    
    //alternate between the two buffers
    bit ? bit = 0 : bit = 1;
    
    indx = bit*CACHE_SIZE;
    for( i = 0; i <= slices; i++ ) {
      
      if( cosCache[i] > fabs(xHigh)/radius ) {
        //calculate cos(gamma)
	sing = -(radius - xHigh/cosCache[i])/inflate;
	cosg = 1.0 - sing*sing;
	cosg > 0 ? cosg=sqrt(cosg) : cosg = 0;
// 	xNormalCache[indx] = 0.;
// 	yNormalCache[indx] = (radius + inflate*cosg) * fabs(cosCache[i]);
// 	zNormalCache[indx] = inflRad * sinCache[i];
	GLdouble normalNorm =   inflRad*inflRad*inflate*inflate*sinCache[i]*sinCache[i]*sinCache[i]*sinCache[i] \
	  + (radius*radius*inflate*inflate*(2*sing*cosg) + inflate*inflate*inflate*inflate)*sinCache[i]*sinCache[i]*cosCache[i]*cosCache[i];
	//normalNorm = sqrt(normalNorm);
	normalNorm = 1;
  	xNormalCache[indx] = -inflRad*inflate*sing*cosCache[i]*fabs(cosCache[i]) / normalNorm;
  	yNormalCache[indx] = -inflRad*inflate*cosg*cosCache[i]*cosCache[i] / normalNorm;
  	zNormalCache[indx] = -((radius+inflate*sing)*inflate*sing*sinCache[i]*fabs(cosCache[i]) + \
	                      (radius+inflate*cosg)*inflate*cosg*fabs(sinCache[i])*cosCache[i]) / normalNorm;
	//	printf("x: %f, y: %f, z: %f\n", xNormalCache[indx], yNormalCache[indx], zNormalCache[indx]);
	yCache[indx]       = (radius + inflate*cosg) * fabs(cosCache[i]);
	zCache[indx]       = inflRad * sinCache[i];
      } else {
	xNormalCache[indx] = 0.;
	yNormalCache[indx] = cosCache[i];
        zNormalCache[indx] = sinCache[i];
        yCache[indx]       = inflRad * cosCache[i];
        zCache[indx]       = inflRad * sinCache[i];
      }
      indx++;
    }

    if( j > 0 ) {

//       if( ! bit )
// 	// glLightningColor( GL_FRONT_AND_BACK, 0xFF0000FF );
// 	glLightningColorOI( 0xFF0000FF, 0x00000000 );
//       else
// 	// glLightningColor( GL_FRONT_AND_BACK, 0x00FF00FF );
// 	glLightningColorOI( 0x00FF00FF, 0x00000000 );
    
      int kk;
      int lindx = ((bit == 0)*CACHE_SIZE);
      int hindx = bit*CACHE_SIZE;
      //the right side of the base cylinder
      glBegin(GL_TRIANGLE_STRIP);
      for( kk = 0; kk <= slices/2; kk++ ) {
	glNormal3d( xNormalCache[hindx+kk], yNormalCache[hindx+kk],  zNormalCache[hindx+kk]);
	glVertex3d( xHigh,                  yCache      [hindx+kk],  zCache      [hindx+kk]);
	glNormal3d( xNormalCache[lindx+kk], yNormalCache[lindx+kk],  zNormalCache[lindx+kk]);
	glVertex3d( xLow,                   yCache      [lindx+kk],  zCache      [lindx+kk]);
      }
      for( kk = slices/2-1; kk >= 0; kk-- ) {
	glNormal3d( xNormalCache[hindx+kk], yNormalCache[hindx+kk], -zNormalCache[hindx+kk]);
	glVertex3d( xHigh,                  yCache      [hindx+kk], -zCache      [hindx+kk]);
	glNormal3d( xNormalCache[lindx+kk], yNormalCache[lindx+kk], -zNormalCache[lindx+kk]);
	glVertex3d( xLow,                   yCache      [lindx+kk], -zCache      [lindx+kk]);
      }
      glEnd();

      //the left side of the base cylinder
      glBegin(GL_TRIANGLE_STRIP);
      for( kk = 0; kk <= slices/2; kk++ ) {
	glNormal3d( -xNormalCache[hindx+kk], yNormalCache[hindx+kk], -zNormalCache[hindx+kk]);
	glVertex3d( -xHigh,                  yCache      [hindx+kk], -zCache      [hindx+kk]);
	glNormal3d( -xNormalCache[lindx+kk], yNormalCache[lindx+kk], -zNormalCache[lindx+kk]);
	glVertex3d( -xLow,                   yCache      [lindx+kk], -zCache      [lindx+kk]);
      }
      for( kk = slices/2-1; kk >= 0; kk-- ) {
	glNormal3d( -xNormalCache[hindx+kk], yNormalCache[hindx+kk],  zNormalCache[hindx+kk]);
	glVertex3d( -xHigh,                  yCache      [hindx+kk],  zCache      [hindx+kk]);
	glNormal3d( -xNormalCache[lindx+kk], yNormalCache[lindx+kk],  zNormalCache[lindx+kk]);
	glVertex3d( -xLow,                   yCache      [lindx+kk],  zCache      [lindx+kk]);
      }
      glEnd();

      //the front of the arm (z>0)
      glBegin(GL_TRIANGLE_STRIP);
      for( kk = 0; kk <= slices/4; kk++ ) {
	glNormal3d( -yNormalCache[hindx+kk], xNormalCache[hindx+kk],  zNormalCache[hindx+kk]);
	glVertex3d( -yCache      [hindx+kk], xHigh                 ,  zCache      [hindx+kk]);
	glNormal3d( -yNormalCache[lindx+kk], xNormalCache[lindx+kk],  zNormalCache[lindx+kk]);
	glVertex3d( -yCache      [lindx+kk], xLow,                    zCache      [lindx+kk]);
      }
      for( kk = slices/4-1; kk >= 0; kk-- ) {
	glNormal3d(  yNormalCache[hindx+kk], xNormalCache[hindx+kk],  zNormalCache[hindx+kk]);
	glVertex3d(  yCache      [hindx+kk], xHigh,                   zCache      [hindx+kk]);
	glNormal3d(  yNormalCache[lindx+kk], xNormalCache[lindx+kk],  zNormalCache[lindx+kk]);
	glVertex3d(  yCache      [lindx+kk], xLow,                    zCache      [lindx+kk]);
      }
      glEnd();

      //the back of the arm (z>0)
      glBegin(GL_TRIANGLE_STRIP);
      for( kk = 0; kk <= slices/4; kk++ ) {
	glNormal3d(  yNormalCache[hindx+kk], xNormalCache[hindx+kk], -zNormalCache[hindx+kk]);
	glVertex3d(  yCache      [hindx+kk], xHigh                 , -zCache      [hindx+kk]);
	glNormal3d(  yNormalCache[lindx+kk], xNormalCache[lindx+kk], -zNormalCache[lindx+kk]);
	glVertex3d(  yCache      [lindx+kk], xLow,                   -zCache      [lindx+kk]);
      }
      for( kk = slices/4-1; kk >= 0; kk-- ) {
	glNormal3d( -yNormalCache[hindx+kk], xNormalCache[hindx+kk], -zNormalCache[hindx+kk]);
	glVertex3d( -yCache      [hindx+kk], xHigh,                  -zCache      [hindx+kk]);
	glNormal3d( -yNormalCache[lindx+kk], xNormalCache[lindx+kk], -zNormalCache[lindx+kk]);
	glVertex3d( -yCache      [lindx+kk], xLow,                   -zCache      [lindx+kk]);
      }
      glEnd();
    }
  }

//   for( j = 0; j < stacks; j++ ) {
//     xLow  = -(radius+inflate/sqrt(2.)) +  j   *inflate/sqrt(2.)/stacks;
//     xHigh = -(radius+inflate/sqrt(2.)) + (j+1)*inflate/sqrt(2.)/stacks;
    
//     if( fmod(j,2) )
//       //      glLightningColor( GL_FRONT_AND_BACK, 0xFF0000FF );
//       glLightningColorOI( 0xFF0000FF, 0x00000000 );
//     else
//       //      glLightningColor( GL_FRONT_AND_BACK, 0x00FF00FF );
//       glLightningColorOI( 0x00FF00FF, 0x00000000 );

//     glBegin(GL_TRIANGLE_STRIP);
//     for( i = slices; i >= 0; i-- ) {
    
//       if( cosCache[i] > fabs(xHigh)/radius ) {
//         //calculate cos(gamma)
// 	cosg = 1.0 - (radius + xHigh/cosCache[i])/inflate*(radius + xHigh/cosCache[i])/inflate;
// 	cosg > 0 ? cosg=sqrt(cosg) : cosg = 0;
// 	glNormal3f(0.,
// 		   (radius + inflate*cosg) * fabs(cosCache[i]),
// 		   inflRad * sinCache[i]);
// 	glVertex3f(xHigh,
// 		   (radius + inflate*cosg) * fabs(cosCache[i]),
// 		   inflRad * sinCache[i]);
//       } else {
//         glNormal3d(0.,
//                    inflRad * cosCache[i],
//                    inflRad * sinCache[i]);
//         glVertex3d(xHigh,
//                    inflRad * cosCache[i],
//                    inflRad * sinCache[i]);
//       }
//       if( cosCache[i] > fabs(xLow)/radius ) {
//         //calculate cos(gamma)
// 	cosg = 1.0 - (radius + xLow/cosCache[i])/inflate*(radius + xLow/cosCache[i])/inflate;
// 	cosg > 0 ? cosg=sqrt(cosg) : cosg = 0;
// 	glNormal3f(0.,
// 		   (radius + inflate*cosg) * fabs(cosCache[i]),
// 		   inflRad * sinCache[i]);
// 	glVertex3f(xLow,
// 		   (radius + inflate*cosg) * fabs(cosCache[i]),
// 		   inflRad * sinCache[i]);
//       } else {
//         glNormal3d(0.,
//                    inflRad * cosCache[i],
//                    inflRad * sinCache[i]);
//         glVertex3d(xLow,
//                    inflRad * cosCache[i], 
//                    inflRad * sinCache[i]);
//       }
//     }
//     glEnd();
//   }
}



void glObjects3d::gluBentCylinder( GLUExtQuadric* extqobj, GLdouble baseRadius, GLdouble topRadius, 
                                   GLdouble height, GLdouble arcLength, GLint slices, GLint stacks,
                                   unsigned long rgba1, unsigned long rgba2, surface surf )
{
  GLint i, j;
  GLfloat sinCache[CACHE_SIZE];
  GLfloat cosCache[CACHE_SIZE];
  GLfloat bentSinCache[CACHE_SIZE]; // cache for the bending angles: Sin(alpha)
  GLfloat bentCosCache[CACHE_SIZE]; // cache for the bending angles: Cos(alpha)
  GLfloat angle;
  GLfloat zLow, zHigh;
  GLfloat deltaRadius, slope;
  GLfloat radiusLow, radiusHigh;
  GLfloat rad, alpha;              // polar coordinates of the cylinders surface points
  GLfloat bentRadius;              // the radius of the banana
  GLfloat normNormal;              // the norm of surface normals
            
  // allocate memory for the enhanced color cache
  static GLfloat* colorCache;
  static int allocated = 0;
  if( allocated < stacks + 1 ) {
    if ( colorCache ) delete[] colorCache;
    allocated  = stacks + 1;
    colorCache = new GLfloat[allocated*4];
  }

  if (slices >= CACHE_SIZE) slices = CACHE_SIZE-1;
  if (stacks >= CACHE_SIZE) stacks = CACHE_SIZE-1;
    
  if (slices < 2 || stacks < 1 || baseRadius < 0.0 || topRadius < 0.0 || height < 0.0) {
    gluExtQuadricError(extqobj, GLU_INVALID_VALUE);
    return;
  }

  /* Compute delta radius */
  deltaRadius = topRadius - baseRadius;

  /* Calculate the slope of the bent cylinder as if it was not bent.
     This is needed for normal calculations later on. If you don't want
     correct normals (e.g. to avoid visible corners in aubergines,
     set slope to 0. */
//  slope = deltaRadius/(2.*arcLength);
  slope = 0;

  /* Calculate the bent radius of the banana */
  bentRadius = height/(2.*arcLength);
    
  /* Cache is the vertex locations cache */
  for (i = 0; i < slices; i++) {
    angle = 2 * PI * i / slices;
    sinCache[i] = SIN(angle);
    cosCache[i] = COS(angle);
  }

  sinCache[slices] = sinCache[0];
  cosCache[slices] = cosCache[0];

  for(i = 0; i <= stacks; i++) {
    alpha = (i * height / stacks) / bentRadius;
    bentSinCache[i] = SIN(alpha);
    bentCosCache[i] = COS(alpha);
  }
  
  //calculate color values for all stacks
  unsigned char rgbaChar1[4];
  unsigned char rgbaChar2[4];
  rgbaChar1[0] = ( rgba1 >> 24 ) & 0xFF;
  rgbaChar1[1] = ( rgba1 >> 16 ) & 0xFF;
  rgbaChar1[2] = ( rgba1 >> 8  ) & 0xFF;
  rgbaChar1[3] =   rgba1         & 0xFF;
  rgbaChar2[0] = ( rgba2 >> 24 ) & 0xFF;
  rgbaChar2[1] = ( rgba2 >> 16 ) & 0xFF;
  rgbaChar2[2] = ( rgba2 >> 8  ) & 0xFF;
  rgbaChar2[3] =   rgba2         & 0xFF;
  GLfloat rstep = (GLfloat)(rgbaChar2[0] - rgbaChar1[0])/(255.*stacks);
  GLfloat gstep = (GLfloat)(rgbaChar2[1] - rgbaChar1[1])/(255.*stacks);
  GLfloat bstep = (GLfloat)(rgbaChar2[2] - rgbaChar1[2])/(255.*stacks);
  GLfloat astep = (GLfloat)(rgbaChar2[3] - rgbaChar1[3])/(255.*stacks);
  for (i = 0; i <= stacks; i++) {
    colorCache[i*4    ] = (GLfloat)rgbaChar1[0]/255. + i*rstep;
    colorCache[i*4 + 1] = (GLfloat)rgbaChar1[1]/255. + i*gstep;
    colorCache[i*4 + 2] = (GLfloat)rgbaChar1[2]/255. + i*bstep;
    colorCache[i*4 + 3] = (GLfloat)rgbaChar1[3]/255. + i*astep;
  }
  glEnable(GL_COLOR_MATERIAL);
  switch (surf) {
    case FRONT: {
      glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);
      break; }
    case BACK:
      glColorMaterial(GL_BACK, GL_AMBIENT_AND_DIFFUSE);
      break;
    case FRONT_AND_BACK:
      glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
      break;
    default:
      MSG.error("gluBentCylinder", "Invalid color surface.");
      break;
  }

  switch (extqobj->drawStyle) {
    case GLU_FILL:
      /* Note:
      ** An argument could be made for using a TRIANGLE_FAN for the end
      ** of the cylinder of either radii is 0.0 (a cone).  However, a 
      ** TRIANGLE_FAN would not work in smooth shading mode (the common 
      ** case) because the normal for the apex is different for every
      ** triangle (and TRIANGLE_FAN doesn't let me respecify that normal).
      ** Now, my choice is GL_TRIANGLES, or leave the GL_QUAD_STRIP and
      ** just let the GL trivially reject one of the two triangles of the
      ** QUAD.  GL_QUAD_STRIP is probably faster, so I will leave this code
      ** alone.
      */
      for (j = 0; j < stacks; j++) {
        zLow       = j * height / stacks;
        zHigh      = (j + 1) * height / stacks;
        radiusLow  = baseRadius + deltaRadius * ((float) j / stacks);
        radiusHigh = baseRadius + deltaRadius * ((float) (j + 1) / stacks);

        glBegin(GL_QUAD_STRIP);
        for (i = 0; i <= slices; i++) {

          /* This is disabled. We always calculate smooth normals. */
          switch(extqobj->normals) {
            case GLU_FLAT:
              break;
            case GLU_SMOOTH:
              break;
            case GLU_NONE:
            default:
              break;
          }
          
          if (extqobj->orientation == GLU_OUTSIDE) {
            rad   = bentRadius + radiusLow*cosCache[i];
            if (extqobj->textureCoords)
              glTexCoord2f(1 - (float)i / slices, (float)j / stacks);
            glColor4fv(colorCache+(j*4));

            normNormal = sqrt(slope*slope + rad*rad);            
            glNormal3f(( rad*sinCache[i])/normNormal,
                       ( slope*bentSinCache[j] + rad*cosCache[i]*bentCosCache[j])/normNormal,
                       (-slope*bentCosCache[j] + rad*cosCache[i]*bentSinCache[j])/normNormal);
            glVertex3f( radiusLow*sinCache[i],
                        rad*bentCosCache[j] - bentRadius,
                        rad*bentSinCache[j]);

            rad   = bentRadius + radiusHigh*cosCache[i];            
            if (extqobj->textureCoords)
              glTexCoord2f(1 - (float)i / slices, (float)(j+1) / stacks);
            glColor4fv(colorCache+((j+1)*4));

            normNormal = sqrt(slope*slope + rad*rad);            
            glNormal3f(( rad*sinCache[i])/normNormal,
                       ( slope*bentSinCache[j+1] + rad*cosCache[i]*bentCosCache[j+1])/normNormal,
                       (-slope*bentCosCache[j+1] + rad*cosCache[i]*bentSinCache[j+1])/normNormal);
            glVertex3f( radiusHigh*sinCache[i],
                        rad*bentCosCache[j+1] - bentRadius,
                        rad*bentSinCache[j+1]);
          } else {
            MSG.error("glBentCylinder:","orientations->GLU_OUTSIDE is not yet implemented!");
          }
        }
        glEnd();
      }
      break;
    case GLU_POINT:
      MSG.error("glBentCylinder:","drawStyle->GLU_POINT is not yet implemented!");
      break;
    case GLU_LINE:
      MSG.error("glBentCylinder:","drawStyle->GLU_LINE is not yet implemented!");
      /* Intentionally fall through here... */
    case GLU_SILHOUETTE:
      MSG.error("glBentCylinder:","drawStyle->GLU_SILHOUETTE is not yet implemented!");
      break;
    default:
      break;
    }
    glDisable(GL_COLOR_MATERIAL);
}


void glObjects3d::banana3D(GLUExtQuadric* extqobj, real length, real capRadius, real arcLength, int depth, int stacks)
//draw a banana in 3D by drawing a bent cylinder and two half spheres as caps
//drawing includes normals for lightning
{  
  char cylDepth    = 1<<depth;
  GLdouble bRadius = length/arcLength;
  
  glPushMatrix();
  glTranslate( bRadius*sin(-arcLength), bRadius*cos(-arcLength) - bRadius, 0 );
  glRotated(arcLength/PI*180, 0, 0, 1.);
  halfSphere3D(capRadius, depth);
  glRotated(90, 0, 1., 0);
  //set the "normal" color to invisible (alpha=00)
  //This has no immediate effect, but since only the front color will be set
  //with glColorMaterial in the bent cylinder, the inside of the banana will
  //be displayed in the color that was last set with any glColor command.
  glColor( 0x00000000 );
  gluBentCylinder(extqobj, capRadius, capRadius, 2.*length, arcLength, 3*cylDepth, stacks, PP.boxcolor, PP.boxcolor, FRONT);
  glPopMatrix();
  
  glPushMatrix();
  glTranslate( -bRadius*sin(-arcLength), bRadius*cos(-arcLength) - bRadius, 0 );
  glRotated((PI-arcLength)/PI*180, 0, 0, 1.);
  halfSphere3D(capRadius, depth);
  glPopMatrix();
}


void glObjects3d::tee3D(GLUExtQuadric* extqobj, GLUquadricObj* qobj, real length, real capRadius, real junction, real armLength, real inflate, int depth, int stacks)
//draw a T in 3D by drawing three cylinders and three half spheres for the arms
//and the intersection area
//drawing includes normals for lightning
{
  //cylDepth must be a multiple of 4
  if( depth < 2 ) {
    MSG.error("glObjects3d::tee3D", "depth must be greater than 1!");
    exit(0);
  }
  
  char cylDepth = 1<<depth;
  real llength  = length+(junction-capRadius);
  real rlength  = length-(junction+capRadius);
  real radius   = capRadius + inflate;
  
  //cylDepth must be a multiple of 4 and
  //3*cylDepth must fit exactly into a certain part of the intersection area (dist)
//  real dist = capRadius + fabs(inflate)*(1./sqrt(2) - 1.);
//  fmod(dist )
  
  //left arm
  glPushMatrix();
  glTranslate(-length, 0, 0);
  halfSphere3D(radius, depth);
  glRotated(90, 0, 1., 0);
  if( llength > 0 ) gluCylinder(qobj, radius, radius, llength, 3*cylDepth, stacks);
  glPopMatrix();
  
  //right arm
  glPushMatrix();
  glTranslate(length, 0, 0);
  glRotated(180, 0, 1., 0);
  halfSphere3D(radius, depth);
  glRotated(90, 0, 1., 0);
  if( rlength > 0 ) gluCylinder(qobj, radius, radius, rlength, 3*cylDepth, stacks);
  glPopMatrix();

  //t arm
  glPushMatrix();
  glTranslate(junction, armLength+capRadius, 0);
  glRotated(270, 0, 0, 1.);
  halfSphere3D(radius, depth);
  glRotated(90, 0, 1., 0);
  if( armLength > 0 ) gluCylinder(qobj, radius, radius, armLength, 3*cylDepth, stacks);
  glPopMatrix();

  //intersection area
  glPushMatrix();
  glTranslate(junction, 0, 0);
  gluTIntersec(extqobj, capRadius, inflate, 3*cylDepth);
  glPopMatrix();
}
