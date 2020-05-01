//RCS: $Id: glextensions.h,v 1.9 2005/03/30 12:30:03 nedelec Exp $

#ifndef GLEXTENSIONS_H
#define GLEXTENSIONS_H

#include "types.h"
#include "vecteur.h"

///a bunch of convenient extensions to OpenGL
namespace glExtensions
{  
  /// checks for any OpenGL error(s), reporting them with printf() on standard output
  void glPrintErrors(const char * name_of_calling_function);
  
  ///display the given text in the given color, in one corner of the display window
  void displayText(const char text[], const long color, const int window_width, const int window_height, const int position = 0 );
  
  ///save a region of the current buffer to the file, in Portable Pixel Map format
  int saveImageAsPPM(FILE * file, const int xpos, const int ypos, const int width, const int height);
  
  //------------------------------------------------------------------------------

  /// call glVertexd() or glVertexf() appropriately, depending on the type of real
  // also works for Vecteur, since they can be converted to ( real * ) automatically
  void glVertex( const real * V );

  /// call glVertex or glVertex2 depending on PSEUDO_YZ
  void glVertex( const real * V, const real pseudoY );
  
  /// call glVertex2d() or glVertex2f() appropriately, depending on the type of real
  void glVertex2( const real x, const real y);
  
  /// call glVertex3d() or glVertex3f() appropriately, depending on the type of real
  void glVertex3( const real x, const real y, const real z );
  
  /// call glTranslated() or glTranslatef() appropriately, depending on the type of real
  void glTranslate( const real x, const real y, const real z );
  
  /// call glTranslated() or glTranslatef() appropriately, depeding on DIM and the type of real
  void glTranslatev(const real* dist);

  //------------------------------------------------------------------------------

  ///set the current color from a long
  void glColor( const long color );
  
  ///set the current color from a long, dimming the color by using transparency
  void glColor( long color, const int dimm_factor );
  
  ///set the current color midway between the two colors, dimming the color by using transparency
  void glColor( const long color1, const long color2, const int dimm_factor );
  
  ///set the current color from a factor in [0,1]
  void glColorRainbow( const real h );
  
  ///set one of the current lightning color
  void glLightningColor( const GLenum face, const long color );

  ///set one of the current lightning color, with a dimm-factor
  void glLightningColor( const GLenum face, const long color, const int dimm_factor );

  ///set one of the current lightning color
  void glLightningColorOI( const long colorOutside, const long colorInside );

  ///set the current lightning color from a factor in [0,1]
  void glLightningColorRainbow( const real h );

  ///set one of the current lightning color
  void gleClearColor( const long color );
  
  
  //------------------------------------------------------------------------------

  ///set a RGB color from a factor in [0, 1], continuously varying through blue, green, red, white
  void setRainbowColor( const real h, real *r, real *g, real *b );
  
  ///conversion function from RGB to HSV color space
  void RGBtoHSV( const real r, const real g, const real b, real *h, real *s, real *v );
  ///conversion functions from HSV to RGB color space
  void HSVtoRGB( const real h, const real s, const real v, real *r, real *g, real *b );
  
  //------------------------------------------------------------------------------

};


#endif
