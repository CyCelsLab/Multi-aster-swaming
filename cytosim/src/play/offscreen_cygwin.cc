//RCS: $Id: offscren_cygwin.cc,v 1.1 2008/01/21 17:14:00 athale Exp $
//RCS: $Id: offscreen_linux.cc,v 1.1 2005/03/13 22:29:03 nedelec Exp $
// ------------------  off-screen rendering -------------------------
// X is not just to confuse you; it supports hardware offscreen rendering.
// We need PBuffers, which are supported from GLX version 1.3 (OpenGL 1.2)
// this is C++ code, but could be translated to C, just by  putting all the
// declaration in the beginning of the function

#include <cstdio>

#ifndef __APPLE__
  #include <X11/Xlib.h>
  #include <GL/glx.h>
#endif

#include <X11/Xlib.h>
#include <GL/glx.h>


Display * dpy;
GLXPbuffer pbuf;
GLXContext glxContext;

// PBuffers are supported from OpenGL 1.3
//----------------------------------------------------------------
int OffScreen::open(const int window_width, const int window_height)
{
  dpy = XOpenDisplay( 0 );

  if ( dpy == 0 ) {
    MSG("play: off-screen rendering aborted: unable to open a connection to the X server\n");
    return 1;
  }

  if ( ! glXQueryExtension( dpy, 0, 0 ) ) {
    MSG("play: off-screen rendering aborted: glX extension not supported\n");
    return 2;
  }

#ifndef GLX_VERSION_1_3
  int minor, major;
  glXQueryVersion( dpy, &major, &minor );
  MSG("GLX version %i.%i detected\n", major, minor );
  MSG("play: off-screen rendering aborted: GLX version 1.3 or better required\n");
  return 3;
#endif

  GLint screen = XDefaultScreen( dpy );

  GLint attribList[]= {GLX_RENDER_TYPE,   GLX_RGBA_BIT,
		     GLX_DRAWABLE_TYPE, GLX_PBUFFER_BIT,
		     GLX_RED_SIZE,   4,
		     GLX_GREEN_SIZE, 4,
		     GLX_BLUE_SIZE,  4,
#if (DIM==3)
		     GLX_DEPTH_SIZE, 8,
#endif
		     None};

  /* config pointer for Frame Buffer */
  GLXFBConfig * FBConfig;
  /* number of FBConfigs returned */
  int FBConfig_Count;
  /* get frame buffer configuration */
  FBConfig = glXChooseFBConfig(dpy, screen, attribList, &FBConfig_Count);

  if (( FBConfig_Count == 0 ) || ( FBConfig == 0 )) {
    MSG("play: off-screen rendering aborted: glXChooseFBConfig returned NULL\n");
    return 4;
  }

  /* PBuffer Creation */
  int attribListP[]= {GLX_PRESERVED_CONTENTS, GL_TRUE,
          GLX_PBUFFER_WIDTH,  window_width,
		      GLX_PBUFFER_HEIGHT, window_height,
		      None };

  pbuf = glXCreatePbuffer(dpy, FBConfig[0], attribListP);

  glxContext = glXCreateNewContext( dpy, FBConfig[0], GLX_RGBA_TYPE, NULL, GL_TRUE);

  XFree( FBConfig );

  if ( glxContext == 0 ) {
    MSG("play: off-screen rendering aborted: glXCreateNewContext returned NULL\n");
    return 5;
  }

  if ( 0 == glXMakeCurrent( dpy, pbuf, glxContext ) ) {
    MSG("play: off-screen rendering aborted: cannot make the PBuffer current\n");
    return 6;
  }

  return 0;
}

//----------------------------------------------------------------
void OffScreen::close()
{
  glXDestroyContext( dpy, glxContext );
  glXDestroyPbuffer( dpy, pbuf );
  XCloseDisplay( dpy );
}
