//CVS: $Id: offscreen.h,v 1.2 2005/03/14 09:35:37 foethke Exp $


///OffScreen encapsulates two functions to open/close an OpenGL off-screen display
namespace OffScreen {

  ///creates an off-screen display of requested (width,height), or return an error code
  int open(const int width, const int height);
  
  ///close the off-screen display, releasing the memory allocated
  void close();

	
};
