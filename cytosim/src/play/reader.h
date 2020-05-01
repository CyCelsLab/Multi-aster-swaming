//RCS: $Id: reader.h,v 2.12 2005/04/18 11:41:48 foethke Exp $

#include <cstdio>

///used by play to easily find a particular frame (eg. frame 10) in cytosim's result.out
/**----------------------------------------------------------------------------
//                               SimReader
//   a helper class to read a simulation state from the file specified in IOWrapper,
//   SimReader has a buffer to remember where each state starts in the file. 
//   The only assumptions is that frames are ordered in the file.
//   Buffering and IO error handling is the main thing SimReader does.
//   It otherwise calls SIM::readState() for the input of the simulation state.
//   nedelec@embl.de,           started october 2002
//----------------------------------------------------------------------------
*/

class SimReader 
{

public:
  
  ///the lowest bound for frame index
  static const int MIN_FRAME = 0;
  
  ///higher bound for frame index
  static const int MAX_FRAME = 8192;

private:

  ///array containing the starting position in the file for each frame
  unsigned long framestart[ MAX_FRAME ]; 

  ///first frame giving eof
  int eofFrame;
  
  ///first known good frame
  int minFrame;
  
  ///last known good frame
  int maxFrame; 

public:

  //These function should be self-explanatory
  //they are used for e.g. in play, to help decide
  //if we are viewing the last frame, to rewind in case loop=1
      
  ///first frame known to be readible
  int firstGoodFrame()       const { return minFrame;     }
  
  ///last frame known to be good
  int lastGoodFrame()        const { return maxFrame;     }
  
  ///first frame leading to and end-of-file message
  int firstEofFrame()        const { return eofFrame;     }
  
  ///last frame before the EOF (not necessarily readible/good)
  int lastFrameBeforeEof()   const { return eofFrame - 1; }
  
  ///clear the buffer
  void clearBuffer();
  
  ///open file for input
  int openFile(const char * filename);  
 
  ///close the current input file
  void closeFile();

  /// return the input file
  FILE *  getInputFile();
  
  ///dummy constructor 
  SimReader() {
    clearBuffer(); 
  }
  
  ///constructor with the name of input file <result.out>
  SimReader(const char * filename) {
    openFile( filename );
  }

  ///dummy destructor
  virtual ~SimReader() {
    //we do not close the input file, this is done in IOWrapper IO
  }
  
  ///return IO.eof()
  bool eof() const;
  
  ///go to the start of the file
  void rewindFile();

  ///find the starting point of frame by brute force !
  int seekFrame( int frame );
 
  /// read the frame specified by an integer:
  /// returns NO_ERROR for success, an error code, or throws an exception
  /// the frame read is stored in the global variable 'sim'
  int readFrameNoCatch( int frame ); 

  /// read the frame specified:
  /// returns NO_ERROR for success, an error code, does not throw exceptions
  /// the frame read is stored in the global variable 'sim'
  int readFrame( int frame ); 
  
  /// read the next frame from the current position in the file
  int readNextFrameNoCatch();
  
  /// read the next frame from the current position in the file
  int readNextFrame();

};


