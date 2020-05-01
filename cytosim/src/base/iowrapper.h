//RCS: $Id: iowrapper.h,v 2.26 2005/04/18 10:22:21 nedelec Exp $

#ifndef  IOWRAPPER_H
#define  IOWRAPPER_H

#include "types.h"
#include <cstdio>
#include <cstdarg>


/// standard name for the output file
static const char RESULT_OUT[] = "result.out";

/// standard name for the input file
const char START_IN[] = "start.in";


///Allows input/output for cytosim's coordinate files such as [result.out]
/** class IOwrapper is made of two parts
simulation state output to a file 'output' by write()
simulation state input from a file 'input' by read()
inputs supports automatic swapping of bytes, for cross-platform compatibility
*/
class IOWrapper
{
  
public:

  /// standard types used in the input/output formats
  typedef          char        int8;
  typedef unsigned char       uint8;
  typedef          short       int16;
  typedef unsigned short      uint16;
  typedef          int         int32;
  typedef unsigned int        uint32;
  typedef          long long   int64;
  typedef unsigned long long  uint64;
  
  ///Constructor base which prepares infile/outfile to be read/written  - sets default values
  void IOWrapperBase();
  
  ///Default constructor, just calls ConstructorBase
  IOWrapper();
  
  ///Destructor which closes the input and output files
  virtual ~IOWrapper();
      
private:
  /// file to write to
  FILE * infile;

  /// The format ID of the input: this allow backward compatibility with old formats
  int inputFormat;   
    
  /// The dimensionality of its vectors: a 2D state can be read from a 3D simulation
  int inputDIM;     
    
/** if the state is stored in a binary format, inputBinary
    is set to 1 or 2. with 2, byte order is swapped automatically
    this occurs for example when reading a simulation calculated 
    on PC from mac, or vice et versa.
*/
  int inputBinary;  
    
  /// Necessary for converting between different binary formats
  void    swapBytes(char*, int);

public:
    
  /// return the input file
  FILE *  getInputFile()                    { return infile; }
	
  /// opens file to be read, by default in binary format (issue in Windows)
  int     openInputFile(const char name[] = START_IN, const bool binary = true);
    
	/// File to be read is set
  void    setInputFile(FILE * f);
	
	/// Makes sure the input file exists and has no errors
  void    checkInputFile();
	
	/// Closes input file if it exists
  void    closeInputFile();
    
	///Sets dimentionnality of vectors
  void    setInputDIM(const int d)       { inputDIM = d; }
	
	/// returns dimentionnality
  int     getInputDIM()       const      { return inputDIM; }
	
  /// returns the type of input
  int     getFileFormat()     const      { return inputFormat; }
	
	/// Returns 1 for native binary format, 2 for non-native binary format, and 0 if not binary
  int     getInputBinary()    const      { return inputBinary; }
	
	/// Returns the End of File of the input file
  bool    eof()               const      { return (infile == 0) || feof(infile); }
	
	///Returns true if the infile is not open or has errors according to ferror()
  bool    error()             const      { return (infile == 0) || ferror(infile); }
	
	///returns current value of file position indicator
  unsigned long   getPos()    const      { return ftell(infile); }
	
	///Sets file position indicator in the input file to position pos, relative to the start of the file
  void    seek(long pos)                 { fseek( infile, pos, SEEK_SET ); }
     
	///Checks way bytes are stored in the binary format
  void    readBinarySignature(FILE * = 0);
	
	///Sets the input format using the input file
  void    readInputFormat(FILE * = 0);
	
	///Sets the vector dimension using the input file
  void    readInputDIM(FILE * = 0);
	
	//Reads a line from the input file. If last argument is false, reads until EOF.
  int     readLine(char *, const int size, FILE * = 0, const bool stopAtNewLine = true );
	
	///skips the next input character of the file, default is the input file.
  void    skipChar(FILE * = 0);
	
	///skips the current line of the file, default is the input file
  void    skipLine(FILE * = 0);
	
	///Searches (incorrectly) in the input file for the string c
  bool    skipUntil(const char *, FILE * = 0);  
    
	///Reads and returns a character from the input file
  char    readEntry(FILE * = 0);
	
	///Reads the identifier (2 chars) that should be at the start of every record. return true if EOF
  int     readRecordTag(char ch[2], FILE * = 0);
	
	///Reads an int 
  int     readStringInt(FILE * = 0);
	
	///Reads a float (ASCII)
  float   readFloatAscii(FILE * = 0);
	
  ///Reads an unsigned short on 8 bits
  char    readCharAscii(FILE * = 0);
  
	///Reads an unsigned int on 8 bits
  uint8   readUInt8(FILE * = 0);
	
	///Reads an unsigned int on 8 bits
  int8    readInt8(FILE * = 0);
	
	///Reads an unsigned int on 16 bits
  int16   readInt16(FILE * = 0);
	
	///Reads a unsigned int on 4 bytes
  uint16  readUInt16(FILE * = 0);
	
	///Reads anunsigned int on 8 bytes
  uint32  readUInt32(FILE * = 0);
	
	///Reads a float (correcting for different binary formats)
  float   readReal32(FILE * = 0);
	
	///Reads an array of floats 
  void    readReal32(real[], const int, FILE * = 0);
	
	///Reads a vector 
  void    readReal32Vect(real *, FILE * = 0);
    
        
        
  /** PART 2: File for simulation state ouput: access by write???()*/
	
private:
		
  ///a file to write to
  FILE * outfile;
    
  // Modifiers for simulation output:
  bool writeBinary;
  
  /** Write mode:
    Normal (0): all the frame are output in the same file [result.out]
           (1): same, but the file is closed and opened for each frame.
    Cut    (2): each frame is output in a separate file [resultXXXX.out]  
    */
  int writeMode;

public:
  /// returns the output file
  FILE *  getOutputFile()                   { return outfile; }
	
	///Opens the output file and prepares for writing
  void    openOutputFile(const char * name = RESULT_OUT, const int append=0);
	
	///Opens the output file and prepares for writing
  void    openOutputFile(const int append);
	
	///Sets the file to write to. 
  void    setOutputFile(FILE * f);
	
	/// Checks if the ouput file exists, and that there are no errors
  void    checkOutputFile() const;
	
	///Closes the output file
  void    closeOutputFile();
	
  ///Sets to write in binary format
  void    produceBinaryOutput(bool m)    { writeBinary = m; }
	  
	///Sets the write mode, determining how the frames are stored
  void    setWriteMode(int m)            { writeMode = m; }
	
	///Gets the write mode
  int     getWriteMode()      const      { return writeMode; }
  
  ///Write the file format version
  void    writeFileFormatVersion(const int, FILE * = 0);
  
  ///Puts a tag to specify a binary file, and the byte order 
  void    writeBinarySignature(FILE * = 0);
	
	///Inserts a return, only in text output mode, for nicer output
  void    writeFormattingNewLine(FILE * out = 0);
	
	///Inserts N space(s), only in text output mode, for nicer output
  void    writeFormattingSpace(const int N = 1, FILE * out = 0);
	
  ///Writes a unsigned char on 8 bit
  void    writeCharAscii(const char, FILE * = 0);
  
	///Writes an unsigned short 
  void    writeUInt8(const int, FILE * = 0);
	
	///Writes an unsigned short 
  void    writeInt8(const int, FILE * = 0);
	
	///Writes an unsigned int 
  void    writeInt16(const int, FILE * = 0);

	///Writes two unsigned int 
  void    writeInt16(const int, const int, FILE * = 0);
	
	///Writes three unsigned int 
  void    writeInt16(const int, const int, const int, FILE * = 0);
  
  ///Writes an unsigned int on 2 bytes = 16 bits, used for Name
  void    writeUInt16(const uint32, FILE * = 0);
	
	///Writes an unsigned int on 4 bytes = 32 bits, used for LongName
  void    writeUInt32(const uint32, FILE * = 0);
	  
	///Writes the real as a float, i.e. on 4 bytes = 32 bits
  void    writeReal32(const real, FILE * = 0);
	
	///Writes an array of float, of given size, each on 32 bits
  void    writeReal32Vect(const int n, const real *, const bool formattingNewLineBefore = false, FILE * = 0);
    
  ///Writes the newline and two chars to identify a new record in the file
  void    writeRecordTag(const char[2], FILE * = 0);

	///Writes an int (ASCII)
  void    writeIntAscii(const int, FILE * = 0);
	
	///Writes 2 ints (ASCII)
  void    writeIntAscii(const int, const int, FILE * = 0);
	
	///Writes a float (ASCII)
  void    writeRealAscii(const real, FILE * = 0);
	
	///Writes a string 
  void    writeString(const char *, FILE * = 0);
  
  ///raw printf on IO's file. use the fprintf() format.
  void    print(const char * fmt, ...);

  ///Flush the IO stream
  void    flush( FILE * = 0 );
  
  
  ///returns true if the first two letters of line match the tag id[]
  inline static bool compareRecordTag(const char * line, const char * id) {
    return ((line[0] == id[0]) && (line[1] == id[1]));
  }
  
  
  ///set the directoryPath that correspond to the given file
  int  resolvePath(char * directoryPath, const int size, const char * fileName);

  ///get the full name for the file, from the given directoryPath
  int  fixPath(char * fileName, const int size, const char * directoryPath);
};

///the global instantiation used for input/output
extern IOWrapper IO;

#endif
