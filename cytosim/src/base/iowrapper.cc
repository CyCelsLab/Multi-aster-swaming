//RCS: $Id: iowrapper.cc,v 2.33 2005/04/08 10:47:12 nedelec Exp $



#include <cstdlib>

#include <sys/param.h>

#include <cstring>

#include <cctype>

#include <libgen.h>



#include "iowrapper.h"

#include "iomessages.h"



#include "assert_macro.h"

#include "main.h"

#include "smath.h"

#include "exceptions.h"





IOWrapper IO;



//------------------------------------------------------------------------



void IOWrapper::IOWrapperBase()

{

    infile  = 0;

    outfile = 0;

    

    //by default, we assume that the input file is going to be of the same

    //version than the current OUTPUT file format. This is only useful in 

    //the case where the input file does not have the fileFormatVersion info.

    //we could for all other cases set here inputFormat = 0;

    inputFormat  = fileFormatVersion;

    inputBinary  = 0;

    inputDIM     = DIM;

    

    writeBinary  = true;

    writeMode    = 1;

    

    //check the size of the type, as we rely on them to write byte-by-byte

    if( sizeof( uint16 ) != 2 ) MSG.error("IOWrapper", "sizeof( uint16 ) != 2");

    if( sizeof( uint32 ) != 4 ) MSG.error("IOWrapper", "sizeof( uint32 ) != 4");

    if( sizeof( uint64 ) != 8 ) MSG.error("IOWrapper", "sizeof( uint64 ) != 8");

    if( sizeof( float  ) != 4 ) MSG.error("IOWrapper", "sizeof( float ) != 4");

}





IOWrapper::IOWrapper()

{

  IOWrapperBase();

}



IOWrapper::~IOWrapper()

{

  closeInputFile();

  closeOutputFile();

}







//============================================================================

//=========                     INPUT                              ===========

//============================================================================





//------------------------------------------------------------------------

void IOWrapper::setInputFile( FILE * file )

{

  closeInputFile();

  if ( file == 0 )

    throw IOException("IOWrapper::setInputFile() file = 0");

  if ( ferror( file ) )

    throw IOException("IOWrapper::setInputFile() ferror(file) returned true");

  

  infile = file;

}



//------------------------------------------------------------------------

/** Sets to read only, and checks if format is binary or not. 

    Default path is given as START_IN. 

    Checks the input file if it is open and has no errors.

*/

int IOWrapper::openInputFile(const char name[], const bool binary)

{

  if ( strlen( name ) <= 0 )

    throw IOException("IOWrapper::openInputFile() empty name");

  

  if ( infile ) 

    closeInputFile();

  

  if ( binary )

    infile = fopen(name, "rb");

  else

    infile = fopen(name, "r");

  

  if ( infile == 0 ) return 1;

    //throw IOException("IOWrapper::openInputFile() file not found");

  if ( ferror( infile ) ) return 2;

    //throw IOException("IOWrapper::openInputFile() file opened with error");

  return NO_ERROR;

}



//------------------------------------------------------------------------

void IOWrapper::closeInputFile()

{

  if ( infile && ( infile != stdin ))

    if ( fclose( infile ) )

      throw IOException("IOWrapper::closeInputFile() fclose() return non-zero");

  infile = 0;

}



//------------------------------------------------------------------------

void IOWrapper::checkInputFile()

{

    if ( infile == 0 )

        throw IOException("IOWrapper:: input file == 0");

    if ( ferror( infile ) )

        throw IOException("IOWrapper:: ferror( input file )");

}





//------------------------------------------------------------------------

void IOWrapper::skipChar(FILE * in)

{

    if (in == 0) in = infile;

    fgetc( in );

}



//------------------------------------------------------------------------

/** The default file read from is the input file. 

    If the input file is a binary file, skipline reads 1 byte at a time 

    but discards the information.

*/

void IOWrapper::skipLine(FILE * in)

{

    if (in == 0) in = infile;

    char ch;

    do {

        

        if ( inputBinary )

            fread( &ch, 1, 1, in );

        else

            ch = getc(in);

        

        if ( ferror(in) )

            throw IOException("IOWrapper::skipLine() failed");

        

    }  while ( !feof(in) && (ch != '\n'));

}



//------------------------------------------------------------------------

/** 

Readline reads a line input file (default) and stores it in the buffer 'line' 

until the end of the line of length 'size', (or until the end of the file). 

*/

int IOWrapper::readLine(char * line, const int size, FILE * in, const bool stopAtNewLine)

{

    if (in == 0) in = infile;

    char ch;

    int c = 0;

    

    do {

        if ( inputBinary )

            fread( &ch, 1, 1, in );

        else

            ch = getc(in);



        if ( ferror(in) ) 

          throw IOException("IOWrapper::readLine() failed");



        if (( ch == '\n' ) && ( stopAtNewLine ))

            break;

        

        if ( c < size-1 )

          line[ c++ ] = ch;

        else

          break;

    }  while ( !feof(in) );

    

    line[ c ] = '\0';

    return c;

}



//------------------------------------------------------------------------

bool IOWrapper::skipUntil(const char * c, FILE * in)

/* First the line is scanned looking for a match for the first character in c, and then checks if

all the characters match until the space. If they all match, it returns true, else false.*/



/** \todo: skipUntil will fail if the pattern has internal repetitions eg. 'aab'. Correct this. 

*/

{

    if (in == 0) in = infile;

    char ch;

    int d = 0;

    

    do {

        

        if ( inputBinary )

            fread( &ch, 1, 1, in );

        else

            ch = getc(in);

        

        if ( ch == c[ d ] ) {

            ++d;

            if ( c[ d ] == '\0' ) {

                fseek( in, -strlen(c), SEEK_CUR );  //to the start of string

                return true;

            }

        } else { 

			d = 0;

        }

        if ( ferror(in) )

			throw IOException("IOWrapper::skipUntil() failed");  

    }  while( !feof(in) );

    

    return false;

}





//------------------------------------------------------------------------

/** Swaps bytes of the first and nth, and the second and n+1th character. 

used to convert between little-endian and big-endian binary files */



void IOWrapper::swapBytes( char * c, int n )

{

    char   y = c[ n-1 ];

    c[ n-1 ] = c[ 0 ];

    c[ 0 ]   = y;

    

    if ( n > 3 ) {

        y = c[ n-2 ];

        c[ n-2 ] = c[ 1 ];

        c[ 1 ]   = y;

    }

}





//------------------------------------------------------------------------

/**writes a short of value 1 and compares if the input file stores bytes in the same way.

If same-endian = 1, opposite-endian = 2

*/

void IOWrapper::readBinarySignature(FILE * in)

{

    if (in == 0) in = infile;

    char ch[4];

    fread( ch, 1, 2, in);

    *( (uint16*) (ch+2) ) = (uint16) 1;//make a short of value 1

    inputBinary = 1 + ( ch[0] != ch[2] );//same-endian = 1, opposite-endian = 2

    MSG(8, "BINARY TAG (%i %i) native (%i %i)\n", ch[0], ch[1], ch[2], ch[3]);

}





//------------------------------------------------------------------------

void IOWrapper::readInputFormat(FILE * in)

{

    if (in == 0) in = infile;

    try {

      inputFormat = readStringInt(in);

      MSG(80, "inputFormat = %i\n", inputFormat );

    } catch( IOException e ) {

      e.addBeforeMessage("readInputFormat():: ");

      throw e;

    }

}



//------------------------------------------------------------------------

void IOWrapper::readInputDIM(FILE * in)

{

    if (in == 0) in = infile;

    inputDIM = readStringInt( in );

}





//----------------------------------------------------------------------------

char IOWrapper::readEntry(FILE * in)

{

  if (in == 0) in = infile;

  char ch;

  do {

    

    ch = fgetc( in );

    

    if ( ferror(in) )

		throw IOException("readEntry() failed");

    

  } while( !feof(in) && isspace(ch) );

  return ch;

}





//----------------------------------------------------------------------------

int IOWrapper::readRecordTag(char id[2], FILE * in)

{

    if (in == 0) in = infile;

    char ch = '#';

    

    //we skip spaces (including newline & return)

    do {

      ch = fgetc( in );

      if ( feof(in) ) return 1;

    } while( isspace(ch) );

    

    //read two characters

    id[0] = ch;

    id[1] = fgetc( in );

    return feof(in);

}





//------------------------------------------------------------------------

int IOWrapper::readStringInt(FILE * in)

{

    if (in == 0) in = infile;

    int res = 0;

    if ( 1 != fscanf(in, " %i", &res) )

        throw IOException("readStringInt() failed");

    return res;

}





//------------------------------------------------------------------------

float IOWrapper::readFloatAscii(FILE * in)

{

    if (in == 0) in = infile;

    float res = 0;

    if ( 1 != fscanf(in, " %f", &res) )

        throw IOException("readFloatAscii() failed");

    return res;

}





//------------------------------------------------------------------------

char IOWrapper::readCharAscii(FILE * in)

{

  if (in == 0) in = infile;

  char ch[2];

  if ( inputBinary ) {

    if ( 1 != fread( ch, 1, 1, in))

      throw IOException("readCharAscii() failed");

  } else {

    if ( 1 != fscanf( in, "%1s", ch ) )

      throw IOException("readCharAscii() failed");

  }

  return ch[0];

}





//------------------------------------------------------------------------

IOWrapper::uint8 IOWrapper::readUInt8(FILE * in)

{

  if (in == 0) in = infile;

  if ( inputBinary ) {

    char ch;

    if ( 1 != fread( &ch, 1, 1, in))

      throw IOException("readUInt8() failed");

    return (uint8) ch;

	}

  

  unsigned int res;

  if ( 1 != fscanf(in, " %u", &res ) )

    throw IOException("readUInt8() failed");

  return res;

}





//------------------------------------------------------------------------

IOWrapper::int8 IOWrapper::readInt8(FILE * in)

{

  if (in == 0) in = infile;

  if ( inputBinary ) {

    char ch;

    if ( 1 != fread( &ch, 1, 1, in))

      throw IOException("readInt8() failed");

    return (int8)ch;

	}

  

  int res;

  if ( 1 != fscanf(in, " %u", &res ) )

    throw IOException("readInt8() failed");

  return res;

}





//------------------------------------------------------------------------

IOWrapper::int16 IOWrapper::readInt16(FILE * in)

{

    if (in == 0) in = infile;

    if ( inputBinary ) {

        char ch[2];

        if ( 2 != fread( ch, 1, 2, in))

            throw IOException("readInt16() failed");

        if ( inputBinary == 2 ) swapBytes( ch, 2 );

        return *( (int16*) ch );

	}

    

    int res;

    if ( 1 != fscanf(in, " %i", &res ) )

        throw IOException("readInt16() failed");

    return res;

}





//------------------------------------------------------------------------

IOWrapper::uint16 IOWrapper::readUInt16(FILE * in)

{

  if (in == 0) in = infile;

  if ( inputBinary ) {

    char ch[2];

    if ( 2 != fread( ch, 1, 2, in))

      throw IOException("readUInt16() failed");

    if ( inputBinary == 2 ) swapBytes( ch, 2 );

    return *( (uint16*) ch );

	}

  

  int res;

  if ( 1 != fscanf(in, " %i", &res ) )

    throw IOException("readUInt16() failed");

  return res;

}



//------------------------------------------------------------------------

IOWrapper::uint32 IOWrapper::readUInt32(FILE * in)

{

    if (in == 0) in = infile;

    if ( inputBinary ) {

        char ch[4];

        if (4 != fread( ch, 1, 4, in))

            throw IOException("readUInt32() failed");

        if ( inputBinary == 2 ) swapBytes( ch, 4 );

        return *( (uint32*) ch );

    }

    

    unsigned long res = 0;

    if ( 1 != fscanf(in, " %lx", &res ) )

        throw IOException("readUInt32() failed");

    return Name( res );

}





//------------------------------------------------------------------------

float IOWrapper::readReal32(FILE * in)

{

    if (in == 0) in = infile;

    if ( inputBinary ) {

        char ch[4];

        if ( 4 != fread( ch, 1, 4, in) )

            throw IOException("readReal32() failed");

        if ( inputBinary == 2 ) swapBytes( ch, 4 );

        return *( (float*) ch );

	}

    

    float res;

    if ( 1 != fscanf(in, " %f", &res ) )

        throw IOException("readReal32() failed");

    return res;

}





//------------------------------------------------------------------------

void IOWrapper::readReal32(real a[], int n, FILE * in)

{

    if (in == 0) in = infile;

    for(int d = 0; d < n; ++d )

        a[ d ] = readReal32( in );

}



//------------------------------------------------------------------------

/**vectors in the file are made of (inputDIM) members

while in the program, they have DIM, hence this code,

to avoid addressing unallocated memory

*/



void IOWrapper::readReal32Vect(real a[], FILE * in)



{

    if (in == 0) in = infile;

    int d;

    if ( inputDIM <= DIM ) {

        

        for(d = 0; d < inputDIM; ++d )

            a[ d ] = readReal32( in );

        for(; d < DIM; ++d )

            a[ d ] = 0;

        

    } else {

        

        for(d = 0; d < DIM; ++d )

            a[ d ] = readReal32( in );

        for(; d < inputDIM; ++d )

            readReal32( in );

        

    }

}





//============================================================================

//=========                    OUTPUT                              ===========

//============================================================================





void IOWrapper::setOutputFile( FILE * f )

{

    closeOutputFile();

    outfile = f;

}



//------------------------------------------------------------------------

/**openOutputFile open the output file and labels it whether it is to be read only, 

and whether to write in binary. By default append is 0, 

so the pointer is set to the beginning of the file, 

otherwise the pointer is set at the end of the current file. 

*/



void IOWrapper::openOutputFile(const char * filename, const int append)

{

    char code[3] = { 'w', '\0', '\0' };

    code[0] = append ? 'a':'w';

    if ( writeBinary )

        code[1] = 'b';

    

    setOutputFile( fopen(filename, code) );

    

    checkOutputFile();

}



//------------------------------------------------------------------------

void IOWrapper::openOutputFile(const int append )

{

    openOutputFile( RESULT_OUT, append );

}



//------------------------------------------------------------------------

void IOWrapper::closeOutputFile()

{

    if ( outfile && ( outfile != stdout ))

        fclose( outfile );

    outfile = 0;

}



//------------------------------------------------------------------------

void IOWrapper::checkOutputFile() const 

{

    if ( outfile == 0 )

        throw IOException("IOWrapper:: output file == 0");

    if ( ferror( outfile ) )

        throw IOException("IOWrapper:: ferror( output file )");

}    





//------------------------------------------------------------------------

void IOWrapper::writeFileFormatVersion(const int version, FILE * out)

{

  if (out == 0) out = outfile;

  writeRecordTag("ht", out);

  writeIntAscii(version, out);

}



//------------------------------------------------------------------------

void IOWrapper::writeBinarySignature(FILE * out)

{

  if (out == 0) out = outfile;

  if ( writeBinary ) {

    writeRecordTag("bs", out);

    uint16 x = 1; 

    if (2 != fwrite( &x, 1, 2, out))

      throw IOException("fwrite() failed in writeBinarySignature");

  }

}



//------------------------------------------------------------------------

void IOWrapper::writeCharAscii(const char ch, FILE * out)

{

  if (out == 0) out = outfile;

  if ( writeBinary ) {

    

    if (1 != fwrite( &ch, 1, 1, out ))

      throw IOException("fwrite() failed in writeCharAscii");

    

  } else {

    

    if ( 2 > fprintf(out, " %c", ch ))

      throw IOException("fprintf() failed in writeCharAscii");

  }

}





//------------------------------------------------------------------------

void IOWrapper::writeUInt8(const int n, FILE * out)

{

    if (out == 0) out = outfile;

    if ( writeBinary ) {

      

      if (( n < 0 ) || ( n > 255 )) 

        MSG.error("writeUInt8", "value out of range");

        

      char ch = (uint8) n;

      if (1 != fwrite( &ch, 1, 1, out ))

        throw IOException("fwrite() failed in writeUInt8");



    } else {

    

      if ( 2 > fprintf(out, " %u", n ))

        throw IOException("fprintf() failed in writeUInt8");

      

    }

}



//------------------------------------------------------------------------

void IOWrapper::writeInt8(const int n, FILE * out)

{

  if (out == 0) out = outfile;

  if ( writeBinary ) {

    

    if (( n < -127 ) || ( n > 128 ))

      MSG.error("writeInt8", "value out of range");

    

    char ch = (int8) n;

    if (1 != fwrite( &ch, 1, 1, out ))

      throw IOException("fwrite() failed in writeInt8");

    

  } else {



    if ( 2 > fprintf(out, " %u", n ))

      throw IOException("fprintf() failed in writeInt8");

  }

}



//------------------------------------------------------------------------

void IOWrapper::writeInt16(const int n, FILE * out)

{

  if (out == 0) out = outfile;

  if ( writeBinary ) {

    

    if (( n < -32767 ) || ( n > 32768 ))

        MSG.error("writeInt16", "value out of range");

    

    char ch[2];

    *( (int16*) ch ) = (int16) n;

    if ( 2 != fwrite( ch, 1, 2, out ))

      throw IOException("fwrite() failed in writeInt16");



  } else {



    if ( 2 > fprintf(out, " %i", n ))

      throw IOException("fprintf() failed in writeInt16");



  }

}



//------------------------------------------------------------------------

void IOWrapper::writeInt16(const int n, const int m, FILE * out)

{

  writeInt16(n, out); 

  writeInt16(m, out); 

}



//------------------------------------------------------------------------

void IOWrapper::writeInt16(const int n, const int m, const int p, FILE * out)

{

  writeInt16(n, out); 

  writeInt16(m, out); 

  writeInt16(p, out); 

}



//------------------------------------------------------------------------

void IOWrapper::writeUInt16(const uint32 n, FILE * out)

{

  if (out == 0) out = outfile;

  if ( writeBinary ) {

    

    if ( n >= 65536 )

      MSG.error("writeUInt16", "value out of range");

    

    char ch[2];

    *( (uint16*) ch ) = (uint16) n;

    if ( 2 != fwrite( ch, 1, 2, out ))

      throw IOException("fwrite() failed in writeUInt16");

    

  } else {

    

    if ( 2 > fprintf(out, " %i", n ))

      throw IOException("fprintf() failed in writeUInt16");

    

  }

}



//------------------------------------------------------------------------

void  IOWrapper::writeRecordTag(const char id[2], FILE * out)

{

  if (out == 0) out = outfile;

  //the identifier is made of three characters: a new line, and two char given in id[]

  char ch[4] = { '\n', id[0], id[1], '\0' };

  if ( writeBinary )

    fwrite( ch, 1, 3, out);

  else

    fprintf(out, ch);

}





//------------------------------------------------------------------------

void IOWrapper::writeUInt32(const uint32 n, FILE * out)

{

  if (out == 0) out = outfile;

  if ( writeBinary ) {



    char ch[4];

    *( (uint32*) ch ) = n;

    if (4 != fwrite( ch, 1, 4, out ))

      throw IOException("fwrite() failed in writeUInt32");

  } else {

    if ( 2 > fprintf(out, " %x", n))

      throw IOException("fprintf() failed in writeUint32");

  }

}





//------------------------------------------------------------------------

void IOWrapper::writeReal32(const real x, FILE * out)

{

    if (out == 0) out = outfile;

    if ( writeBinary ) {

      char ch[4] = { 0, 0, 0, 0 };

      *( (float*) ch ) = x;

      if (4 != fwrite( ch, 1, 4, out ))

        throw IOException("fwrite() failed in writeReal32");

    } else {

      if( 2 > fprintf(out, " %.4f", x ))

        throw IOException("fprintf() failed in writeReal32");

    }

}



//------------------------------------------------------------------------

void IOWrapper::writeReal32Vect(const int n, const real * a, const bool formattingNewLineBefore, FILE * out)

{

  if (out == 0) out = outfile;

  if ( formattingNewLineBefore && ! writeBinary)

    fprintf(out, "\n");

  for(int d = 0; d < n; ++d )

    writeReal32( a[ d ], out );

}





//------------------------------------------------------------------------

void IOWrapper::writeFormattingNewLine(FILE * out)

{

    if (out == 0) out = outfile;

    if (!writeBinary)

        fprintf(out, "\n");

}



//------------------------------------------------------------------------

void IOWrapper::writeFormattingSpace(const int N, FILE * out)

{

  if (out == 0) out = outfile;

  if (!writeBinary) 

    switch( N ) {

      case 0: break;

      case 1: fprintf(out, " "); break;

      case 2: fprintf(out, "  "); break;

      case 3: fprintf(out, "   "); break;

      case 4: fprintf(out, "    "); break;

      case 5: fprintf(out, "     "); break;

      default:

        for(int ii = 0; ii<N; ++ii )

          fprintf(out, " ");

        break;

    }

}





//------------------------------------------------------------------------

void IOWrapper::writeString(const char * S, FILE * out)

{

    if (out == 0) out = outfile;

    if (1 > fprintf(out, "%s", S))

        throw IOException("writeString() failed");

}



//------------------------------------------------------------------------

void IOWrapper::writeIntAscii(const int a, FILE * out)

{

    if (out == 0) out = outfile;

    if (2 > fprintf(out, " %i", a))

        throw IOException("fprintf() failed in writeString");

}



//------------------------------------------------------------------------

void IOWrapper::writeIntAscii(const int a, const int b, FILE * out)

{

    if (out == 0) out = outfile;

    if (4 > fprintf(out, " %i %i", a, b))

        throw IOException("fprintf() failed in writeString");

}





//------------------------------------------------------------------------

void IOWrapper::writeRealAscii(const real a, FILE * out)

{

    if (out == 0) out = outfile;

    if (6 > fprintf(out, " %.4f", a))

        throw IOException("fprintf() failed in writeString");    

}



//------------------------------------------------------------------------

void IOWrapper::print(const char * fmt, ...)

{

  if( outfile == 0 )

    throw IOException("IO.print() failed because outfile==0");    

  va_list args;

  va_start(args, fmt);

  vfprintf(outfile, fmt, args);

  va_end(args);

}





//------------------------------------------------------------------------

void IOWrapper::flush( FILE * out )

{

    if (out == 0) out = outfile;

    fflush( out );      

}





//------------------------------------------------------------------------

int  IOWrapper::resolvePath(char * directoryPath, const int ssize, const char * fileName)

{

  char * dname, fullpath[PATH_MAX];



  if ( 0 == realpath( fileName, fullpath )) {

    MSG.warning("IOWrapper::resolvePath: realpath(%s) failed", fileName);

    return 1;

  }

  

  //printf("resolvePath fileName: %s, fullname: %s", fileName, fullpath);



  if ( 0 == strchr( fullpath, '/' )) {

    snprintf(directoryPath, ssize, "." );

    return NO_ERROR;

  }



  //get the directory part

  dname = dirname( fullpath );

  //printf(" dirname: %s\n", dname);



  //store this in directoryPath:

  if ( ssize <= snprintf(directoryPath, ssize, "%s", dname)) {

    MSG.warning("IOWrapper::resolvePath: ", "pathname too long for given string");

    return 2;

  }

  return NO_ERROR;

}



//------------------------------------------------------------------------

int  IOWrapper::fixPath(char * fileName, const int ssize, const char * directoryPath)

{

  //get the file part

  char filepart[PATH_MAX], * dname = basename( fileName );

  //we copy, because we are going to write to fileName:

  snprintf(filepart, sizeof(filepart), "%s", dname ); 



  //fix the trailling '/' if present

  if ((*filepart != '\0') && (filepart[ strlen(filepart) -1 ] == '/'))

    filepart[ strlen(filepart) -1 ] = '\0';

  

  //printf("fixPath fileName: %s, directory: %s, basename: %s", fileName, directoryPath, dname);

  //recompose the fileName:

  if ( ssize <= snprintf(fileName, ssize, "%s/%s", directoryPath, filepart)) {

    MSG.warning("IOWrapper::fixPath: ", "fullname too long for given string");

    return 2;

  }

  

  //printf("  new: %s\n", fileName);

  return NO_ERROR;

}

