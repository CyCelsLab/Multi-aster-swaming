//RCS: $Id: parameter_typed.cc,v 2.2 2005/03/04 09:56:44 nedelec Exp $

#include "parameter_typed.h"

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//all print codes should start with a space, except for 'char'
//this space ensures separation of multiple entries when writing.
template<> char * ParameterTyped<int>::typeCode()                    const { return "int"; }
template<> char * ParameterTyped<int>::scanCode()                    const { return "%i"; }
template<> char * ParameterTyped<int>::printCode()                   const { return " %-9i"; }
template<> bool   ParameterTyped<int>::shouldBeInteger()             const { return 1; }
template<> bool   ParameterTyped<int>::shouldBePositive()            const { return 0; }

template<> char * ParameterTyped<unsigned int>::typeCode()           const { return "uint"; }
template<> char * ParameterTyped<unsigned int>::scanCode()           const { return "%u"; }
template<> char * ParameterTyped<unsigned int>::printCode()          const { return " %-9u"; }
template<> bool   ParameterTyped<unsigned int>::shouldBeInteger()    const { return 1; }
template<> bool   ParameterTyped<unsigned int>::shouldBePositive()   const { return 1; }

template<> char * ParameterTyped<long>::typeCode()                   const { return "long"; }
template<> char * ParameterTyped<long>::scanCode()                   const { return "%li"; }
template<> char * ParameterTyped<long>::printCode()                  const { return " %-9li"; }
template<> bool   ParameterTyped<long>::shouldBeInteger()            const { return 1; }
template<> bool   ParameterTyped<long>::shouldBePositive()           const { return 0; }

template<> char * ParameterTyped<unsigned long>::typeCode()          const { return "ulong"; }
template<> char * ParameterTyped<unsigned long>::scanCode()          const { return "%lx"; }
template<> char * ParameterTyped<unsigned long>::printCode()         const { return " 0x%-9lx"; }
template<> bool   ParameterTyped<unsigned long>::shouldBeInteger()   const { return 1; }
template<> bool   ParameterTyped<unsigned long>::shouldBePositive()  const { return 1; }

template<> char * ParameterTyped<float>::typeCode()                  const { return "float"; }
template<> char * ParameterTyped<float>::scanCode()                  const { return "%f"; }
template<> char * ParameterTyped<float>::printCode()                 const { return " %-9.4g"; }
template<> bool   ParameterTyped<float>::shouldBeInteger()           const { return 0; }
template<> bool   ParameterTyped<float>::shouldBePositive()          const { return 0; }

template<> char * ParameterTyped<double>::typeCode()                 const { return "double"; }
template<> char * ParameterTyped<double>::scanCode()                 const { return "%lf"; }
template<> char * ParameterTyped<double>::printCode()                const { return " %-9.4lg"; }
template<> bool   ParameterTyped<double>::shouldBeInteger()          const { return 0; }
template<> bool   ParameterTyped<double>::shouldBePositive()         const { return 0; }

template<> char * ParameterTyped<char>::typeCode()                   const { return "char"; }
template<> char * ParameterTyped<char>::scanCode()                   const { return "%c"; }
template<> char * ParameterTyped<char>::printCode()                  const { return "%c"; }
template<> bool   ParameterTyped<char>::shouldBeInteger()            const { return 0; }
template<> bool   ParameterTyped<char>::shouldBePositive()           const { return 0; }
