/****************************************************************************\
SeSiMCMC. Looking - for - motifs by MCMC project. (c) A. Favorov 2001-2013
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
For general describtion of the classes declared in the header, see headers.txt 
$Id: Exception.hpp 1904 2013-07-10 01:45:50Z favorov $
\****************************************************************************/

#ifndef _EXCEPTION_HPP
#define _EXCEPTION_HPP

#include <iostream>

using namespace std;

struct DumbException
{
	const char * info;
	DumbException(const char * str=""):info(str){};
};

inline
ostream & operator<< (ostream & o, const DumbException & de)
{
	o<<de.info;
	return o;
}

struct AtgcException : public DumbException
{
	AtgcException(const char * str=""):DumbException(str){};
};

struct IOStreamException : public DumbException
{
	IOStreamException(const char * str=""):DumbException(str){};
};


#endif //_EXCEPTION_HPP
