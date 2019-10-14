//$Id$
#ifndef _CONFREAD_H_
#define _CONFREAD_H_

#include <stdio.h>

#define Conf_File_Line_Len 512
#define Int_String_Len 10
#define Double_String_Len 30
#define Long_String_len 20
#define Boolean_String_Len 5
#define List_String_Len 450
#define Max_Int_List_Len 300

int ReadConfigString(FILE *ifile, const char *tag, 
		char *value, unsigned int allocated_value_len,
		unsigned int is_obligatory);
int ReadConfigBoolean(FILE *ifile,const char *tag,unsigned int *value,
		unsigned int is_obligatory);
int ReadConfigInt(FILE *ifile, const char *tag, int *value,
		unsigned int is_obligatory);
int ReadConfigLong(FILE *ifile, const char *tag, long *value,
		unsigned int is_obligatory);
int ReadConfigUnsignedInt(FILE *ifile, const char *tag, unsigned int *value,
		unsigned int is_obligatory);
int ReadConfigUnsignedLong(FILE *ifile, const char *tag, unsigned long *value,
		unsigned int is_obligatory);
int ReadConfigDouble(FILE *ifile, const char *tag, double *value,
		unsigned int is_obligatory);
int ReadConfigIntList(FILE *ifile, const char *tag, int **list, int* list_lenght,
		unsigned int is_obligatory);
//it allocates the list and returns its length in list_lenght. 

//1 - OK
//0 - do not find. (mute)
//-1 - the tag is but it is not of the desired type (mute)
//-10 - file format error
//-100 - file open eror
#endif // _CONFREAD_H_
   
