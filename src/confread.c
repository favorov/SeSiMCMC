//$Id: confread.c 1014 2009-03-01 16:50:36Z favorov $
//
#include <stdlib.h>
#include <string.h>
#include <stdio.h>  
#include "confread.h"

//1 - OK
//0 - do not find. (mute)
//-5 - the tag is but it is not of the desired type
//-10 - file format error
//-100 - file open eror

int ReadConfigString(FILE *ifile, const char *tag, 
		char *value, unsigned int allocated_value_len,unsigned int is_obligatory)
{
	char buf[Conf_File_Line_Len+1];
	char string[Conf_File_Line_Len];
	char *src,*token;
	if (!ifile)
	{
		fprintf(stderr,"NULL file to read %s from\n",tag);
		return -100;
	}
  fseek(ifile, 0, SEEK_SET);
	while (fgets(buf,Conf_File_Line_Len,ifile))//Main cycle 
	{
		src=buf+strlen(buf)-1;
		*src=(*src==13)?0:*src; //to fgets()'s trailing endline
		*src=(*src==10)?0:*src;
		src--;										
		*src=(*src==13)?0:*src;                  	
		*src=(*src==10)?0:*src; //removing any kind of endline
		buf[Conf_File_Line_Len-1]='\0';   
		if(*buf=='\0') continue;//empty line

		src = strchr(buf,'#');
		if (src) *src='\0';  //we has removed comment

		src = buf+strspn(buf," \t"); //remove starting spaces
		if (*src=='\0') continue;//empty line
		token = strchr(src,'=');
		if (!token)
		{
			fprintf(stderr,"Cannot parse a line in config file for it has no \'=\':\n%s\n.\n",buf);
			return -10;
		}
		*(token++)='\0';
		strcpy(string,token+strspn(token," \t")); //we've processed symbols after = 
		token=strtok(src," \t");
		if (!token)
		{
			fprintf(stderr,"Empty tag line in config file\n");
			return -10;
		}
		if (strcmp(token,tag)) continue; //we do need this line
		strncpy(value,string,allocated_value_len-1);
		value[allocated_value_len-1]='\0';
		return 1;
	}
	if (is_obligatory) fprintf (stderr,"Cannot find tag %s in config file.\n",tag);
	return 0; //ne nashli
}

int ReadConfigBoolean(FILE *ifile, const char *tag,unsigned int *value,unsigned int is_obligatory)
{
	char *str_value;
	int read_result;
	str_value=(char*)calloc(Boolean_String_Len+1,sizeof(char));
	read_result=ReadConfigString(ifile,tag,str_value,
																Boolean_String_Len+1,is_obligatory); 
	if (read_result<=0) return read_result;
	*value=2; //nonsense nalue
	if (!strncmp(str_value,"true",Boolean_String_Len+1))*value=1;	
	if (!strncmp(str_value,"yes",Boolean_String_Len+1))*value=1;	
	if (!strncmp(str_value,"on",Boolean_String_Len+1))*value=1;	
	if (!strncmp(str_value,"1",Boolean_String_Len+1))*value=1;	
	if (!strncmp(str_value,"false",Boolean_String_Len+1))*value=0;	
	if (!strncmp(str_value,"no",Boolean_String_Len+1))*value=0;	
	if (!strncmp(str_value,"off",Boolean_String_Len+1))*value=0;	
	if (!strncmp(str_value,"0",Boolean_String_Len+1))*value=0;	
	if (*value==2) 
	{
		fprintf(stderr,"Tag %s is not Boolean.\n",tag);
		return -5;
	}
	free(str_value);
	return 1;
}

int ReadConfigInt(FILE *ifile, const char *tag, int *value,unsigned int is_obligatory)
{
	char *str_value;
	int read_result;
	char * ptr;
	str_value=(char*)calloc(Int_String_Len+1,sizeof(char));
	read_result=ReadConfigString(ifile,tag,str_value,
																Int_String_Len+1,is_obligatory); 
	if (read_result<=0) return read_result;
	*value=strtol(str_value,& ptr,0);
	if (ptr!=(str_value+strlen(str_value)) ) 
	{
		fprintf(stderr,"Tag %s is not integer.\n",tag);
		return -5;
	}
	free(str_value);
	return 1;
}

int ReadConfigLong(FILE *ifile, const char *tag, long *value,unsigned int is_obligatory)
{
	char *str_value;
	int read_result;
	char * ptr;
	str_value=(char*)calloc(Int_String_Len+1,sizeof(char));
	read_result=ReadConfigString(ifile,tag,str_value,
																Int_String_Len+1,is_obligatory); 
	if (read_result<=0) return read_result;
	*value=strtol(str_value,& ptr,0);
	if (ptr!=(str_value+strlen(str_value)) ) 
	{
		fprintf(stderr,"Tag %s is not integer.\n",tag);
		return -5;
	}
	free(str_value);
	return 1;
}

int ReadConfigUnsignedInt(FILE *ifile, const char *tag, unsigned int *value,
		unsigned int is_obligatory)
{
	char *str_value;
	int read_result;
	char * ptr;
	str_value=(char*)calloc(Int_String_Len+1,sizeof(char));
	read_result=ReadConfigString(ifile,tag,str_value,
																Int_String_Len+1,is_obligatory); 
	if (read_result<=0) return read_result;
	*value=strtoul(str_value,& ptr,0);
	if (ptr!=(str_value+strlen(str_value)) ) 
	{
		fprintf(stderr,"Tag %s is not integer.\n",tag);
		return -5;
	}
	free(str_value);
	return 1;
}

int ReadConfigUnsignedLong(FILE *ifile, const char *tag, unsigned long *value,
		unsigned int is_obligatory)
{
	char *str_value;
	int read_result;
	char * ptr;
	str_value=(char*)calloc(Int_String_Len+1,sizeof(char));
	read_result=ReadConfigString(ifile,tag,str_value,
																Int_String_Len+1,is_obligatory); 
	if (read_result<=0) return read_result;
	*value=strtoul(str_value,& ptr,0);
	if (ptr!=(str_value+strlen(str_value)) ) 
	{
		fprintf(stderr,"Tag %s is not integer.\n",tag);
		return -5;
	}
	free(str_value);
	return 1;
}

int ReadConfigDouble(FILE *ifile, const char *tag, double *value,unsigned int is_obligatory)
{
	char *str_value;
	int read_result;
	char * ptr;
	str_value=(char*)calloc(Int_String_Len+1,sizeof(char));
	read_result=ReadConfigString(ifile,tag,str_value,
																Int_String_Len+1,is_obligatory); 
	if (read_result<=0) return read_result;
	*value=strtod(str_value,& ptr);
	if (ptr!=(str_value+strlen(str_value)) ) 
	{
		fprintf(stderr,"Tag %s is not double.\n",tag);
		return -5;
	}
	free(str_value);
	return 1;
}

int ReadConfigIntList(FILE *ifile, const char *tag, int **list, int* list_lenght,
		unsigned int is_obligatory)
//it allocates the list and returns its length in list_lenght.
{
	char *str_value,*token;
	int read_result;
	char * ptr;
	int * list_value;
	int read=0;

	str_value=(char*)calloc(List_String_Len+1,sizeof(char));
	read_result=ReadConfigString(ifile,tag,str_value,
																List_String_Len+1,is_obligatory);
	list_value=(int*)calloc(Max_Int_List_Len,sizeof(int));
	if (read_result<=0) return read_result;
	while((token=strtok(read?NULL:str_value,";:,")))
	{
		list_value[read]=strtol(token,& ptr,0);
		if (ptr!=(token+strlen(token)) ) 
		{
			fprintf(stderr,"Tag %s is not an integer list.\n",tag);
			return -5;
		}
		read++;
	}
	*list=(int*)calloc(read,sizeof(int));
	*list_lenght=read;
	for(read=0;read<*list_lenght;read++)
		(*list)[read]=list_value[read];
	free(str_value);
	free(list_value);
	return 1;
}

