/****************************************************************************\
SeSiMCMC. Looking - for - motifs by MCMC project. (c) A. Favorov 2001-2019
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
For general describtion of the classes declared in the header, see headers.txt 
$Id$
\****************************************************************************/

#ifndef _LOGGER_HPP
#define _LOGGER_HPP

#include <time.h>
#include <iomanip>
#include <string>
#include <fstream>
#include <sstream>
//the class has only constructor and destructor
//it tries to append a start line to the logfile on constructing and
//a finish line on destructing.
//

using namespace std;

#define WAIT_FOR_OPEN 3 //seconds
class LogRecorder
{
private:
	LogRecorder();
	LogRecorder(const LogRecorder &);
	LogRecorder & operator=(const LogRecorder &);
	string logfilename,Id,IP;
	time_t starttime;
public:
	LogRecorder(const string & logfilename_par, 
							const string & Id_par="undefined",
							const string & IP_par="unknown");
	double how_long() const {return difftime(time(NULL),starttime);};
	~LogRecorder();
};

inline 
LogRecorder::LogRecorder
	(
		const string & logfilename_par,
		const string & Id_par, 
		const string & IP_par
	):logfilename(logfilename_par),Id(Id_par),IP(IP_par)
{
	starttime=time(NULL);
	if (!(logfilename=="/dev/null"))
	{
		ostringstream output;
		ofstream logfile;
		char tim[80];tim[79]=0;
		strftime(tim,79,"%Y %b %d %X %Z",localtime(&starttime));
		output<<"+Id:"<<Id<<" started from IP:"<<IP<<" at "<<tim<<"."<<endl;
		time_t time_ch=time(NULL);
		do logfile.open(logfilename.c_str(),ios::app);
		while ( !logfile && difftime(time(NULL),time_ch)<WAIT_FOR_OPEN);
		if (logfile) 
		{
			logfile<<output.str();
			logfile.close();
		}
	}
}

inline 
LogRecorder::~LogRecorder()
{
	ostringstream output;
	ofstream logfile;
	time_t finishtime;
	finishtime=time(NULL);
	char tim[80];tim[79]=0;
	strftime(tim,79,"%Y %b %d %X %Z",localtime(&finishtime));
	output<<"-Id:"<<Id<<" (IP:"<<IP<<") finished at "<<tim<<"."<<endl
			<<"%Id:"<<Id<<" worked "<<difftime(finishtime,starttime)<<" seconds."<<endl;
	time_t time_ch=time(NULL);
	do logfile.open(logfilename.c_str(),ios::app);
	while ( !logfile && difftime(time(NULL),time_ch)<WAIT_FOR_OPEN);
	if (logfile) 
	{
		logfile<<output.str();
		logfile.close();
	}
}
#endif //_LOGGER_HPP

