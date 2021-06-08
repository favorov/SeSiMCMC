#include <vector>
#include <string>

#include "../Atgc.hpp"

int main()
{
	vector<unsigned short> atgc,atgc1;
	string str;
	str="attttaaaagcccc";
	Atgc::string2atgc(str,atgc);
	cout<<"Original string:\n"<<str<<endl
			<<"Back converted:"<<Atgc::atgc2string(atgc,str)<<endl;
	Atgc::complement(atgc,atgc1);
	cout<<"Copy complement:"<<Atgc::atgc2string(atgc1,str)<<endl;
	Atgc::complement(atgc);
	cout<<"Inplace complement:"<<Atgc::atgc2string(atgc,str)<<endl;
	Atgc::complement(atgc,atgc);
	cout<<"One more inplace :"<<Atgc::atgc2string(atgc,str)<<endl;
}
