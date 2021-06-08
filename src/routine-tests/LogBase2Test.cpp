#include <iomanip>
#include "SymbolsCounter.hpp"


int main()
{
	double q=1.;
	for (int i=1;i<=40;i++) q*=.5;
	for (;q<=4.001;q*=2)
		cout<<"log_2("<<q<<")="<<log_2(q)<<endl;
	return 0;
}
