#include <iostream>
#include <iomanip>

#include "../Sequences.hpp"
#include "../SymbolsCounter.hpp"
#include "../MarkovChainState.hpp"

extern "C"
{
	#include "../Random.h"
}

unsigned int fake_data_positions[20]={84,81,78,75,72,69,66,63,60,50,
													47,44,41,38,35,32,29,26,23,20};

//here, we omit 5(4) and complement 9(8).
													
void createFakeExperimentalBunch(SequencesPile &sp,double noise_prob=0);

void createFakeExperimentalBunch(SequencesPile &sp,double noise_prob)
//*sp is already created
{
	//20 sequnces, 105 random letters of ATGC (1234),
	// pattern ATGGCCACTT with random change in one letter with probability prob
	// is included in positions: 95,90,85......0 of corresponding pattern.
	unsigned int seq_no=20;
	unsigned int seq_length=120;
	unsigned int p_leng=10;
	unsigned short tpattern[10]={1,2,3,3,4,4,1,4,2,2};
	unsigned short cpattern[10]={1,1,3,2,3,3,4,4,1,2};
	unsigned short *pattern;
	cerr<<"Positions of motifs in test bunch are:\n";
	for (unsigned int i=0;i<seq_no;i++)
	{
		unsigned short * a=(unsigned short*)calloc(seq_length+1,sizeof(unsigned short));
		a[seq_length]=0;
		pattern=(i==8)?cpattern:tpattern; //9 has complemeted seq
		//we just alloc it, it we be deallocated by 
		//SequencesBunch::~SequencesBunch
		for (unsigned int j=0;j<seq_length;j++) a[j]=1+(int)floor(4*uni());
		//random seq is ready. Adding pattern
		if (i!=4) //omitting 5-th sequence
		{
			if (uni()<noise_prob) //with noise or not
			{
				unsigned int ch_pos=(int)floor(p_leng*uni());
				unsigned short ch_let=1+(int)floor(4*uni());
				unsigned short old_let=pattern[ch_let];
				pattern[ch_pos]=ch_let;
				for (unsigned int j=0;j<p_leng;j++) a[fake_data_positions[i]+j]=
					pattern[j];
				pattern[ch_pos]=old_let;
			}
			else
			{
				for (unsigned int j=0;j<p_leng;j++) a[fake_data_positions[i]+j]=
					pattern[j];
			}
			cerr<<fake_data_positions[i]<<" ";
		}
		else
		{
			for (unsigned int j=0;j<seq_length;j++) a[j]=1;
			cerr<<fake_data_positions[i]<<"** ";
		}
		string name="Fake|seq";
		char str_num[5];
		sprintf(str_num,"%02i",i);
		name.append(str_num);
		//the sequence is prepared.
		sp.add(a,seq_length,name);
	}
	cerr<<endl;
}
int main()
{
	
}
