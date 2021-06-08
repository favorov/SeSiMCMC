#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>

using namespace std;

#include "../Exception.hpp"
#include "../MarkovChainState.hpp"
#include "../Sequences.hpp"
#include "../MCMC.hpp"
#include "../ResultsSet.hpp"


void howtouse()
{
	cout<<"Usage: ResultsSetTest FastAfilename <MarkovChainState\n"
			<<"The param is oblitatory.\n";
	
};


int main(int argc, char ** argv)
{
	SequencesPile sp;
	string FastAName;
	AcgtLetterOrder order;
	Diagnostics diagnostics;
	if (argc==2)
	{
		FastAName=argv[1];
	}
	else
	{
		howtouse();
		return 10;
	}

	ifstream fastainput(FastAName.c_str());
	if (! fastainput.good()) 
	{
		fprintf(stderr,"Input file %s cannot be opened for read.\n",FastAName.c_str());
		return 10;
	}

	try{
		fastainput>>sp;
	} catch(DumbException de) {cout<<de;return 5;}
	
	MarkovChainState state(sp.size());

	try{
	cin>>state;
	} catch(DumbException de) {cerr<<de;return 1;}
/*
	(
		const SequencesPile & Sequences, 
		double motif_absence_prior,
		unsigned long local_step_cycles_between_adjustments,
		unsigned short if_common_background,
		const vector<double> & background_opt,
		unsigned short be_quiet,
		ostream & logstream=cerr
	) throw (DumbException);
 */ 
	LookingForAtgcMotifsInTwoThreadsMultinomialGibbs sampler
	(
		sp,
		0.01,
		0,
		off,
		off,
		off,
		0,
		0,
		*new vector<double>(4),
		1);

/*
	Profile (const SequencesPile & sp, const MarkovChainState & mcs,
		const LookingForAtgcMotifsMultinomialGibbs & sampler,
		unsigned short was_length_adjusted_p,
		unsigned int minimal_motif_length,
		unsigned int maximal_motif_length,
		unsigned short was_mode_slow_p,
		unsigned int retrieve_other_sites_threshold, //o if not ot retrieve
		unsigned int flanks_len=5,
		unsigned short show_lettrs_in_gap=0,
		unsigned int all_chains_failed_n=0,
		const set <pair<unsigned int,double> > & G_values_p=
			*(new const set <pair<unsigned int,double> >),
		unsigned int output_complements_as_in_seq=1);
 */ 

	sampler.the_state=state;
	Profile profile
			( sp, state, 
				sampler,
				1,
				5,
				45,
				0,
				0,
				1,
				0,
				diagnostics,
				order
				);
	output_html_header(cout,"Test for ResultsSetTest");
	profile.html_output(cout);
	profile.short_html_output(cout);
	output_html_footer(cout);
	return 0;
}
