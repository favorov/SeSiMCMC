/****************************************************************************\
SeSiMCMC. Looking - for - motifs by MCMC project. (c) A. Favorov 2001-2021
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
For general describtion of the classes declared in the header, see headers.txt
$Id$
\****************************************************************************/

#ifndef _SYMBOLS_COUNTER_HPP
#define _SYMBOLS_COUNTER_HPP

#include <math.h>

#include <vector>
#include <iostream>
#include <iomanip>

using namespace std;

#include "Exception.hpp"
#include "Sequences.hpp"
#include "MarkovChainState.hpp"

#define EPSILON 1E-12

const long double ln2=log(2.);

inline double log_2(double x)
{
	return (log(x)/ln2);
}

/*
As usual, the letters are 1..letters_in_alphabet
the positions in pattern and the sequence number are 0-based
here, we maintain:

We use SymbolsCounter::letters or complete dataset letters count.
It is accesible via complete_count(letter). For AtgcSymbolsCounter,
we use private *r vector of the same structure. We need it because
every sequence can be treated as complement.

its relative is double * psc - a vector of pseudocounts, i.e
b_i=B*pho_i , where pho_i is r_i/sum(r_i) and B=sqrt(sequence_number).
The vector is accessible via pseudocount(letter).

Their sum B is to be accessed via pseudocounts_sum().

Every call of calculate() updates the set of members described below.

c  maintains positional information about counts of letters by position in
the pattern. Th counts are collected from all patterns except current
it can be accessed as pattern_count(m,n) for occurences of letter n in
position m. The complement pattern are cpunted backwards.

s nonpositional information about all letters occurences in nonpattern
positions the information is gathered from all sequences except current
It is accessed as nonpattern(n) for letter n.

There is a special call calculate_complete(vector<unsigned int> positions)
that counts c and s form complete dataset, without omitting current sequence.

For a spaced motif, we still support full length, but the gap positions have
background probabilities, so almost all the procedurec are unchanged.
*/

template <class T>
inline void print_foreground_probablities(const T & inst, ostream & o )
{
	o<<"PL="<<inst.pattern_length<<endl;
	for (unsigned short n=1;n<=inst.letters_in_alphabet;n++)
	{
		for (unsigned int m=0;m<inst.pattern_length;m++)
			o<<setprecision(2)<<inst.foreground_probability(m,n)<<" ";
		o<<endl;
	};
	o<<endl<<"++++++++++++++++"<<endl<<flush;
};


class SymbolsCounter
{
//we allocate all counters with letters_in_alphabet+1
//to use 1-based letters and to reserve 0 for control
//purposes.
	SymbolsCounter(const SymbolsCounter &);
	SymbolsCounter & operator= (const SymbolsCounter &);
friend void print_foreground_probablities<SymbolsCounter>(const SymbolsCounter & inst,ostream & o);
friend class KullbakCounter;
protected:
	const SequencesPile & Sequences;
	const unsigned int letters_in_alphabet;
	unsigned int pattern_length;
	unsigned int spacer_5_end;
	unsigned int spacer_3_end;
	unsigned int is_spaced;
	unsigned int informative_part_length;

	unsigned short if_common_background;

	//the three fields below make a kind of statistics of letters in the bunch

	unsigned long total_letters;
	//it is amount af all letters in all slots

	unsigned long current_total_letters;
	//it is amount af all letters in all slots except excluded sequence

	unsigned long motifs_count;
	//it is a count for all motifs counted in pattern counters

	vector <unsigned int> letters;
	//counts for letters 1..letters_in_alphabet in complete data set,
	//count for n-th letter is r[n], r[0] is for masked positions.
	//it is not the same thing with Sequences::letters for here we
	//take into account the fact that some sequences are read as complement.

	vector < vector<unsigned int> > letters_per_sequence;

	// residues_per_sequence[1][5] is count af A's in sequence 5(0-based)
	// residues are 1 - based, so we extract 1 from address.
	// format is letters_per_sequence[letter][sequence]

	//the strange letters for definitions below are taken from Lawrence93 paper.

	vector <double> b;
	//pseudocounts for letters 1..letters_in_alphabet in complete data set,
	//pseudocount for n-th letter is b[n] (b[0]) is to match counters,
	//which all have 0 element allocated, but not used (for masked)
	//we will recount them every time after r[] changes

	double B;
	//sum of all pseudocounts b[]

	vector <unsigned int> c;
	//counts of pattern positions, n-th letter in m-th (0-based) position is
	//c[m*letters_in_alphabet+n]
	//position is 0-based; letter is 1-based
	//(0 - th letter is x, it is allocated to simplify the program).
	//see headers.txt about the pattern position representation

	//we do not want to organise it in [][] way because it will be too
	//time - hungry to allocate/deallocate such a thing every count.

	vector <unsigned int> c5prev;

	//it is a vector [letters+1] which carry the first motif position weight
	//vector before the last 3'-shift (which is shifted away from c).
	//we need it to recalculate the entropy incrementally.
	//It makes sense only for calculate_incrementally_after_3_end_shift.
	vector <unsigned int> c5prev_second;
	//the same, but for the second half of a spaced site

	// old position  ..***..***...
	// new (current) ...***..***..
	//               0123456789
	// c5prev is count of [2], c5prev_second - of [7]
	//
	vector <unsigned int> s;
	//counts for letters 1..letters_in_alphabet in non-pattern part of set,
	//count for n-th letter is r[n]

	vector <unsigned int> s5prev;
	//the same thing as c5prev, but for nonpattern counters.
	//

	unsigned int complete_count(unsigned int n) const
	//return complete dataset count for letter n
	{
		return letters[n];
	}
public:
	double pseudocount(unsigned int n) const
	//return pseudocount for letter n; n is 1 - based!
	{
		return b[n];
	}
//protected:
	unsigned int pattern_count(unsigned int m,unsigned int n) const
	//counts of pattern positions, n-th letter in m-th position is
	//c[m*letters_in_alphabet+n]
	//position is 0-based; letter is 1-based.
	{
		return c[m*letters_in_alphabet+n];
	}

	unsigned int nonpattern_count(unsigned int n) const
	//counts for letters 1..letters_in_alphabet in full dataset,
	//count for n-th letter is r[n-1]
	{
		return s[n];
	}
public:

	double pseudocounts_sum() const
	//The pseudocounts sum is accesed via pseudocounts_sum()
	{
		return B;
	}

	double change_pseudocounts_sum (double sum);

	double restore_pseudocounts_sum ();  //restores it to sqrt((double)Sequences.size())

	struct OtherLengthPWMException : public DumbException
	{
		OtherLengthPWMException(const char * str=""):DumbException(str){};
	};

	class PWM
	{
	private:
		vector<double> the_PWM;
	public:
		unsigned int letters_in_alphabet;
		unsigned int pattern_length;
		PWM(const SymbolsCounter & sc);
		PWM & operator = (const SymbolsCounter & sc);
		double foreground_probability(unsigned int m,unsigned int n) const;
	};
	//PWM
	SymbolsCounter
	(
		const SequencesPile & sp,
		unsigned int letters_in_alphabet_par,
		unsigned short if_common_background=0,
		const vector<double> & background_opt =*new vector<double>
	);

	//standard ?
	friend class SymbolsCounter::PWM;

	unsigned int motifs() const
	{
		return motifs_count;
	}

	unsigned long total() const
	{
		return current_total_letters;
	}

	unsigned long all_total() const
	{
		return total_letters;
	}

	virtual void calculate (const MarkovChainState & theState);
	//gets the pattern_count and nonpattern_count values from the full dataset
	//to get the pattern_count and nonpattern_count for the dataset
	//excluding current sequence use exclude_sequence.
	//Backwards, use include_sequence.
	//
	//
	virtual void exclude_sequence
			(
			 	const MarkovChainState & theState,
				unsigned int sequence
			);
	//All we do is removal of sequence from the sum.
	//Works only if the pattern_length is the same for
	//the state and the SymbolsCounter.
	//We remove the sequence gain from the counter
	//as it was done in accordance with theState.
	//If it is untruth, the result can be stange
	//
	virtual void include_sequence
			(
			 	const MarkovChainState & theState,
				unsigned int sequence
			);
	//All we do is addition of current_sequence from the sum in
	//accordance with theState.
	//Works only if the pattern_length is the same for
	//the state and the SymbolsCounter.
	//

	virtual void calculate_incrementally_after_3_end_shift
								(const MarkovChainState & theState);
	//an analog of calculate, but supposes that the last action was
	//calculate for preceding theState (or any equvalent, like
	//calculate+exclude_sequence+include_sequence)
	//


	virtual double foreground_probability(unsigned int m,unsigned int n)
				const;
	//n-th letter in m-th motif position probability
	//position is 0-based; letter is 1-based.

	virtual double foreground_weight(unsigned int m,unsigned int n)
				const;
	//n-th letter in m-th motif position weight
	//position is 0-based; letter is 1-based.

	virtual double background_probability(unsigned int n) const;
	//n-th letter background probability
	//position is 0-based; letter is 1-based.

	virtual double long NegativeEntropy ();
	//F in Lawrence93
	//
	virtual double long NegativeEntropy_change_after_3_end_shift() ;

	virtual void print(ostream & o) const;

	bool operator== (const SymbolsCounter & sc1) const;

	virtual ~SymbolsCounter(void){};
};

inline
ostream & operator << (ostream & o,const SymbolsCounter & sc)
{
	sc.print(o);
	return o;
}

inline
void SymbolsCounter::print(ostream & o) const
{
	o<<"Total counters:"<<endl;
	for(unsigned int l=1;l<=letters_in_alphabet;l++)
		o<<letters[l]<<((l==letters_in_alphabet)?"=":"+");
	o<<total_letters<<endl;
	o<<"Nonpattern counters:"<<endl;
	for(unsigned int l=1;l<=letters_in_alphabet;l++)
		o<<nonpattern_count(l)<<"  ";
	o<<endl<<"Pattern counters gathered from "<<motifs()<<" motifs:";
	for(unsigned int pos=0;pos<pattern_length;pos++)
	{
		o<<endl;
		for(unsigned int l=1;l<=letters_in_alphabet;l++)
			o<<pattern_count(pos,l)<<" ";
	}
	o<<endl;
	return;
}

inline
bool SymbolsCounter::operator== (const SymbolsCounter & sc1) const
{
	if (sc1.letters_in_alphabet!=letters_in_alphabet) return 0;
	if (sc1.pattern_length!=pattern_length) return 0;
	if (sc1.motifs_count!=motifs_count) return 0;
	if (sc1.total_letters!=total_letters) return 0;
	if (sc1.letters!=letters) return 0;
	if (sc1.letters_per_sequence!=letters_per_sequence) return 0;
	if (sc1.b!=b) return 0;
	if (sc1.B!=B) return 0;
	if (sc1.c!=c) return 0;
	if (sc1.s!=s) return 0;
	return 1;
}


class SymmetricSymbolsCounter:public SymbolsCounter
{
private:
	SymmetricSymbolsCounter(const SymmetricSymbolsCounter &);
	SymmetricSymbolsCounter & operator= (const SymmetricSymbolsCounter &);
	long double NegativeEntropyValue;
public:
	SymmetricSymbolsCounter
	(
		const SequencesPile & sp,
		unsigned int letters_in_alphabet_par,
		unsigned short if_common_background,
		const vector<double> & background_opt
	);
	virtual double foreground_probability(unsigned int m,unsigned int n)
				const;
	//n-th letter in m-th motif position probability
	//position is 0-based; letter is 1-based.

	virtual double foreground_weight(unsigned int m,unsigned int n)
				const;
	//n-th letter in m-th motif position weight
	//position is 0-based; letter is 1-based.
	virtual double long NegativeEntropy () ;
	//F in Lawrence93
	//
	virtual double long NegativeEntropy_change_after_3_end_shift() ;
	//we do not use any incremantal form here, but we
	//must provide a unified interface for all counters.
	bool operator== (const SymmetricSymbolsCounter & sc1) const;
};

inline
bool SymmetricSymbolsCounter::operator== (const SymmetricSymbolsCounter & sc1) const
{
	if (sc1.letters_in_alphabet!=letters_in_alphabet) return 0;
	if (sc1.pattern_length!=pattern_length) return 0;
	if (sc1.motifs_count!=motifs_count) return 0;
	if (sc1.total_letters!=total_letters) return 0;
	if (sc1.letters!=letters) return 0;
	if (sc1.letters_per_sequence!=letters_per_sequence) return 0;
	if (sc1.b!=b) return 0;
	if (sc1.B!=B) return 0;
	if (sc1.c!=c) return 0;
	if (sc1.s!=s) return 0;
	return 1;
}

class DoubletSymbolsCounter:public SymbolsCounter
{
private:
	long double NegativeEntropyValue;
public:
	DoubletSymbolsCounter
	(
		const SequencesPile & sp,
		unsigned int letters_in_alphabet_par,
		unsigned short if_common_background,
		const vector<double> & background_opt
	);

	virtual double foreground_probability(unsigned int m,unsigned int n)
				const;
	//n-th letter in m-th motif position probability
	//position is 0-based; letter is 1-based.

	virtual double foreground_weight(unsigned int m,unsigned int n)
				const;
	//n-th letter in m-th motif position weight
	//position is 0-based; letter is 1-based.
	virtual double long NegativeEntropy () ;
	//F in Lawrence93
	//
	virtual double long NegativeEntropy_change_after_3_end_shift() ;
	//we do not use any incremantal form here, but we
	//must provide a unified interface for all counters.

	bool operator== (const DoubletSymbolsCounter & sc1) const;
};

inline
bool DoubletSymbolsCounter::operator== (const DoubletSymbolsCounter & sc1) const
{
	if (sc1.letters_in_alphabet!=letters_in_alphabet) return 0;
	if (sc1.pattern_length!=pattern_length) return 0;
	if (sc1.motifs_count!=motifs_count) return 0;
	if (sc1.total_letters!=total_letters) return 0;
	if (sc1.letters!=letters) return 0;
	if (sc1.letters_per_sequence!=letters_per_sequence) return 0;
	if (sc1.b!=b) return 0;
	if (sc1.B!=B) return 0;
	if (sc1.c!=c) return 0;
	if (sc1.s!=s) return 0;
	return 1;
};

class AtgcSymbolsCounter:public SymbolsCounter
{
private:

	AtgcSymbolsCounter(const AtgcSymbolsCounter &);
	AtgcSymbolsCounter & operator= (const AtgcSymbolsCounter &);
	vector<unsigned int> counted_as_complement;

public:
	AtgcSymbolsCounter
	(
		const SequencesPile & sp,
		const MarkovChainState & mcs,
		unsigned short if_common_background_par,
		const vector<double> & background_opt
	);

	void current_sequence_was_complemented(unsigned int seq_no);
	//
	//one more iterative thing:
	//if a local step has changed a sequence
	//to comlement, it is enough to recount the
	//common counters (they are non pattern, no nonpattern)
	//and the counters devoted to the current (i.e. complemented)
	//sequence
	//

//All the four operations below  are rewritten
//to take into account the fact that the
//sequence can be read is complement,
//All of them suppose that r and psc
//counerts are already updated when we want to register
//that the a sequence is started to be treat as complemented
//by call of current_sequence_was_complemented after
//the MarkovChainState change.
	virtual void calculate (const MarkovChainState & theState);
	virtual void exclude_sequence
			(
			 	const MarkovChainState & theState,
				unsigned int sequence
			);
	virtual void include_sequence
			(
			 	const MarkovChainState & theState,
				unsigned int sequence
			);
	virtual void calculate_incrementally_after_3_end_shift
								(const MarkovChainState & theState);
	bool operator== (const AtgcSymbolsCounter & sc1) const;

	virtual void print(ostream & o) const;
};

inline
void AtgcSymbolsCounter::print(ostream & o) const
{
	o<<"Total counters:"<<endl;
	for(unsigned int l=1;l<=letters_in_alphabet;l++)
		o<<letters[l]<<((l==letters_in_alphabet)?"=":"+");
	o<<total_letters<<endl;
	o<<"Nonpattern counters:"<<endl;
	for(unsigned int l=1;l<=letters_in_alphabet;l++)
		o<<nonpattern_count(l)<<"  ";
	o<<endl<<"Pattern counters gathered from "<<motifs()<<" motifs:";
	for(unsigned int pos=0;pos<pattern_length;pos++)
	{
		o<<endl;
		for(unsigned int l=1;l<=letters_in_alphabet;l++)
			o<<pattern_count(pos,l)<<" ";
	}
	o<<endl;
	o<<"Direction of sequences accounted: "<<endl;
	for(unsigned int p=0;p<Sequences.size();p++)
		o<<counted_as_complement[p]<<" ";
	o<<endl;
	return;
}

inline
bool AtgcSymbolsCounter::operator== (const AtgcSymbolsCounter & sc1) const
{
	if (sc1.letters_in_alphabet!=letters_in_alphabet) return 0;
	if (sc1.pattern_length!=pattern_length) return 0;
	if (sc1.motifs_count!=motifs_count) return 0;
	if (sc1.total_letters!=total_letters) return 0;
	if (sc1.letters!=letters) return 0;
	if (sc1.letters_per_sequence!=letters_per_sequence) return 0;
	if (sc1.counted_as_complement!=counted_as_complement) return 0;
	if (sc1.b!=b) return 0;
	if (sc1.B!=B) return 0;
	if (sc1.c!=c) return 0;
	if (sc1.s!=s) return 0;
	return 1;
}

class DoubletAtgcSymbolsCounter:public AtgcSymbolsCounter
{
private:
	DoubletAtgcSymbolsCounter (const DoubletAtgcSymbolsCounter &);
	DoubletAtgcSymbolsCounter & operator= (const DoubletAtgcSymbolsCounter &);
	long double NegativeEntropyValue;
public:
	DoubletAtgcSymbolsCounter
	(
		const SequencesPile & sp,
		const MarkovChainState & mcs,
		unsigned short if_common_background_par,
		const vector<double> & background_opt
	);

	virtual double foreground_probability(unsigned int m,unsigned int n)
				const;
	//n-th letter in m-th motif position probability
	//position is 0-based; letter is 1-based.
	virtual double foreground_weight(unsigned int m,unsigned int n)
				const;
	//n-th letter in m-th motif position weight
	//position is 0-based; letter is 1-based.
	virtual double long NegativeEntropy () ;
	//F in Lawrence93
	//
	virtual double long NegativeEntropy_change_after_3_end_shift() ;
	//we do not use any incremantal form here, but we
	//must provide a unified interface for all counters.
	bool operator== (const DoubletAtgcSymbolsCounter & sc1) const;
};

inline
bool DoubletAtgcSymbolsCounter::operator== (const DoubletAtgcSymbolsCounter & sc1) const
{
	return AtgcSymbolsCounter::operator==(sc1);
}



#endif //_SYMBOLS_COUNTER_HPP
