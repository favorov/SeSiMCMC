/****************************************************************************\
SeSiMCMC. Looking - for - motifs by MCMC project. (c) A. Favorov 2001-2019
$Id$
\****************************************************************************/

#ifndef _MCMC_HPP
#define _MCMC_HPP

//#define __DEBUGLEVEL__ 2
//

#include <math.h>

#include <algorithm>
#include <vector>
#include <set>

#include <iostream>

using namespace std;

#include "MarkovChainState.hpp"
#include "Logger.hpp"
#include "Exception.hpp"
#include "Sequences.hpp"
#include "SymbolsCounter.hpp"
#include "Diagnostics.hpp"



const double motif_significance_level=0.05;

typedef enum {unknown_mode=-1,one_thread=1,two_threads} complement_mode;
typedef enum {unknown_smode=-1,no=0,repeats,palindromes} symmetry_mode;
typedef enum {unknown_ctmode=-1,site_positions=0,information=1,correlation=2} change_test_mode;
typedef enum {unknown_ptmode=-1, permut_never=0,permut_chainwise=1,permut_adjustmentwise=2,permut_cyclewise=3} permutation_mode;


inline ostream & operator<<(ostream &o, caps_mode cm)
{
	switch (cm)
	{
	case off: o<<"off";break;
	case one: o<<"one";break;
	case all: o<<"all";break;
	default: o<<"unknown";
	}
	return o;
}

/*
 * All the mathemtical idea of the sampler is taken from
 * Lawrence93, i.e
 *
 * Lawrence, Altschul, Boguski, Liu, Neuwald and Wootton (1993)
 * Detecting Subtle Sequence Signals: A Gibbs Sampling
 * Strategy for Multiple Alignment, Science 262:208-214
 *
 * There are two kinds of sampler activity: LocalGibbsStep and PositionAdjustment.
 * Here is the description of LocalGibbsStep.
 * One pattern is current of every step.
 * From all other patterns, we count probability for
 * a letter j to be in position i in a pattern and
 * probability p for it to be in nonpattern space.
 *
 * Than, we count a_i=mult(q_i_j/p_i) along all letters of the pattern.
 *
 * Than, we create y_i=non-normalised vector of a_i.
 *
 * Calculation of Y_i is made by the routine calculate_position_weights()
 *
 * Then, we generate the new pattern i in current sequence from y_i vector by
 * get_int_draw_from_a_ditribution().
 *
 * One of PositionAdjustment or PositionAndLengthAdjustment (we will call it
 * Adjustment) is made once in  :
 * 		(sequence_number*local_cycles_between_global_adjustments),
 * i.e after PositionAdjustment we make local_cycles_between_global_adjustments
 * cycles of LocalGibbsStep and then make the Adjustment.
 * Both Adjustment and LocalGibbsStep shifts the current pattern,
 * so Adjustment happens  every time with new current_sequence.
 *
 * local_cycles_between_global_adjustments=0 in constructor means that we
 * count automatically.
 *
 * Is counts the weight af all possible shifts of all patterns except current
 * sticked in a as a solid bunch and then generates the shift of the bunch from
 * the weights. We just maximise the weight by shifting the bunch.
 * It is more effective, but we cannot gather any statistics other
 * than between two PositionAdjustment. When we adjust length, we count the
 * best shift for every length and than compare Information_per_pattern_position
 * for every length, finding the best.
 *
 * Here, we write sampler devoted to DNA text sampling (we always suggest that
 * compements do make sense). To test protein sequences, we can put line
 * complements_are_possible=0 to configuration file.
 *
 *
 * The sampler has two stages of work : annealing (looking for smth sesible)
 * which is done only by local steps and maximising where we suppose
 * nothins about the motifs, they can find their lenght and presence them-
 * selves.
 * The test for annealing is over is: scalar products of the states before
 * and after  annealing_is_over_test_cycles (default 2) sequential cycles
 * of local samling are not less than for sequential_states_product
 * (default 0.5).
 */

/*
 * LookingForAtgcMotifsMultinomialGibbs is an abstract base
 * LookingForAtgcMotifsInOneThreadMultinomialGibbs
 * LookingForAtgcMotifsInTwoThreadsMultinomialGibbs,
 * LookingForAtgcPalindromesMultinomialGibbs
 * which are slightly different in their numerical
 * functions and works with generic motifs in one direction,
 * generic motifs in both directions and palindromic motifs,
 * correspondingly. They all are based on
 * LookingForAtgcMotifsMultinomialGibbs
 *
 * LookingForAtgcRepeatsInOneThreadMultinomialGibbs and
 * LookingForAtgcRepeatsInTwoThreadsMultinomialGibbs are based on
 * LookingForAtgcMotifsInOneThreadMultinomialGibbs
 * LookingForAtgcMotifsInTwoThreadsMultinomialGibbs, correspondingly.
 * They looks for paired repeats instead of asymmetrical motifs
 */

class LookingForAtgcMotifsMultinomialGibbs
{
private:
	LookingForAtgcMotifsMultinomialGibbs
			(const LookingForAtgcMotifsMultinomialGibbs &);
	//copy constructor is forbidden
	LookingForAtgcMotifsMultinomialGibbs & operator=
			(const LookingForAtgcMotifsMultinomialGibbs &);
	//copy operator is forbidden


protected:
	vector <double> Y;  //the vector to keep position weigths.

	const SequencesPile & Sequences;

	unsigned int letters_in_alphabet;

	SymbolsCounter * Symbols;

	unsigned int sequences_count;

	complement_mode mode;
	symmetry_mode smode;
	//it carries information about what kind of sampler do we derive,
	//sometimes it is useful not to duplicate code which is almost
	//identical.

	unsigned short last_adjustment_has_changed_state;

	unsigned int current_sequence_index;
  //0..Sequences.count

  unsigned short be_quiet;

	caps_mode current_caps_mode;

	unsigned int get_int_draw_from_a_ditribution
			(
				unsigned int least_possible_draw,
				unsigned int most_possible_draw,
				vector<double> y
				//get the distribution from here.
				//the vector it to be at least
				//last_possible_position-first_possible_position+1
			);

	virtual void calculate_position_weights
			(
				const MarkovChainState & state,
			 	unsigned int current_sequence_no,
				unsigned int motif_positions,
				SymbolsCounter & symbols,
				//use this SymbolsRate for count normalised weights,
				vector<double> & y
				//place the results here.
				//the vector it to be at least
				//last_possible_position-first_possible_position+1
			 )=0;

	virtual void LocalGibbsStep(unsigned int sequence_no)=0;


	unsigned int ShiftOptimisation
			(
				MarkovChainState & state,
				unsigned int shift_magnitude,
				unsigned short & success,
				//success == 0 if all possibilities are masked
				long double & max_entropy,
				long double & G
			);
	//we return 1 if we have changed the state

	unsigned int ShiftAndGapOptimisation
			(
				MarkovChainState & state,
				unsigned int shift_magnitude,
				unsigned short & success,
				//success == 0 if all possibilities are masked
				long double & F,
				long double & best_G
			);
	//we return 1 if we have changed the state

	long double PositionAdjustment
			(
				unsigned long local_steps_made,
				unsigned long adjustments_made,
				unsigned short test_gaps
			);

	long double PositionAndLengthAdjustment
			(
				unsigned int minimal_length,
				unsigned int maximal_length,
				unsigned short test_gaps,
				unsigned long local_steps_made,
				unsigned long adjustments_made
			);
	//they return the resulting
	//Negative_Entropy_with_defined_patterns_per_pattern_position.
	//It is useful not to recount it.
	virtual void NormaliseState()=0;
	//it is necessary only for TwoThread,
	//but we describe it as an empty one here to
	//compile code with its call.
public:
	class SlowSearchFailed
	{
	public:
		SlowSearchFailed(){};
	};

	class LostInSpaceOnSecondaryAnnealing  //exception thrown if max_steps is exceeded
	{
	public:
		unsigned long max_steps;
		unsigned int length;
		LostInSpaceOnSecondaryAnnealing(unsigned long ms,unsigned int len):
			max_steps(ms),length(len){};
	};

	class AllInitialAnnealingAttemptsFailed  //exception thrown if max_steps is exceeded
	{
	public:
		unsigned int attempts,length;
		AllInitialAnnealingAttemptsFailed(unsigned int att,unsigned int len):
			attempts(att),length(len){};
	};

	class LostInSpaceOnTracing  //exception thrown if max_steps is exceeded
	{
	public:
		unsigned long max_steps;
		unsigned int lenght;
		double long maxG;
		LostInSpaceOnTracing(unsigned long ms, unsigned int len, double long g):
			max_steps(ms),lenght(len),maxG(g){};
	};

	class TooMuchChainsFailed  //exception thrown if max_steps is exceeded
	{
	public:
		unsigned int chains,lenght;
		double long maxG;
		TooMuchChainsFailed(unsigned int ch,unsigned int len, double long g):
			chains(ch),lenght(len),maxG(g){};
	};

	class AllShiftsAreMasked  //exception thrown if all shifts are masked
	{
	public:
		const char * info;
		AllShiftsAreMasked(const char* s=""):info(s){};
	};

	class TimeLimitException
	{
	public:
		unsigned int time_limit;
		TimeLimitException(unsigned int limit):time_limit(limit){};
	};

	class TooResrtrictiveCapsMode
	{
	public:
		unsigned int sequences_OK,sequences;
		TooResrtrictiveCapsMode(unsigned int seqs_OK,unsigned int seqs):
			sequences_OK(seqs_OK),sequences(seqs){};
	};

	class LenghtRequirementsFailed
	{
	public:
		unsigned int sequences_OK,sequences;
		LenghtRequirementsFailed(unsigned int seqs_OK,unsigned int seqs):
			sequences_OK(seqs_OK),sequences(seqs){};
	};

	MarkovChainState the_state;
	//we need it to be public to create additional counters if we need
	virtual SymbolsCounter * create_new_counter
			(
			 	const SequencesPile & sp,
				const MarkovChainState & mcs
			) const =0;
	double motif_absence_prior;

	caps_mode initial_caps;
	caps_mode annealing_caps;
	caps_mode refinement_caps;
	unsigned short if_adaptive_pseudocounts;
	unsigned short if_common_background;
	vector <double> background;
	unsigned int default_pattern_lenght;
	unsigned int shortest_significant_motif;
	unsigned int longest_sensible_motif;

	unsigned long local_step_cycles_between_adjustments;
	//we mean count of times we hit one sequence locally

	unsigned long overall_local_steps_made;
	unsigned long overall_adjustments_made;
	//it is not the best idea to make these two pubclic,
	//but we never use it inside, it is necessary only for outworld.
	//So, if outworld changes it, it does it for itself.

	ostream & log_stream;  //ouptut log file, currently cerr.
	const MarkovChainState & theState;
	                       //still all the state information is the
												 //positions, it is a const ref to the_state
												 //

	LookingForAtgcMotifsMultinomialGibbs
	(
		const SequencesPile & Sequences,
		double motif_absence_prior,
		unsigned long local_step_cycles_between_adjustments,
		caps_mode initial_caps,
		caps_mode annealing_caps,
		caps_mode refinement_caps,
		unsigned short if_adaptive_pseudocounts,
		unsigned short if_common_background,
		const vector<double> & background_opt,
		unsigned short be_quiet,
		ostream & logstream
	);

	unsigned int reinit(unsigned int len);

	double long find_maximum
			(
				unsigned int pattern_length,
				unsigned int minimal_pattern_length,
				unsigned int maximal_pattern_length,
				unsigned short adjust_pattern_length,
				unsigned short spaced,
			 	unsigned long steps_number_maximum_is_to_be_global_for,
				unsigned int annealings_number_maximum_is_to_be_global_for,
				unsigned int cycles_with_minor_change_to_say_cold,
			 	double minor_change_level,
				change_test_mode ctmode,
				unsigned int adjustments_during_annealing,
				unsigned int chain_fails_after,
				unsigned int chains_to_try,
				unsigned long cycles_per_annealing_attempt,
				unsigned int annealing_attempts,
				//unsigned int annealings_with_different_length,
				unsigned long time_limit,
				const LogRecorder & timer,
				//unsigned int reanneal_if_all_chains_fail=0,
				unsigned long max_steps=0ul-1
			);

	double find_maximum_slowly
			(
				unsigned int minimal_pattern_length,
				unsigned int maximal_pattern_length,
				unsigned short spaced,
				set <pair<unsigned int,double> > & G_values,
			 	unsigned long steps_number_maximum_is_to_be_global_for,
				unsigned int annealings_number_maximum_is_to_be_global_for,
				unsigned int cycles_with_minor_change_to_say_cold,
			 	double minor_change_level,
				change_test_mode ctmode,
				unsigned int adjustments_during_annealing,
				unsigned int chain_fails_after,
				unsigned int chains_to_try,
				unsigned long cycles_per_annealing_attempt,
				unsigned int annealing_attempts,
				unsigned long time_limit,
				const LogRecorder & timer,
				Diagnostics & diags,
				unsigned long max_steps=0ul-1
			);

	double long Negative_Entropy_with_defined_patterns_per_pattern_position
	//it never changes Symbols, but it supposes it to be
	//up-to-date (for Symbols.NegativeEntropy() call
	(
		const MarkovChainState & state,
		double long & NegativeEntropy
		//it is. to return NegativeEntropy of Symbols not to recount it.
	);
	void SayLengthBounds
		 (
				unsigned int & conf_minimal_pattern_length,
				unsigned int & conf_pattern_length,
				unsigned int & conf_maximal_pattern_length
		 );
	//we return the bound which will be actually used
	//when calling a procedure with
	//minimal_pattern_length, maximal_pattern_length
	//parameters.
	unsigned int permutation_forbidden; //test workaaround

	virtual ~LookingForAtgcMotifsMultinomialGibbs(){};
};

inline
void LookingForAtgcMotifsMultinomialGibbs::SayLengthBounds
		 (
				unsigned int & conf_minimal_pattern_length,
				unsigned int & conf_pattern_length,
				unsigned int & conf_maximal_pattern_length
		 )
{
	if (conf_minimal_pattern_length==0) //use default
		conf_minimal_pattern_length=shortest_significant_motif;
	else
		conf_minimal_pattern_length=
				max(shortest_significant_motif,conf_minimal_pattern_length);
	if (conf_maximal_pattern_length==0) //use default
		conf_maximal_pattern_length=longest_sensible_motif;
	else
		conf_maximal_pattern_length=
				min(longest_sensible_motif,conf_maximal_pattern_length);

	unsigned int length=
		conf_pattern_length?conf_pattern_length:default_pattern_lenght;

	length=min(length,conf_maximal_pattern_length);
	length=max(length,conf_minimal_pattern_length);

	conf_pattern_length=length;
	default_pattern_lenght=length;
}


inline
ostream & operator<<
		(ostream & o, const LookingForAtgcMotifsMultinomialGibbs::SlowSearchFailed & lost)
{
	o<<"The slow mode search failed.\n";
	return o;
}

inline
ostream & operator<<
		(ostream & o, const
				LookingForAtgcMotifsMultinomialGibbs::LostInSpaceOnSecondaryAnnealing & lost)
{
	o<<"The sampler was not secondary annealed after "<<
		lost.max_steps<< " steps made. The initial length was "<<lost.length<<".\n";
	return o;
}

inline
ostream & operator<<
		(ostream & o, const LookingForAtgcMotifsMultinomialGibbs::AllInitialAnnealingAttemptsFailed & lost)
{
	o<<"All of "<<lost.attempts<<" initial annealing attempts failed, initial pattern length was "<<
			lost.length<<".\n";
	return o;
}

inline
ostream & operator<<
		(ostream & o, const LookingForAtgcMotifsMultinomialGibbs::LostInSpaceOnTracing & lost)
{
	o<<"The sampler did not find a good maximum after "<<
			lost.max_steps<< " steps made. Initial pattern length was "<<lost.lenght<<".\n";
	return o;
}


inline
ostream & operator<<
		(ostream & o, const LookingForAtgcMotifsMultinomialGibbs::TooMuchChainsFailed & toomuch)
{
	o<<"The sampler has failed by failing all the "<<
			toomuch.chains<< " chains. Initial pattern length was "<<toomuch.lenght<<".\n";
	return o;
}

inline
ostream & operator<< (ostream & o, const LookingForAtgcMotifsMultinomialGibbs::AllShiftsAreMasked & asam)
{
	o<<asam.info;
	return o;
}

inline
ostream & operator<< (ostream & o, const LookingForAtgcMotifsMultinomialGibbs::TimeLimitException & tle)
{
	o<<"The program was halted because\n"<<
	"the given time limit of "<<tle.time_limit<<" seconds is used.\n"<<
	 "You can use the standalone command-line version with\n"<<
	 "--time-limit 0 option that switch off time limits.\n\n";
	return o;
}

inline
ostream & operator<< (ostream & o, const LookingForAtgcMotifsMultinomialGibbs::LenghtRequirementsFailed & lrf)
{
	o<<"The program was stopped because\nthere are ";
	if (lrf.sequences_OK)
		o<<"only "<<lrf.sequences_OK;
	else
		o<<"no";
	o<<" sequences that can meet the caps requirements"<<endl<<"among "
   <<lrf.sequences<<" input sequences.\n"<<"It looks like something is wrong.\n";
	return o;
}

inline
ostream & operator<< (ostream & o, const LookingForAtgcMotifsMultinomialGibbs::TooResrtrictiveCapsMode & trcm)
{
	o<<"The program was stopped because\nthere are ";
	if (trcm.sequences_OK)
		o<<"only "<<trcm.sequences_OK;
	else
		o<<"no";
	o<<" sequences that can meet the lenght requirements"<<endl<<"among "
   <<trcm.sequences<<" input sequences.\n"<<"It looks like something is wrong.\n";
	return o;
}


class LookingForAtgcMotifsInOneThreadMultinomialGibbs:
	public LookingForAtgcMotifsMultinomialGibbs
{

private:

	virtual void calculate_position_weights
			(
				const MarkovChainState & state,
			 	unsigned int sequence_no,
				unsigned int motif_positions,
				SymbolsCounter & symbols,
				//use this SymbolsRate for count normalised weights,
				vector<double> & y
				//place the results here.
				//the vector it to be at least
				//last_possible_position-first_possible_position+1
			 );

//the result is written to y;
//y[0] is propotional to prob. of motif absence
//y[1..motif_positions] are proportinal to
//probabilities of motif present in the position

	virtual void LocalGibbsStep(unsigned int sequence_no);
	//they return the resulting
	//Negative_Entropy_with_defined_patterns_per_pattern_position.
	//It is useful not to recount it.
	virtual void NormaliseState();
public:
	virtual SymbolsCounter * create_new_counter
			(
			 	const SequencesPile & sp,
				const MarkovChainState & mcs
			) const ;

	LookingForAtgcMotifsInOneThreadMultinomialGibbs
	(
		const SequencesPile & Sequences,
		double motif_absence_prior,
		unsigned long local_step_cycles_between_adjustments,
		caps_mode initial_caps,
		caps_mode annealing_caps,
		caps_mode refinement_caps,
		unsigned short if_adaptive_pseudocounts,
		unsigned short if_common_background,
		const vector<double> & background_opt,
		unsigned short be_quiet,
		ostream & logstream=cerr
	);


};

class LookingForAtgcMotifsInTwoThreadsMultinomialGibbs:
	public LookingForAtgcMotifsMultinomialGibbs
{

private:

	virtual void calculate_position_weights
			(
				const MarkovChainState & state,
			 	unsigned int sequence_no,
				unsigned int motif_positions,
				SymbolsCounter & symbols,
				//use this SymbolsRate for count normalised weights,
				vector<double> & y
				//place the results here.
				//the vector it to be at least
				//last_possible_position-first_possible_position+1
			 );
//the result is written to y;
//y[0] is propotional to prob. of motif absence
//y[1..motif_positions] are proportinal to
//probabilities of motif present in the position+1 (pos is 0-based)
//y[motif_positions+1..motif_positions+motif_positions] - to
//probabilities of the comlementary presence of the motif
//is position (the corrresponding index is motif_positions+position)
//and y[2*motif_positions+1] is prob. of complementary reading of
//the sequence and motif absence in it.
/*
    nc   pos+1      |          mp+pos+1             comp
0   | 1   .....  mp |  mp+1 | mp+2 |         | mp+mp | mp+mp+1
abs      present
*/

	virtual void LocalGibbsStep(unsigned int sequence_no);
	//they return the resulting
	//Negative_Entropy_with_defined_patterns_per_pattern_position.
	//It is useful not to recount it.
	void BalanceComplements();
	//it is necessary only for TwoThread,
	//but we describe it as an empty one here to
	//compile code with its call.
	virtual void NormaliseState();

public:
	virtual SymbolsCounter * create_new_counter
			(
			 	const SequencesPile & sp,
				const MarkovChainState & mcs
			) const ;

	LookingForAtgcMotifsInTwoThreadsMultinomialGibbs
	(
		const SequencesPile & Sequences,
		double motif_absence_prior,
		unsigned long local_step_cycles_between_adjustments,
		caps_mode initial_caps,
		caps_mode annealing_caps,
		caps_mode refinement_caps,
		unsigned short if_adaptive_pseudocounts,
		unsigned short if_common_background,
		const vector<double> & background_opt,
		unsigned short be_quiet,
		ostream & logstream=cerr
	);

};

class LookingForAtgcPalindromesMultinomialGibbs:
	public LookingForAtgcMotifsMultinomialGibbs
{

private:
	virtual SymbolsCounter * create_new_counter
			(
			 	const SequencesPile & sp,
				const MarkovChainState & mcs
			) const;

	virtual void calculate_position_weights
			(
				const MarkovChainState & state,
			 	unsigned int current_sequence_no,
				unsigned int motif_positions,
				SymbolsCounter & symbols,
				//use this SymbolsRate for count normalised weights,
				vector<double> & y
				//place the results here.
				//the vector it to be at least
				//last_possible_position-first_possible_position+1
			 );

//the result is written to y;
//y[0] is propotional to prob. of motif absence
//y[1..motif_positions] are proportinal to
//probabilities of motif present in the position

	virtual void LocalGibbsStep(unsigned int sequence_no);
  //they return the resulting
	//Negative_Entropy_with_defined_patterns_per_pattern_position.
	//It is useful not to recount it.
	virtual void NormaliseState();
public:
	LookingForAtgcPalindromesMultinomialGibbs
	(
		const SequencesPile & Sequences,
		double motif_absence_prior,
		unsigned long local_step_cycles_between_adjustments,
		caps_mode initial_caps,
		caps_mode annealing_caps,
		caps_mode refinement_caps,
		unsigned short if_adaptive_pseudocounts,
		unsigned short if_common_background,
		const vector<double> & background_opt,
		unsigned short be_quiet,
		ostream & logstream=cerr
	);

};

class LookingForAtgcRepeatsInOneThreadMultinomialGibbs:
	public LookingForAtgcMotifsInOneThreadMultinomialGibbs
{

private:

public:
	virtual SymbolsCounter * create_new_counter
			(
			 	const SequencesPile & sp,
				const MarkovChainState & mcs
			) const ;

	LookingForAtgcRepeatsInOneThreadMultinomialGibbs
	(
		const SequencesPile & Sequences,
		double motif_absence_prior,
		unsigned long local_step_cycles_between_adjustments,
		caps_mode initial_caps,
		caps_mode annealing_caps,
		caps_mode refinement_caps,
		unsigned short if_adaptive_pseudocounts,
		unsigned short if_common_background,
		const vector<double> & background_opt,
		unsigned short be_quiet,
		ostream & logstream=cerr
	);
};

class LookingForAtgcRepeatsInTwoThreadsMultinomialGibbs:
	public LookingForAtgcMotifsInTwoThreadsMultinomialGibbs
{

private:

public:
	virtual SymbolsCounter * create_new_counter
			(
			 	const SequencesPile & sp,
				const MarkovChainState & mcs
			) const ;

	LookingForAtgcRepeatsInTwoThreadsMultinomialGibbs
	(
		const SequencesPile & Sequences,
		double motif_absence_prior,
		unsigned long local_step_cycles_between_adjustments,
		caps_mode initial_caps,
		caps_mode annealing_caps,
		caps_mode refinement_caps,
		unsigned short if_adaptive_pseudocounts,
		unsigned short if_common_background,
		const vector<double> & background_opt,
		unsigned short be_quiet,
		ostream & logstream=cerr
	);
};


#endif //_MCMC_HPP
