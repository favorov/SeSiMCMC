/****************************************************************************\
SeSiMCMC. Looking - for - motifs by MCMC project. (c) A. Favorov 2001-2013
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
For general describtion of the classes declared in the header, see headers.txt 
$Id: Sequences.hpp 1904 2013-07-10 01:45:50Z favorov $
\****************************************************************************/

#ifndef _SEQUENCES_HPP
#define _SEQUENCES_HPP

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include <vector>
#include <string>
#include <algorithm>

using namespace std;

typedef enum {unknown_capsmode=-1,off=0,one=1,all=2} caps_mode;

class Profile;

#include "Exception.hpp"
#include "Atgc.hpp"
//The sequence is a set unsigned shorts in range 1..letters_in_alphabet.
//It is capped with \0 at the end.
//

struct SequencesPile: public vector < vector<unsigned short> >
{
private:
	void recount_min_and_max_length();
public:
	
	string Read_Diagnostics;
		
	vector <string> names;

	vector <string> nucleotides;

	vector < vector<unsigned short> > mask;

	//we can avoid the mask array in current version,
	//because we have 0 for any masked position,
	//but it holds only for current file representation,
	//so we want to leave possibility to mask 
	//two strands independently or to use more sophisticated mask,
	//in current version, we can
	//use "mask[seq][pos]!=0" as well as for the mask test
	//*this[seq][pos]==0
	
	//it is !=if the position is frobidden to
	//be contained in a motif
	
	//we do need any construstor for actually we just
	//make the default:
	//
	/*
	SequencesPile::SequencesPile():sequences(),names(),mask{};
	*/
	vector < vector<unsigned short> > caps;
	//organised in the same way as mask,
	//nesseccary for obligatory_caps

	
	unsigned int max_length;
	unsigned int min_length;

	void add (unsigned short * seq, unsigned int len=0, string name="");
	//this routine is for fake (test) sequences
	//addition.
	//
	//seq is the sequence, len is its length, len=0 means that the routine
	//is to count the length itself and seq is 0-terminated
	//We maintain the sequence itself, not its copy!!!
	//the sequence is supposed to be allocated with allocs, not new.
	unsigned int find_median_length() const;
	unsigned int remove_garbage(unsigned int min_length_to_remain,caps_mode caps);
	//removes all sequences which length is less or equal to
	//min_length_to_remain; returns the count of the removed sequences
	void put_mask(const Profile & result,double masked_part)
		throw (DumbException)
			;
	//we rely of fact thet 0<masked_part<1
};


ostream & operator<< (ostream & os, const SequencesPile & sp) 
												throw (AtgcException); 

istream & operator>> (istream & is, SequencesPile & sp) 
												throw (AtgcException, IOStreamException); 


inline
void SequencesPile::add (unsigned short * seq, unsigned int len, string name)
//this routine is for fake (test) sequences
//addition.
//
//So, it is supposed to use void mask.
//It can be setup later by PutMask
//
//seq is the sequence, len is its length, len=0 means that the routine
//is to count the length itself and seq is 0-terminated
//We maintain the sequence itself, not its copy!!!
//the sequence is supposed to be allocated with allocs, not new.
{
	unsigned int length=0;

	string nucl_text;
	if (len) length=len;
	else 
		while(seq[length]) 
			length++ ;
	//we got length in the same sense as in Sequence.
	
	push_back(*new vector<unsigned short>);
	//new sequence added
	copy(seq,seq+length,back_inserter(back()));
	//the line was copied to new sequence
	Atgc::atgc2string(back(),nucl_text);
	nucleotides.push_back(nucl_text);
	names.push_back(name);
	mask.push_back(* new vector<unsigned short>(length));
	for (unsigned int p=0;p<length;p++)
		mask.back()[p]=0;
	caps.push_back(* new vector<unsigned short>(length));
	for (unsigned int p=0;p<length;p++)
		caps.back()[p]=1;
	//it is because it is used only for fake sequence pile
	if (size()>1)
	{
		if (length>max_length) max_length=length;
		if (length<min_length) min_length=length;
	}
	else
		min_length=max_length=length;
}
	
inline
unsigned int SequencesPile::find_median_length() const
{
	if (size()==0) return 0;
	vector<unsigned int> lengthes;
	for (vector < vector<unsigned short> >::const_iterator seqit=begin();
			seqit!=end();seqit++)
		lengthes.push_back(seqit->size());
	sort(lengthes.begin(),lengthes.end());
	if (size()%2==1)	return lengthes[(size()+1)/2];
	return (lengthes[size()/2]);
}

inline
void SequencesPile::recount_min_and_max_length()
{
	if (size()==0) {min_length=max_length=0;return;}
	iterator it=begin();
	min_length=max_length=it->size();
	while (it!=end())
	{
		min_length=min((size_t)min_length,it->size());
		max_length=max((size_t)max_length,it->size());
		it++;
	}
	return;
}

inline
unsigned int 
SequencesPile::remove_garbage(unsigned int min_length_to_remain, caps_mode capsmode)
//removes all sequences which length is less than
//min_length_to_remain; returns the count of the removed sequences
{
	unsigned int erased=0;
	vector<vector<unsigned short> >::iterator main_it=begin();
	vector<string>::iterator names_it=names.begin();
	vector<string>::iterator nucleotides_it=nucleotides.begin();
	vector<vector<unsigned short> >::iterator mask_it=mask.begin();
	vector<vector<unsigned short> >::iterator caps_it=caps.begin();
	do
	{
		unsigned int ok=1;
		if (capsmode)
		{
			ok=0;
			for (unsigned int pos=0;pos<caps_it->size();pos++)
				if ((*caps_it)[pos]) ok=1;
//if ok == 0 here, the sequense has no caps
		}
		if ((main_it->size()<min_length_to_remain) || !ok)
		{
			main_it=erase(main_it);	
			names_it=names.erase(names_it);
			nucleotides_it=nucleotides.erase(nucleotides_it);
			mask_it=mask.erase(mask_it);
			caps_it=caps.erase(caps_it);
			erased++;
		}
		else
		{
			main_it++;
			names_it++;
			nucleotides_it++;
			mask_it++;
			caps_it++;
		}
	}
	while (main_it!=end());
	recount_min_and_max_length();
	return erased;
		
}

	
#endif // _SEQUENCES_HPP
