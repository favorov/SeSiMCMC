/****************************************************************************\
SeSiMCMC. Looking - for - motifs by MCMC project. (c) A. Favorov 2001-2019
APSampler project. (c) A. Favorov 1999-2007
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
For general describtion of the classes declared in the header, see headers.txt 
$Id$
\****************************************************************************/

#ifndef _RANDOM_MAPPINGS_HPP
#define _RANDOM_MAPPINGS_HPP

using namespace std; 

#include <vector>
#include <iostream>
#include <iomanip>
#include <iterator>
#include <string>
#include <math.h>

extern "C" {
	#include "Random.h"
}


template <class T> //T is any container of unsigned int with [] operation
inline const T & RandomPermutationMapping (T & container, unsigned int start,
		unsigned int end)
//we will permutted integers in range [start..end-1]
//and put it to first end-start elements in the container
{
	if (container.size()<end-start)
	{
		container.resize(end-start);
	}
	
	T scattering(end-start);

	for (unsigned int drawn_item_no=0;
			drawn_item_no<end-start;
			drawn_item_no++)
	{
		unsigned int u;
		unsigned int outcome=
				(unsigned int)floorf((float)(end-start-drawn_item_no)*uni());
		//int in [0,end-start-item_no-1]

		//it is mapping from [0,end-start-drawn_item_no-1]
		//to ( [0,end-start] \ {set of all previously drawn values}
		//the value is i,
		//the map(i) is i+n,
		//where n is the number of already 
		//drawn elements that are less than map(i) 
		//(looks like an equation :) )
	    //the mapping requires that all the the already drawn 
		//items are sorted in a list. Then
		//we start to look for the place in the list is the usual sorting
	  	//manner, but invreasing the value each time the pointer moves up

		
		for (u=0;u<drawn_item_no;u++) 	
			//the previous are ascendingly sorted
			if(scattering[u]<=outcome) outcome++; 
			else 
			//we found the place for the (mapped) outcome in the
			//sorted set of values 
			{
				for(unsigned u_plus=drawn_item_no;u_plus>u;u_plus--)
					scattering[u_plus]=scattering[u_plus-1]; 
				//shift all the previous values that are > outcome 
				//to 1 position right
				scattering[u]=outcome;
				break;
			}
		if (u==drawn_item_no)	
		{
			scattering[u]=outcome; //it was the largest.
		}
		container[drawn_item_no]=outcome+start;
	}
/*	for (unsigned int item_no=0;
			item_no<end-start;
			item_no++)
	{
		unsigned int u;
		unsigned int outcome=
				(unsigned int)floorf((float)(end-start-item_no)*uni());
		//int in [0,end-start-item_no-1]
		outcome+=start;
		unsigned int mapping_shift=0;
		for (u=0;u<item_no;u++)
		//the previous are ascendingly sorted
			if(container[u]<=outcome) mapping_shift++;
			//for each prevuiosly drawn value that is < outcome,
			//we increase the oucome for 1
			//it is mapping from [0,range-drawn_item_no-1]
			//to ( [0,range-1] \ {set of all previously drawn values}
		container[item_no]=outcome+mapping_shift; //it was the largest.
	}*/
	return container;
};

template <class T> //T is any container of unsigned int with [] operation
inline const T & RandomPermutationMapping (T & container, unsigned int end)
//we will permutted integers in range [0..end-1]
//and put it to first end-1 elements in the container
{
	return RandomPermutationMapping(container,0,end);
};

template <class T> //T is any container of unsigned int with [] and size ()operation
inline const T & RandomPermutationMapping (T & container)
//we will permut integers in range [0..container.size()-1]
//and put it in the container
{
	return RandomPermutationMapping(container,0,container.size());
};

template <class V> //V is any container of anything with [] and size ()operation
inline const V & RandomPermutation (V & container)
//we will permut the container
{
	vector<unsigned int> mapping(container.size());
	RandomPermutationMapping (mapping);
	V container_copy(container);
	cout<<"MAP:\n";
	copy(mapping.begin(), mapping.end(), ostream_iterator<unsigned int>(cout," "));
	cout<<"\n";
	for (unsigned int i=0;i<container.size();i++)
		container.at(i)=container_copy[mapping[i]];
	return container;
};


template <class T> //T is any container of unsigned int with [] operation
inline const T & RandomDrawMapping (T & scattering, unsigned int range,
		unsigned int items_number)
//we will put in T[0]..T[items_number-1]
//different random draws from [0..range-1]
//and sort them
//range is WHERE FROM (values range) 
//items_number is WHERE TO (destination capacity)
//normally, range>items_number
//range=items_number is a coorect operation, but it is insensible
{
	if (scattering.size()<items_number)
	{
		scattering.resize(items_number);
	}
	
	unsigned int items_to_draw=0;
	items_to_draw=range>items_number?items_number:range;
	
	for (unsigned int drawn_item_no=0;
			drawn_item_no<items_to_draw;
			drawn_item_no++)
	{
		unsigned int u;
		unsigned int outcome=
				(unsigned int)floorf((float)(range-drawn_item_no)*uni());
		//int in [0,range-drawn_item_no-1]

		//it is mapping from [0,range-drawn_item_no-1]
		//to ( [0,range-1] \ {set of all previously drawn values}
		//the value is i,
		//the map(i) is i+n,
		//where n is the number of already 
		//drawn elements that are less than map(i) 
		//(looks like an equation :) )
	    //the mapping requires that all the the already drawn 
		//items are sorted in a list. Then
		//we start to look for the place in the list is the usual sorting
	  	//manner, but invreasing the value each time the pointer moves up

		
		for (u=0;u<drawn_item_no;u++) 	
			//the previous are ascendingly sorted
			if(scattering[u]<=outcome) outcome++; 
			else 
			//we found the place for the (mapped) outcome in the
			//sorted set of values 
			{
				for(unsigned u_plus=drawn_item_no;u_plus>u;u_plus--)
					scattering[u_plus]=scattering[u_plus-1]; 
				//shift all the previous values that are < outcome 
				//to 1 position right
				scattering[u]=outcome;
				break;
			}
		if (u==drawn_item_no)	
		{
			scattering[u]=outcome; //it was the largest.
		}
	}
	if (items_number>items_to_draw)
		for(unsigned int item=items_to_draw+1;item<items_number;item++) 
			scattering[item]=range; //impossible draw
	return scattering;
};

template <class T> //T is any container of unsigned int with [] and size() operation
inline const T & RandomDrawMapping (T & scattering,
		unsigned int items_number)
//we will put in T[0]..T[T.size()-1]
//different random draws from [0..range-1]
{
	return RandomDrawMapping(scattering,items_number,scattering.size());
};

template <class T> //T is any container of unsigned int with [] operation
inline const T & RandomDrawWithReturnMapping (T & scattering,
		unsigned int items_number, unsigned int range)
//we will put in T[0]..T[num-1]
//random draws from [0..range-1]
{
	if (scattering.size()<items_number)
	{
		scattering.resize(items_number);
	}
	for (unsigned int drawn_item_no=0;
			drawn_item_no<items_number;
			drawn_item_no++)
	{
		unsigned int u;
		unsigned int outcome=
				(unsigned int)floorf((float)(range)*uni());
		//int in [0,range-1]
		for (u=0;u<drawn_item_no;u++) //the previous are ascendingly
																								//sorted
			if(scattering[u]>outcome)
			{
				for(unsigned u_plus=drawn_item_no;u_plus>u;u_plus--)
					scattering[u_plus]=scattering[u_plus-1];
				scattering[u]=outcome; //put it sorted
				break;
			}
		if (u==drawn_item_no)	
		{
			scattering[u]=outcome; //it was the largest.
		}
	}
	return scattering;
};

template <class T> //T is any container of unsigned int with [] and size() operation
inline const T & RandomDrawWithReturnMapping (T & scattering,
		unsigned int items_number)
//we will put in T[0]..T[num-1]
//random draws from [0..range-1]
{
	return RandomDrawWithReturnMapping(scattering,items_number,scattering.size());
};

#endif
