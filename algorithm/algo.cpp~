/*
arguments:
1 -> formated file for first segmented protein sequence to align
2 -> formated file for second segmented protein sequence to align
3 -> threshold for important segments
4 -> word length


*/

#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <vector>
#include <algorithm>
#include "align.h"

using namespace std;

int main(int argc, char* argv[])
{
	ifstream input1 (argv[1]);
	ifstream input2 (argv[2]);
	float threshold = atof(argv[3]);
	//Setting k-word length to 2
	int k = atoi(argv[4]);
	
	//Temperory variables to read input file
	int temp_seg_no1, temp_seg_no2, temp_start, temp_end;
	float temp_weight;
	char temp_type;
	string temp_segment;
	
	vector<int> seg_no1;
	vector<char> type1;
	vector<float> weight1;
	vector<int> seg_no2;
	vector<char> type2;
	vector<float> weight2;
	vector<int> hits_seg_no1;
	vector<int> hits_seg_no2;
	vector<int> hit_size;
	vector<float> hit_weight;

	
	if (input1.is_open() && input2.is_open())
	{
		//Read first line in each line which is not the data
		getline(input1, temp_segment);
		getline(input2, temp_segment);
		
		while (1)
		{
			//Read file1
			if(!(input1 >> temp_seg_no1)) break;
			input1 >> temp_type;
			input1 >> temp_weight;
			input1 >> temp_start;
			input1 >> temp_end;
			input1 >> temp_segment;
			
			if(temp_weight >= threshold)
			{
				seg_no1.push_back(temp_seg_no1);
				type1.push_back(temp_type);
				weight1.push_back(temp_weight);
			}
			
		}
		while (1)
		{
			//Read file1
			if(!(input2 >> temp_seg_no2)) break;
			input2 >> temp_type;
			input2 >> temp_weight;
			input2 >> temp_start;
			input2 >> temp_end;
			input2 >> temp_segment;
			
			if(temp_weight >= threshold)
			{
				seg_no2.push_back(temp_seg_no2);
				type2.push_back(temp_type);
				weight2.push_back(temp_weight);
			}
			
		}

		int size1 = seg_no1.size();
		int size2 = seg_no2.size();
		int i;
				
		//Match k-words or larger
		for(temp_seg_no1 = 1; temp_seg_no1 <= size1-k+1; temp_seg_no1++)
		{

			for(temp_seg_no2 = 1; temp_seg_no2 <= size2-k+1; temp_seg_no2++)
			{
				for(i=1; i<=k; i++)
				{
					if(type1[temp_seg_no1 + i-2] != type2[temp_seg_no2 + i-2])
					{
						break;
					}
					else if (i==k)
					{
						//Storing the starting segment number of matching k-words
						hits_seg_no1.push_back(seg_no1[temp_seg_no1-1]);
						hits_seg_no2.push_back(seg_no2[temp_seg_no2-1]);
						hit_size.push_back(k);
					}
				}
			}
		}


		//merge hits. By merging we do not need to extend the hits on both sides
		for(temp_seg_no1 = 0; (temp_seg_no1) < hits_seg_no1.size(); temp_seg_no1++)
		{
			for(i=1; ((hits_seg_no1[temp_seg_no1+i]-hits_seg_no1[temp_seg_no1] <= hit_size[temp_seg_no1]) || (hits_seg_no2[temp_seg_no1+i]-hits_seg_no2[temp_seg_no1] <= hit_size[temp_seg_no1])) && (i+temp_seg_no1 < hits_seg_no1.size()); i++)
			{
				if(hits_seg_no1[temp_seg_no1+i]-hits_seg_no1[temp_seg_no1] == hits_seg_no2[temp_seg_no1+i]-hits_seg_no2[temp_seg_no1])
				{
					hit_size[temp_seg_no1] = hit_size[i + temp_seg_no1] + hits_seg_no1[temp_seg_no1+i]-hits_seg_no1[temp_seg_no1];
					hits_seg_no1.erase(hits_seg_no1.begin()+i+temp_seg_no1);
					hits_seg_no2.erase(hits_seg_no2.begin()+i+temp_seg_no1);
					hit_size.erase(hit_size.begin()+i+temp_seg_no1);
					i--;
				}
			}
		}

		size1 = hit_size.size();
		
		//Calculate weights for each segment
		for(temp_seg_no1=0; temp_seg_no1 < size1; temp_seg_no1++)
		{
			hit_weight.push_back(0);
			for(i=0; i<hit_size[temp_seg_no1]; i++)
			{
				hit_weight[temp_seg_no1] += (weight1[hits_seg_no1[temp_seg_no1]+i-1]*weight2[hits_seg_no2[temp_seg_no1]+i-1]);
			}
//						cout << hits_seg_no1[temp_seg_no1] << "\t" << hits_seg_no2[temp_seg_no1] << "\t" << hit_size[temp_seg_no1] << "\t" << hit_weight[temp_seg_no1] << endl;
		}
	

		//Decide the group of segments to be aligned. Hit with highest score or weight is to be aligned
		int final = distance(hit_weight.begin(), max_element (hit_weight.begin(), hit_weight.end()) );
		
		//Final Alignment
		cout << hits_seg_no1[final] << "\t" << hits_seg_no2[final] << "\t" << hit_size[final] << endl;
		if(hit_size.size() != 0)
		{
			align ("seq1.fasta", "seq2.fasta", hits_seg_no1[final], hits_seg_no2[final], hit_size[final], threshold);
		}
		else
		{
			cout << "Apply normal algo" << endl;
		}

	}
	else
	{
		cout << "Error in input files" << endl;
	}
	return 0;
}
