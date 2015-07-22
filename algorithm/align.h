#include <iostream>
#include <fstream>
#include <string>
#include <unistd.h>
#include <cmath>

#define lambda 0.267
#define K 0.041

using namespace std;

void align(char const* fasta_file1, char const* fasta_file2, int seq1_segment, int seq2_segment, int segment_length, float threshold)
{
	int i, j, start1, start2, end1, end2, pre_end1, pre_end2;
	string partial_seq1, partial_seq2, temp, s1, s2, a1, a2;
	float weight, score = 0.0, temp_score;
	
	ifstream input_file1 (fasta_file1);
	ifstream input_file2 (fasta_file2);
	getline(input_file1, a1);
	getline(input_file2, a2);
	input_file1 >> a1;
	input_file2 >> a2;
//	cout << a1 << endl;
//	cout << seq2 << endl;
	
	//read the segment file
	ifstream seg_file1 ("seq1.txt");
	ifstream seg_file2 ("seq2.txt");
	
	
	if (seg_file1.is_open() && seg_file2.is_open())
	{
		//Read first line in each line which is not the data
		getline(seg_file1, temp);
		getline(seg_file2, temp);


		//Sequence1
		i=1;
		while (i <= seq1_segment-1)
		{
			getline(seg_file1, temp);
			i++;
		}
		
		//Reach the intrested part in sequence2
		i=1;
		while (i <= seq2_segment-1)
		{
			getline(seg_file2, temp);
			i++;
		}
		

		ofstream align_final ("final_align.txt");
		
		for(i=1;i<=segment_length; i++)
		{
			ofstream temp_seq1 ("temp_seq1.txt");
			ofstream temp_seq2 ("temp_seq2.txt");
			if (temp_seq1.is_open() && temp_seq2.is_open())
			{
				weight = 0.0;
				s1 = "\0";
				while(weight <= threshold)
				{
					seg_file1 >> temp;
					seg_file1 >> temp;
					seg_file1 >> weight;
					seg_file1 >> start1;
					seg_file1 >> end1;
					seg_file1 >> temp;
					if(weight < threshold)
					{
						s1 = s1 + temp;
					}
				}
				cout << temp << endl;
				temp_seq1 << temp;
		
				//Sequence2
				weight = 0.0;
				s2 = "\0";
				while(weight <= threshold)
				{
					seg_file2 >> temp;
					seg_file2 >> temp;
					seg_file2 >> weight;
					seg_file2 >> start2;
					seg_file2 >> end2;
					seg_file2 >> temp;
					if(weight < threshold)
					{
						s2 = s2 + temp;
					}
				}
				cout << temp << endl;
				temp_seq2 << temp;
				
				
				//Aligning in between sequences
				if(s1.length() != 0 && s2.length() == 0)
				{
					for(j=1; j<=s1.length(); j++)
					{
						s2 = s2 + "-";
					}
					align_final << "\n" << pre_end1+1 <<  "\t" << s1 << "\t" << start1-1 << endl;
					align_final << "\t" << s2 << "\t" << endl;
				}
				else if(s1.length() == 0 && s2.length() != 0)
				{
					for(j=1; j<=s2.length(); j++)
					{
						s1 = s1 + "-";
					}
					align_final << "\n" <<  "\t" << s1 << "\t" << endl;
					align_final << pre_end2+1 << "\t" << s2 << "\t" << start2-1 << endl;
				}
				
				else if (s1.length() != 0 && s2.length() != 0)
				{
					ofstream a1 ("temp1.txt");
					ofstream a2 ("temp2.txt");
					a1 << s1;
					a2 << s2;
					a1.close();
					a2.close();
					ifstream fs1 ("temp1.txt");
					ifstream fs2 ("temp2.txt");
					system("/usr/local/emboss/bin/needle -asequence temp1.txt -bsequence temp2.txt -gapopen 10.0 -gapextend 0.5 -outfile alignment -datafile EBLOSUM62 -awidth 200");
					
					//Reading the output file
					ifstream align_temp1 ("alignment");
					for(j=1;j<=32;j++)
					{
						getline(align_temp1, temp);
					}
					
					align_temp1 >> temp;
					align_temp1 >> temp;
					align_final << "\n" << pre_end1+1 << "\t" << temp << "\t" << start1-1 << endl;
					
					//For sequence2
					getline(align_temp1, temp);
					getline(align_temp1, temp);
					align_temp1 >> temp;
					align_temp1 >> temp;
					align_final << pre_end2+1 << "\t" << temp << "\t" << start2-1 << endl;

					align_temp1.close();
					fs1.close();
					fs2.close();
				}

			
				temp_seq1.close();
				temp_seq2.close();
			}
			
			ifstream seq1 ("temp_seq1.txt");
			ifstream seq2 ("temp_seq2.txt");
			if (seq1.is_open() && seq2.is_open())
			{
				
				system("/usr/local/emboss/bin/needle -asequence temp_seq1.txt -bsequence temp_seq2.txt -gapopen 10.0 -gapextend 0.5 -outfile alignment -datafile EBLOSUM62 -awidth 200");
				//j = execl("/usr/local/emboss/bin/needle", "needle", "-asequence", "temp_seq1.txt", "-bsequence", "temp_seq2.txt", "-gapopen", "10.0", "-gapextend", "0.5", "-outfile", "alignment", "-datafile", "EBLOSUM62", NULL);

				ifstream align_temp ("alignment");
//				getline(align_temp, temp);
				for(j=1;j<=32;j++)
				{
					if(j == 29)
					{
						align_temp >> temp;
						align_temp >> temp;
						align_temp >> temp_score;
						score += temp_score;
						
					}
					getline(align_temp, temp);
				}
					
				align_temp >> temp;
				align_temp >> temp;
				align_final << "\n" << start1 << "\t" << temp << "\t" << end1 << endl;
				
				//For sequence2
				getline(align_temp, temp);
				getline(align_temp, temp);
				align_temp >> temp;
				align_temp >> temp;
				align_final << start2 << "\t" << temp << "\t" << end2 << endl;
				
				pre_end1 = end1;
				pre_end2 = end2;
				
				align_temp.close();
			}
			seq1.close();
			seq2.close();
			

		}		
		align_final.close();
		cout << "Score for the alignment is " << score << endl;
		cout << a1.length() << "\t" << a2.length() << endl;

		//Bit score calculation
		cout << "Bit-Score for the alignment is " << (lambda*score - log(K))/log(2) << endl << endl;

		//Calculation of E-value
		cout << "E-value is " << (K * a1.length() * a2.length() * exp(-lambda*score)) << endl;
		cout << "E-value from Bit score is " << (a1.length() * a2.length()) / pow(2, (lambda*score - log(K))/log(2)) << endl << endl;

	}
	else
	{
		cout << "Error in input files in align.h" << endl;
	}
	
	seg_file1.close();
	seg_file2.close();

//	ofstream needleman_input1 ("");
//	ofstream needleman_input2 ("");
	
}
