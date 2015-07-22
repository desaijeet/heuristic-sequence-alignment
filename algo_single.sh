mv seq1.fasta psipred/psipred/
mv seq2.fasta psipred/psipred/
cd psipred/psipred
./runpsipred_single seq1.fasta
./runpsipred_single seq2.fasta
cd ..
cd ..
cd algorithm/
echo "Formatting the input files of the algorithm..."
./file_psipred seq1.ss2 seq1.txt
./file_psipred seq2.ss2 seq2.txt
echo "Aligning the segmented sequences..."
mv ../psipred/psipred/seq1.fasta ../../Thesis/algorithm/
mv ../psipred/psipred/seq2.fasta ../../Thesis/algorithm/
./myalgo seq1.txt seq2.txt $2 $1
echo "removing temperory files...\n"
rm temp_seq1.txt
rm temp_seq2.txt
rm alignment
mv seq1.fasta ../../Thesis/
mv seq2.fasta ../../Thesis/
mv final_align.txt ../../Thesis/
echo "final_align.txt is the output file"
