# SRA Submission

# Copy all relevant files in submission folder
find $(/usr/bin/ls -d /group/pathogens/Alexp/Metabarcoding/dros_metabarcoding/data/*/ | grep 'CB3DR\|HLVKYDMXX' ) -maxdepth 1 -name *.fastq.gz -type f \
| grep '_D[0-9]\|_T[0-9]\|_M[0-9]\|_DM[0-9]\|_DLarv\|_NTC\|_POS\|-D[0-9]\|-T[0-9]\|-M[0-9]\|-DM[0-9]\|-DLarv\|-NTC\|-POS\|-pcrblank\|-extblank' \
| sort | uniq > to_copy.txt

# Check the difference between the copied and uncopied files
find $(/usr/bin/ls -d /group/pathogens/Alexp/Metabarcoding/dros_metabarcoding/data/*/ | grep 'CB3DR\|HLVKYDMXX' ) -maxdepth 1 -name *.fastq.gz -type f \
| sort | uniq > all_files.txt
diff --new-line-format="" --unchanged-line-format=""  all_files.txt to_copy.txt

# Copy files
cat to_copy.txt | xargs cp -t /group/pathogens/Alexp/Metabarcoding/dros_metabarcoding/data/submit

# get a list of filenames
cat to_copy.txt | grep '_R1_' |  sed 's!.*/!!' | sort | uniq > R1.txt
cat to_copy.txt | grep '_R2_' |  sed 's!.*/!!' | sort | uniq > R2.txt

# Get list of sample names & library ids
cat R1.txt | sed 's/_S[0-9].*$//g' | sed 's/^.*_//g' > sample_names.txt
cat R1.txt | sed 's/_S[0-9].*$//g' > library_ids.txt

paste -d'\t' sample_names.txt library_ids.txt R1.txt R2.txt > sra_submission_details.txt

# cleanup
rm all_files.txt 
rm to_copy.txt
rm R1.txt
rm R2.txt
rm sample_names.txt
rm library_ids.txt

# move to folder containing fastqs
/group/pathogens/Alexp/Metabarcoding/dros_metabarcoding/data/submit

# Initiate FTP
ftp -i
open ftp-private.ncbi.nlm.nih.gov

# Asks for username
subftp

# Password
w4pYB9VQ
#<password provided for me in the SRA submission portal>

# Move to the account folder (XXX part given on SRA submission portal)
cd uploads/alexpiperdesigns_gmail.com_LG3Vb4vu

# Create subfolder
mkdir drosophila_metabarcoding_submission_2

# Move there
cd drosophila_metabarcoding_submission_2

# Put files there (mput is for multiple files)
mput *.fastq.gz