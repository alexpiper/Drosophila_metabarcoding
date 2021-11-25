# Remove the unused miseq run
rm -r CJKFJ/


####### Renaming for CB3DR
cd /group/pathogens/Alexp/Metabarcoding/dros_metabarcoding/data/CB3DR

# Remove any Carpophilus samples
rm *-CM[1-9]_*.fastq.gz

# Rename sach
for file in *.fastq.gz ; do mv $file ${file//DM0Sach/DM0-Syn} ; done

# Rename all DC to Syn
for file in *.fastq.gz ; do mv $file ${file//-DC_/-Syn_} ; done

# Rename samples in samplesheet
cat SampleSheet.csv | grep -v '\-CM[1-9]' \
| sed 's/DM0Sach/DM0-Syn/g' | sed 's/-DC_/-Syn_/g' \
| sed 's/fwhF2/CB3DR_fwhF2/g' \
| sed 's/BF1/CB3DR_BF1/g' \
| sed 's/SauronS878/CB3DR_SauronS878/g' \
> SampleSheet_CB3DR.csv

####### Renaming for CK3HD
cd /group/pathogens/Alexp/Metabarcoding/dros_metabarcoding/data/CK3HD

# Remove any Carpophilus samples
rm *-CM[1-9]_*.fastq.gz

# Rename sach
for file in *.fastq.gz ; do mv $file ${file//DM0Sach/DM0-Syn} ; done

# Rename all DC to Syn
for file in *.fastq.gz ; do mv $file ${file//-DC_/-Syn_} ; done

# Rename samples in samplesheet
cat SampleSheet.csv | grep -v '\-CM[1-9]' \
| sed 's/DM0Sach/DM0-Syn/g' | sed 's/-DC_/-Syn_/g' \
| sed 's/fwhF2/CK3HD_fwhF2/g' \
| sed 's/BF1/CK3HD_BF1/g' \
| sed 's/SauronS878/CK3HD_SauronS878/g' \
> SampleSheet_CK3HD.csv

####### Renaming for HLVKYDMXX
cd /group/pathogens/Alexp/Metabarcoding/dros_metabarcoding/data/HLVKYDMXX

# Remove any Carpophilus samples
rm *_CM[1-9]-ex*.fastq.gz
rm *_CM[1-9]-s-ex*.fastq.gz
rm *_CM[1-9][0-9]-ex*.fastq.gz
rm *_CM[1-9][0-9]-s-ex*.fastq.gz
rm *_CML[1-9]-ex*.fastq.gz
rm *_CT[1-9]-ex*.fastq.gz
rm *_CT[1-9][0-9]-ex*.fastq.gz
rm *_CT[1-9]-s-ex*.fastq.gz
rm *_CT[1-9]dup*.fastq.gz

# Rename DM
for file in *.fastq.gz ; do mv $file ${file//_DM/_D100M} ; done

# Rename all DC to Syn
for file in *.fastq.gz ; do mv $file ${file//DC-/Syn-} ; done

#Rename individual samples that were mixed up in submitted samplesheet
for file in *.fastq.gz ; do mv $file ${file//D250M1-/D250M4REP-} ; done # Works
for file in *.fastq.gz ; do mv $file ${file//D250M4-/D250M2REP-} ; done # Works
for file in *.fastq.gz ; do mv $file ${file//D250M5-/D250M3REP-} ; done #FAILED library
for file in *.fastq.gz ; do mv $file ${file//D250M3-/D250M1REP-} ; done #FP suzukii - low reads
for file in *.fastq.gz ; do mv $file ${file//D250M2-/D250M5REP-} ; done #Works

for file in *.fastq.gz ; do mv $file ${file//D500M1-/D500M4REP-} ; done #Works
for file in *.fastq.gz ; do mv $file ${file//D500M4-/D500M1REP-} ; done #FP Suzukii
for file in *.fastq.gz ; do mv $file ${file//D500M5-/D500M2REP-} ; done #Works
for file in *.fastq.gz ; do mv $file ${file//D500M3-/D500M3REP-} ; done #Works but low reads for Suz + Biarmipes 
for file in *.fastq.gz ; do mv $file ${file//D500M2-/D500M5REP-} ; done #Works

for file in *.fastq.gz ; do mv $file ${file//D1000M1-/D1000M3REP-} ; done #Works
for file in *.fastq.gz ; do mv $file ${file//D1000M4-/D1000M1REP-} ; done #Works
for file in *.fastq.gz ; do mv $file ${file//D1000M5-/D1000M2REP-} ; done #Works
for file in *.fastq.gz ; do mv $file ${file//D1000M3-/D1000M5REP-} ; done #Works
for file in *.fastq.gz ; do mv $file ${file//D1000M2-/D1000M4REP-} ; done #Works

for file in *.fastq.gz ; do mv $file ${file//REP-/-} ; done


# Rename samples in samplesheet
cat SampleSheet.csv |grep -v '_CM[1-9]-ex\|_CM[1-9]-s-ex\|_CM[1-9][0-9]-ex\|_CM[1-9][0-9]-s-ex\|_CML[1-9]-ex\|_CT[1-9]-ex\|_CT[1-9][0-9]-ex\|_CT[1-9]-s-ex\|_CT[1-9]dup' \
| sed 's/_DM/_D100M/g' | sed 's/DC-/Syn-/g' \
| sed 's/D250M1-/D250M4REP-/g' \
| sed 's/D250M4-/D250M2REP-/g' \
| sed 's/D250M5-/D250M3REP-/g' \
| sed 's/D250M3-/D250M1REP-/g' \
| sed 's/D250M2-/D250M5REP-/g' \
| sed 's/D500M1-/D500M4REP-/g' \
| sed 's/D500M4-/D500M1REP-/g' \
| sed 's/D500M5-/D500M2REP-/g' \
| sed 's/D500M3-/D500M3REP-/g' \
| sed 's/D500M2-/D500M5REP-/g' \
| sed 's/D1000M1-/D1000M3REP-/g' \
| sed 's/D1000M4-/D1000M1REP-/g' \
| sed 's/D1000M5-/D1000M2REP-/g' \
| sed 's/D1000M3-/D1000M5REP-/g' \
| sed 's/D1000M2-/D1000M4REP-/g' \
| sed 's/REP-/-/g' \
| sed 's/Sample_HLVKYDMXX/HLVKYDMXX/g' \
> SampleSheet_HLVKYDMXX.csv

##### Rename mixed samples
#rownames(seqtab)  <- rownames(seqtab) %>%
#  str_replace_all("D250M1-", "D250M4REP-") %>% # Works
#  str_replace_all("D250M4-", "D250M2REP-") %>% # Works
#  str_replace_all("D250M5-", "D250M3REP-") %>% #FAILED library
#  str_replace_all("D250M3-", "D250M1REP-") %>% #FP suzukii - low reads
#  str_replace_all("D250M2-", "D250M5REP-") %>% #Works

#  str_replace_all("D500M1-", "D500M4REP-") %>% #Works 
#  str_replace_all("D500M4-", "D500M1REP-") %>% #FP Suzukii
#  str_replace_all("D500M5-", "D500M2REP-") %>% #Works
#  str_replace_all("D500M3-", "D500M3REP-") %>% #Works but low reads for Suz + Biarmipes 
#  str_replace_all("D500M2-", "D500M5REP-") %>% #Works

#  str_replace_all("D1000M1-", "D1000M3REP-") %>% #Works
#  str_replace_all("D1000M4-", "D1000M1REP-") %>% #Works
#  str_replace_all("D1000M5-", "D1000M2REP-") %>% #Works
#  str_replace_all("D1000M3-", "D1000M5REP-") %>% #Works
#  str_replace_all("D1000M2-", "D1000M4REP-") %>% #Works
#  str_replace_all("CM10-", "CM9REP-") %>%
#  str_replace_all("CM11-", "CM10REP-") %>%
#  str_replace_all("CM9-", "CM11REP-") %>%
#  str_replace_all("CML2-", "CML6REP-")%>%
#  str_replace_all("CML3-", "CML2REP-")%>%
#  str_replace_all("CML4-", "CML3REP-")%>%
#  str_replace_all("CML5-", "CML4REP-")%>%
#  str_replace_all("CML6-", "CML5REP-")%>%
#  str_replace_all("CT5-", "CT4REP-")%>%
#  str_replace_all("CT4-", "CT5REP-") %>%
#  str_replace_all("REP", "")

