#!/bin/csh

###################################################################
#
# This script was written by Garrett A. Meek,
#  (garrett.a.meek@gmail.com), in February/March of 2019,
# and reproduces the analyses described in:
#
# Carmody LA, Caverly LJ, Foster BK, Rogers MAM, Kalikin LM, et al. (2018) Fluctuations in airway bacterial communities associated with clinical states and disease stages in cystic fibrosis. PLOS ONE 13(3): e0194060.
# URL: https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0194060
#
# The protocol is divided into the following sections:
#
# 1) User input
# 2) Download project-associated data from NCBI
# 3) Prepare for and run mothur analysis
#    A) Download current copy of mothur analysis standard operating procedure (SOP),
#       which is cited as the analytical protocol followed in Carmody LA, et al.
#    B) Modify mothur analysis script with input for the current project
#    C) Run mothur analysis
# 4) Assemble datasets from mothur output and NCBI metadata
# 5) Reproduce figures 1 and 2
# 6) Reproduce generalized estimating equations analyses
#
###################################################################

###################################################################
#
# 1) User input
#
###################################################################

# Run-time options (used for debugging)
set num_processors = "6"
set download_fastq_files = "False"
set get_new_data = "False"
set make_mothur_input = "False"
set make_fasta_files = "False"
set run_mothur = "True"
set remove_old_output = "False"
set make_new_datasets = "True"
set current_dir = "`pwd`"
set data_dir = "$current_dir/data"
set mothur_input_dir = "$current_dir/mothur_input"
set mothur_output_dir = "$current_dir/mothur_output"
set figures_dir = "$current_dir/figures"
mkdir -p $mothur_input_dir
mkdir -p $mothur_output_dir
mkdir -p $data_dir
mkdir -p $figures_dir
# NCBI database Project ID:
set project_id = "PRJNA423040" 
# Filename for the mothur analysis script
set mothur_batch = "$mothur_input_dir/mothur_steps.batch"
# Download destination for fastq files
set fastq_data_dir = "/mnt/d/data/PLOS_One_Carmody_LiPuma_2018_data/fastq"
mkdir -p $fastq_data_dir
set sra_download_dir = "/home/gmeek/ncbi/public/sra"

###################################################################
#
# 2) Download project-associated data from NCBI
#
###################################################################

if ( "$remove_old_output" == "True") then
 rm $mothur_output_dir/*
endif

# If $get_data = True we start this protocol from scratch, 
# downloading all project-associated files from NCBI
if ( "$get_new_data" == "True") then
# Make an accession list from the Project ID:
# If the accession list already exists, remove it
 if ( -f $data_dir/SRR_Acc_list.txt ) then
  rm $data_dir/SRR_Acc_list.txt
  rm $data_dir/full_dataset.txt
  rm $data_dir/Sample_name_list.txt
  rm $data_dir/stability.files
 endif

 echo "Getting accession list from NCBI"
# Get the full accession list for the project
 esearch -db sra -query $project_id | efetch --format runinfo | cut -d ',' -f 1 | grep SRR > $data_dir/SRR_Acc_List.txt

# Get the sample name list
 esearch -db sra -query $project_id | efetch --format runinfo | cut -d ',' -f 12 | grep SP > $data_dir/Sample_name_list.txt

# Get all metadata for each biosample
 esearch -db sra -query $project_id | efetch --format native > $data_dir/metadata.xml
 esearch -db sra -query $project_id | efetch --format runinfo > $data_dir/full_dataset.txt

 if ( "$download_fastq_files" == "True") then
# Download the SRA files in the accession list:
  prefetch --option-file $data_dir/SRR_Acc_List.txt
# Convert the SRA files to fasta files:
  cat $data_dir/SRR_Acc_List.txt | xargs -I sra_id fastq-dump --fasta $sra_download_dir/sra_id.sra --outdir $fastq_data_dir
 endif

endif

# Now we have fetched and formatted all data for this project
# and are ready to create input files for mothur analysis

###################################################################
#
# 3) Prepare for and run mothur analysis
#
###################################################################

###################################################################
#
#    A) Download current copy of mothur analysis standard operating procedure (SOP),
#       which is cited as the analytical protocol followed in Carmody LA, et al.
#
###################################################################

if ( "$make_mothur_input" == "True") then

 echo "Making input files for mothur analysis."

# Record the date that we are taking a snapshot of the Mothur 454 SOP.
 set date = "`date +'%m-%d-%y'`"

 echo "Downloading a current copy of the mothur 454 sequencing platform standard operating procedure"
 echo "from: https://www.mothur.org/wiki/454_SOP"

# Download an .html copy of the mothur 454 SOP
 wget --html-extension -O $mothur_input_dir/454_SOP_mothur_website_$date.html "https://www.mothur.org/wiki/454_SOP"

# Extract, from the .html file, all mothur command line input, and save to "mothur_steps.batch"
 grep "<pre>mothur &gt;" $mothur_input_dir/454_SOP_mothur_website_$date.html > $mothur_batch

# Reformat the mothur command lines so that the syntax of mothur_steps.batch is mothur-compatible
 sed -i 's/<pre>//;s/&gt;/>/' $mothur_batch
# Make a copy of the original 454_SOP for protocol comparison
 cp $mothur_batch $mothur_input_dir/454_SOP_pure.batch

###################################################################
#
#    B) Modify mothur analysis script with input for the current project
#
###################################################################

# We start our protocol with '.fastq' files, and the sequences for this project
# only contain forward reads.  Thus, we follow the recommended mothur protocol discussed here:
# https://github.com/mothur/mothur/issues/396).
#
# Notably, this protocol is different from the 454 SOP that Carmody, et al., cite in Supplementary Information, since the SOP requires paired reads,
# and the sequences for this project are un-paired.
#
# Consequently, we modify the 454 SOP steps during the stage where sequences are read, trimmed, and evaluated for their quality,
# including skipping the creation of 'stability.files', which is only used to make contigs with paired reads.

 echo "Modifying the mothur analysis batch script"

##############################################################
#
# Remove all steps until .fastq files are converted to .fasta and .qual files
#
##############################################################

 sed -i '/mothur > sffinfo(sff=GQY1XT001.sff, flow=T)/d' $mothur_batch
 sed -i '/mothur > summary.seqs(fasta=GQY1XT001.fasta)/d' $mothur_batch 
 sed -i '/mothur > trim.flows(flow=GQY1XT001.flow, oligos=GQY1XT001.oligos, pdiffs=2, bdiffs=1, processors=2)/d' $mothur_batch
 sed -i '/mothur > shhh.flows(file=GQY1XT001.flow.files, processors=2)/d' $mothur_batch
 sed -i '/mothur > trim.seqs(fasta=GQY1XT001.shhh.fasta, name=GQY1XT001.shhh.names, oligos=GQY1XT001.oligos, pdiffs=2, bdiffs=1, maxhomop=8, minlength=200, flip=T, processors=2)/d' $mothur_batch
 sed -i '/mothur > summary.seqs(fasta=GQY1XT001.shhh.trim.fasta, name=GQY1XT001.shhh.trim.names)/d' $mothur_batch
 sed -i '/mothur > summary.seqs/d' $mothur_batch
 sed -i '/mothur > get.current/d' $mothur_batch
 sed -i '/mothur > trim.seqs(fasta=GQY1XT001.fasta, oligos=GQY1XT001.oligos, qfile=GQY1XT001.qual, maxambig=0, maxhomop=8, flip=T, bdiffs=1, pdiffs=2, qwindowaverage=35, qwindowsize=50, processors=2)/d' $mothur_batch
 sed -i '/mothur > trim.seqs(fasta=GQY1XT001.fasta, oligos=GQY1XT001.oligos, maxambig=0, maxhomop=8, flip=T, bdiffs=1, pdiffs=2, keepfirst=200, processors=2)/d' $mothur_batch

##############################################################
#
# We skip sequence analyses that aren't part of the PLOS One paper
#
##############################################################

 sed -i '/mothur > clearcut(phylip=final.phylip.dist)/d' $mothur_batch
 sed -i '/mothur > collect.single(shared=final.an.shared, calc=chao-invsimpson, freq=100)/d' $mothur_batch
 sed -i '/mothur > rarefaction.single(shared=final.an.shared, calc=sobs, freq=100)/d' $mothur_batch
 sed -i '/mothur > summary.single(calc=nseqs-coverage-sobs-invsimpson, subsample=4419)/d' $mothur_batch
 sed -i '/mothur > heatmap.bin(shared=final.an.shared, scale=log2, numotu=50)/d' $mothur_batch
 sed -i '/mothur > heatmap.sim(calc=jclass-thetayc)/d' $mothur_batch
 sed -i '/mothur > venn(groups=F003D000-F003D002-F003D004-F003D006)/d' $mothur_batch
 sed -i '/mothur > summary.shared(calc=sharedchao, groups=F003D000-F003D002-F003D004-F003D006, all=T)/d' $mothur_batch
 sed -i '/mothur > tree.shared(calc=thetayc-jclass, subsample=4419)/d' $mothur_batch
 sed -i '/mothur > parsimony(tree=final.an.thetayc.0.03.ave.tre, group=mouse.sex_time.design, groups=all)/d' $mothur_batch
 sed -i '/mothur > unifrac.weighted(tree=final.an.thetayc.0.03.ave.tre, group=mouse.sex_time.design, random=T)/d' $mothur_batch
 sed -i '/mothur > unifrac.unweighted(tree=final.an.thetayc.0.03.ave.tre, group=mouse.sex_time.design, random=T, groups=all)/d' $mothur_batch
# sed -i '/mothur > dist.shared(shared=final.an.shared, calc=thetayc-jclass, subsample=4419)/d' $mothur_batch
 sed -i '/mothur > pcoa(phylip=final.an.thetayc.0.03.lt.ave.dist)/d' $mothur_batch
 sed -i '/mothur > nmds(phylip=final.an.thetayc.0.03.lt.ave.dist)/d' $mothur_batch
 sed -i '/mothur > nmds(phylip=final.an.thetayc.0.03.lt.ave.dist, mindim=3, maxdim=3)/d' $mothur_batch
 sed -i '/mothur > amova(phylip=final.an.thetayc.0.03.lt.ave.dist, design=mouse.sex_time.design)/d' $mothur_batch

# Replace SOP filenames with filenames specific to the current project
 sed -i 's/fasta=.*,/fasta=current,/g' $mothur_batch
 sed -i 's/file=.*,/file=current,/g' $mothur_batch
 sed -i 's/shared=.*,/shared=current,/g' $mothur_batch
 sed -i 's/name=.*,/name=current,/g' $mothur_batch
 sed -i 's/name=.*)/name=current)/g' $mothur_batch
 sed -i 's/accnos=.*,/accnos=current,/g' $mothur_batch
 sed -i 's/processors=2/processors=6/' $mothur_batch
 sed -i 's/mothur > //' $mothur_batch

# We use a renamed version of the Silva V4 fasta database, silva.v4.fasta
# for reference-based sequencing
 sed -i 's/silva.bacteria.fasta/silva.v4.fasta/' $mothur_batch

# The analyses required to reproduce the results of the PLOS One paper include
# the following additional modifications to teh beginning of the 454 SOP:

##############################################################
#
# Process '.fastq' files for sequence analysis
#
##############################################################

# sed -i "1iset.dir(input=$mothur_input_dir)" $mothur_batch
 sed -i "1iscreen.seqs(fasta=current, group=current, maxambig=0, maxlength=275,processors=$num_processors)" $mothur_batch # Screen based on sequence length
 sed -i '2icount.seqs(name=current, group=current)' $mothur_batch # Count the number of sequences in each unique group

 if ( "$make_fasta_files" == "True") then
  set temp_batch = "$mothur_input_dir/mothur_batch_temp"
  touch $temp_batch
  set i = 1
  set num_samples = `wc -l < $data_dir/Sample_name_list.txt`
  while ($i <= $num_samples)
# Get the sample names
   set sample_name = "`sed -n '$i""p' $data_dir/Sample_name_list.txt`"
# Convert '.fastq' files to '.fasta' and '.qual' files:
   echo "fastq.info(fastq=$fastq_data_dir/$sample_name.fastq)" >> $temp_batch # Deconstruct the .fastq file
   @ i = $i + 1
  end

  sed -i "1iset.dir(input=$fastq_data_dir)" $temp_batch
  sed -i "2iset.dir(output=$mothur_output_dir)" $temp_batch
  echo "Converting fastq files to fasta files"
  mothur $temp_batch > $mothur_output_dir/mothur.out

  cat $temp_batch $mothur_batch > $mothur_output_dir/cat_batch
  mv $mothur_output_dir/cat_batch $mothur_batch
  rm $temp_batch

  set temp_batch = "$mothur_input_dir/mothur_batch_temp"
  touch $temp_batch
  set i = 1
  set num_samples = `wc -l < $data_dir/Sample_name_list.txt`
  while ($i <= $num_samples)
# Get the sample names
   set sample_name = "`sed -n '$i""p' $data_dir/Sample_name_list.txt`"
   if ( -f "$mothur_output_dir/$sample_name.fasta") then
    set fasta_file = "$mothur_output_dir/$sample_name.fasta"
    set qual_file = "$mothur_output_dir/$sample_name.qual"
   endif
   if ( -f "$mothur_output_dir/$sample_name.trim.fasta") then
    set fasta_file = "$mothur_output_dir/$sample_name.trim.fasta"
    set qual_file = "$mothur_output_dir/$sample_name.trim.qual"
   endif
   echo "trim.seqs(fasta=$mothur_output_dir/$fasta_file,qfile=$mothur_output_dir/$qual_file,processors=$num_processors)" >> $temp_batch # Trim based on sequence quality
   @ i = $i + 1
  end

  sed -i "1iset.dir(input=$fastq_data_dir)" $temp_batch
  sed -i "2iset.dir(output=$mothur_output_dir)" $temp_batch
  echo "Trimming sequences..."
  mothur $temp_batch >> $mothur_output_dir/mothur.out

  cat $temp_batch $mothur_batch > $mothur_output_dir/cat_batch
  mv $mothur_output_dir/cat_batch $mothur_batch
  rm $temp_batch
 endif

 set temp_batch = "$mothur_input_dir/mothur_batch_temp"
 rm $temp_batch
 touch $temp_batch

###
#
# Merge fasta files
#
###

 set i = 1
 set sample_list = ""
 set max_list_length = 10
 set current_list_length = 0
 set merge_file_index = 1
 set num_samples = `wc -l < $data_dir/Sample_name_list.txt`
 while ($i <= $num_samples)
# Get the sample names
  set sample_name = "`sed -n '$i""p' $data_dir/Sample_name_list.txt`"
# Create a concatenated list of the sample names in order to trim and merge them
  if ( -f "$mothur_output_dir/$sample_name.fasta") then
   set fasta_file = "$mothur_output_dir/$sample_name.fasta"
  endif
  if ( -f "$mothur_output_dir/$sample_name.trim.fasta") then
   set fasta_file = "$mothur_output_dir/$sample_name.trim.fasta"
  endif
  if ("$sample_list" != "") then
   set sample_list = "$sample_list-$fasta_file"
  else
   set sample_list = "$fasta_file"
  endif
# Convert '.fastq' files to '.fasta' and '.qual' files:
# Trim based on sequence quality
  if ($current_list_length == $max_list_length) then
   echo "merge.files(input=$sample_list, output=merge_$merge_file_index.fasta)" >> $temp_batch
   set current_list_length = 0
   set sample_list = ""
   @ merge_file_index = $merge_file_index + 1
  else
   @ current_list_length = $current_list_length + 1
  endif
  @ i = $i + 1
 end

###
#
# Merge fasta files again
#
###
 set i = 1
 set sample_list = ""
 set max_list_length = 10
 set current_list_length = 0
 set total_merge_files = `expr $merge_file_index - 1`
 while ($i <= $total_merge_files)
  set fasta_file = "$mothur_output_dir/merge_$i.fasta"
  if ("$sample_list" != "") then
   set sample_list = "$sample_list-$fasta_file"
  else
   set sample_list = "$fasta_file"
  endif
  if ($current_list_length == $max_list_length) then
   echo "merge.files(input=$sample_list, output=merge_$merge_file_index.fasta)" >> $temp_batch
   set current_list_length = 0
   set sample_list = ""
   @ merge_file_index = $merge_file_index + 1
  else
   @ current_list_length = $current_list_length + 1
  endif
  @ i = $i + 1
 end
 echo "merge.files(input=$sample_list, output=merge_$merge_file_index.fasta)" >> $temp_batch

###
#
# Consolidate all fasta files
#
###

 set sample_list = ""
 set merge_files = `expr $merge_file_index`
 set merge_file_index = `expr $total_merge_files + 1`
 echo "$merge_files"
 echo "$merge_file_index"
 while ($merge_file_index <= $merge_files)
  set fasta_file = "$mothur_output_dir/merge_$merge_file_index.fasta"
  if ("$sample_list" != "") then
   set sample_list = "$sample_list-$fasta_file"
  else
   set sample_list = "$fasta_file"
  endif
  @ merge_file_index = $merge_file_index + 1
 end 

 echo "merge.files(input=$sample_list, output=full.fasta)" >> $temp_batch
 sed -i "1iset.dir(input=$fastq_data_dir)" $temp_batch
 sed -i "2iset.dir(output=$mothur_output_dir)" $temp_batch

 echo "Merging fasta files"
 mothur $temp_batch >> $mothur_output_dir/mothur.out

 cat $temp_batch $mothur_batch > $mothur_output_dir/cat_batch
 mv $mothur_output_dir/cat_batch $mothur_batch
 rm $temp_batch

 sed -i "1iset.dir(input=$data_dir)" $mothur_batch
 sed -i "2iset.dir(output=$mothur_output_dir)" $mothur_batch

endif

###################################################################
#
#    C) Run mothur analysis
#
###################################################################

if ("$run_mothur" == "True") then
 echo "Running mothur"
 mothur $mothur_batch >> $mothur_output_dir/mothur.out
endif

###################################################################
#
# 4) Assemble datasets from mothur output and biosample metadata
#
###################################################################

# Assemble datasets for each figure:
set raw_metadata = "$data_dir/metadata.xml"
set fev1_data = "$data_dir/fev1_data.txt"
set age_data = "$data_dir/age_data.txt"
set clinical_state_data = "$data_dir/clinical_state_data.txt"
touch $fev1_data
touch $age_data
touch $clinical_state_data

set i = 1
set num_samples = `wc -l < $data_dir/Sample_name_list.txt`
while ($i <= $num_samples)
# Get the sample names
#   echo "$i"
   set sample_name = "`sed -n '$i""p' $data_dir/Sample_name_list.txt`"
# Get the accession names
   set accession_id = "`sed -n '$i""p' $data_dir/SRR_Acc_List.txt`"
#   echo "grep sample_name="'"'$sample_name'"'" $raw_metadata | grep accession="'"'$accession_id'"'
   grep sample_name='"'$sample_name'"' $raw_metadata | grep accession='"'$accession_id'"' | grep -o -P '(?<=host_fev1</TAG><VALUE>).*(?=</VALUE></SAMPLE_ATTRIBUTE><SAMPLE_ATTRIBUTE><TAG>host_disease_aggressiveness)' >> $fev1_data
   grep sample_name='"'$sample_name'"' $raw_metadata | grep accession='"'$accession_id'"' | grep -o -P '(?<=host_age</TAG><VALUE>).*(?=</VALUE></SAMPLE_ATTRIBUTE><SAMPLE_ATTRIBUTE><TAG>host_fev1)' >> $age_data
   grep sample_name='"'$sample_name'"' $raw_metadata | grep accession='"'$accession_id'"' | grep -o -P '(?<=host_age</TAG><VALUE>).*(?=</VALUE></SAMPLE_ATTRIBUTE><SAMPLE_ATTRIBUTE><TAG>host_fev1)' >> $clinical_state_data
   @ i = $i + 1
end

# Combine counts from anaerobic genera:
# Actinomyces, Fusobacterium, Gemella, Granulicatella,
# Porphyromonas, Prevotella, Rothia, Streptococcus, Veillonella
set anaerobe_data = "$data_dir/anaerobe_data.txt"

# Figure 1-A
#
# Contains the relative anaerobe abundance for all samples
# when grouped by clinical state (B,E,T,R) 
set fig_1_A_data = "$data_dir/fig_1_A_data.txt"
set fig_1_B_data = "$data_dir/fig_1_B_data.txt"
set fig_1_C_data = "$data_dir/fig_1_C_data.txt"
set fig_1_D_data = "$data_dir/fig_1_D_data.txt"
set fig_2_A_1_data = "$data_dir/fig_2_A_1_data.txt"
set fig_2_A_2_data = "$data_dir/fig_2_A_2_data.txt"
set fig_2_A_3_data = "$data_dir/fig_2_A_3_data.txt"
set fig_2_B_1_data = "$data_dir/fig_2_B_1_data.txt"
set fig_2_B_2_data = "$data_dir/fig_2_B_2_data.txt"
set fig_2_B_3_data = "$data_dir/fig_2_B_3_data.txt"

###################################################################
#
# 5) Reproduce figures 1 and 2
#
###################################################################

###################################################################
#
# 6) Reproduce generalized estimating equations results
#
###################################################################
