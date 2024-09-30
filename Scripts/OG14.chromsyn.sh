#!/bin/bash

####################################### ::: SETUP ::: ############################################
date=$(date +%y%m%d)
#---------------
module add singularity/4.1.0-nompi
SLIMSUITE=/software/projects/pawsey0812/rjedwards/slimsuite
THREADS=16
echo "THREADS: $THREADS"
#---------------
# Run code

ABASE=OG14.haps
HAP1=OG14_v230209.hic1.3.curated.hap1.chr_level.fa
HAP2=OG14_v230209.hic1.3.curated.hap2.chr_level.fa
PREF1=OG14.hap1; PREF2=OG14.hap2

#i# Select a BUSCO Database and path to the Compleasm datasets
LPATH=/software/projects/pawsey0812/rjedwards/compleasmdb
ODB10=actinopterygii_odb10

####################################### ::: INTRO ::: #############################################
#i# Input is the two haplotypes, named $ABASE.hap[12].fasta by default
#i# Option $PREF1 and $PREF2 can set the ChromSyn genome labels for the two haplotypes (else $HAPBASE)
HAPBASE1=$(basename ${HAP1%.*})
HAPBASE2=$(basename ${HAP2%.*})
#i# Output is named with the prefix $NEWBASE
if [ -z "$NEWBASE" ]; then
	NEWBASE=$(basename $ABASE)
fi

####################################### ::: MAIN CODE ::: #############################################

### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
### ~~~~ STEP 1: Telociraptor generation of input files. ~~~~ ###
## Goal: Telomere prediction for each haplotype.
## Key inputs: Genome Assembly haplotypes.
## Key outputs: Telomeres and gaps files *.telomeres.tdt *.gaps.tdt
### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
STEP=Telociraptor
for GENOME in $HAP1 $HAP2; do
	GENBASE=$(basename ${GENOME%.*})
	echo "$GENOME -> $GENBASE.* ..."
	KEYOUT=$GENBASE.telomeres.tdt
	if [ ! -f "$KEYOUT" ]; then
		python3 $SLIMSUITE/tools/telociraptor.py -seqin $GENOME -basefile $GENBASE tweak=F chromsyn=T i=-1 v=0
	else
		echo "[$(date)] File found: $KEYOUT - skipping $STEP step"
	fi
done
touch $STEP.done

### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
### ~~~~ STEP 2: TIDK Telomere Repeat prediction. ~~~~ ###
## Goal: Generate additional locations of telomeric repeats.
## Key inputs: Genome Assembly haplotypes.
## Key outputs: TIDK repeat windows (*.tidk.tsv)
### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
STEP=TIDK
for GENOME in $HAP1 $HAP2; do
	GENBASE=$(basename ${GENOME%.*})
	echo "$GENOME -> $GENBASE.* ..."
	KEYOUT=$GENBASE.tidk.tsv
	if [ ! -f "$KEYOUT" ]; then
		singularity run $SING/tidk:0.2.31.sif tidk search -o $GENBASE -s AACCCT -d ./ -e tsv $GENOME
		awk '$3 > 9 || $4 > 9' ${GENBASE}_telomeric_repeat_windows.tsv > $GENBASE.tidk.tsv
	else
		echo "[$(date)] File found: $KEYOUT - skipping $STEP step"
	fi
done
touch $STEP.done

### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
### ~~~~ STEP 3: Run Compleasm on two haplotypes. ~~~~ ###
## Goal: Generate BUSCO Complete and Duplicated sets for synteny and plotting
## Key inputs: Genome Assembly haplotypes.
## Key outputs: BUSCO-format full tables (*.full_table.tsv)
### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
STEP=Compleasm
for GENOME in $HAP1 $HAP2; do
	GENBASE=$(basename ${GENOME%.*})
	echo "$GENOME -> $GENBASE.* ..."
	KEYOUT=$GENBASE.$ODB10.full_table.tsv
	if [ ! -f "$KEYOUT" ]; then
		# Run Compleasm
		RUN=run_$GENBASE.compleasm
		echo singularity run $SING/compleasm:v0.2.5.sif compleasm run -a $GENOME -o $RUN -t $THREADS -l ${ODB10/_odb10/} -L $LPATH
		singularity run $SING/compleasm:v0.2.5.sif compleasm run -a $GENOME -o $RUN -t $THREADS -l ${ODB10/_odb10/} -L $LPATH
		# Tidy Compleasm 
		echo "Cleanup $ODB Compleasm run $RUN"
		CBASE=$GENBASE.$ODB10
		cp -v $RUN/summary.txt $CBASE.compleasm.txt
		cp -v $RUN/$ODB10/full_table.tsv $CBASE.compleasm.tsv
		cp -v $RUN/$ODB10/full_table_busco_format.tsv $CBASE.full_table.tsv
		tar -c $RUN | pigz > $RUN.tar.gz && rm -rf $RUN
	else
		echo "[$(date)] File found: $KEYOUT - skipping $STEP step"
	fi
done
touch $STEP.done

### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
### ~~~~ STEP 4: barrnap rRNA predictions for ChromSyn plots and rRNA contigs   ~~~~ ###
## Goal: Generate rRNA predictions to identify rDNA repeats.
## Key inputs: Genome Assembly haplotypes.
## Key outputs: full-length rRNA predictions reformatted as a features CSV.
### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
STEP=Barrnap
for GENOME in $HAP1 $HAP2; do
	GENBASE=$(basename ${GENOME%.*})
	echo "$GENOME -> $GENBASE.* ..."
	KEYOUT=$GENBASE.rrna.csv
	if [ ! -f "$KEYOUT" ]; then
		sed 's/^[ACGT][ACGT][ACGT][ACGT]/GATC/' $GENOME > $GENBASE.tmp.fasta
		singularity run $SING/barrnap:0.9.sif barrnap --kingdom euk -o $GENBASE.rrna.fa --threads $THREADS < $GENBASE.tmp.fasta | tee $GENBASE.rrna.gff
		rm -v $GENBASE.tmp.fasta
		#i# Make a GFF for plotting with ChromSyn and analysis with DepthKopy
		grep -v 5S $GENBASE.rrna.gff | grep -v partial | tee $GENBASE.rdna.gff
		#i# Find full-length rDNA repeats
		touch $GENBASE.full.ctg.txt
		for S in 5.8 18 28; do
		  echo $GENBASE ${S}S "=>" $(grep "${S}S" $GENBASE.rrna.gff | grep -v partial | awk '{print $1;}' | sort | uniq -c ) | tee -a $GENBASE.full.ctg.txt
		done
		#i# Make a feature file for ChromSyn
		echo SeqName,Start,End,Strand,Col,Shape | tee $GENBASE.rrna.csv
		grep RNA $GENBASE.rrna.gff | grep -v partial | sed 's/=/ /g' | awk '{print $1,$4,$5,$7,$11;}' | sed "s/ /,/g" | sed 's/5S/"white",23/' | sed 's/5.8S/"pink",22/' | sed 's/18S/"skyblue",24/' | sed 's/28S/"seagreen",25/' | tee -a $GENBASE.rrna.csv
	else
		echo "[$(date)] File found: $KEYOUT - skipping $STEP step"
	fi
done
touch $STEP.done

### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
### ~~~~ STEP 5: Simple ChromSyn Compleasm synteny plot. ~~~~ ###
## Goal: Generate a simple ChromSyn Plot of the two haplotypes.
## Key inputs: key outputs from previous steps
## Key outputs: ChromSyn plot of two haplotypes
### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
STEP=ChromSynA
CHROMSYNSET=pdfwidth=20 orphans=FALSE minregion=0 minbusco=1 ticks=1e7
SEQORDER=SUPER_1,SUPER_2,SUPER_3,SUPER_4,SUPER_5,SUPER_6,SUPER_7,SUPER_8
SYNBASE=A1-8; RESTRICT=$SEQORDER

KEYOUT=$NEWBASE.$SYNBASE.pdf
if [ ! -f "$KEYOUT" ]; then

	# Generate the FOFN files from a template:
	echo $PREF1 $HAPBASE1 | tee template.fofn
	echo $PREF2 $HAPBASE2 | tee -a template.fofn

	for FOFNDATA in telomeres.tdt-sequences $ODB10.full_table.tsv-busco gaps.tdt-gaps tidk.tsv-tidk rrna.csv-ft; do
		EXT=$(echo $FOFNDATA | awk -F '-' '{print $1;}')
		FOFN=$(echo $FOFNDATA | awk -F '-' '{print $2;}').fofn
		ls -lh $HAPBASE1.$EXT $HAPBASE2.$EXT && sed "s/$/.$EXT/" template.fofn | tee $FOFN
	done

	for FOFN in sequences.fofn busco.fofn; do
		if [ ! -f "$FOFN" ]; then
			echo "File not found: $FOFN"
			exit 1
		fi
	done

	# Generate the ChromSyn plot, using hap1 as the focus.
	FOCUS=$PREF1
	SETTINGS="$CHROMSYNSET focus=$FOCUS"
	RUN=$NEWBASE
	BASE=$SYNBASE
	singularity run $SING/depthsizer:v1.9.0.sif Rscript $SLIMSUITE/libraries/r/chromsyn.R $SETTINGS basefile=$RUN.$BASE | tee $RUN.$BASE.log
	touch $STEP.done

else
	echo "[$(date)] File found: $KEYOUT - skipping $STEP step" | tee -a $LOG
fi

# Check for key output file
if [ ! -f "$KEYOUT" ]; then
	echo "[$(date)] Failed to generate $KEYOUT" | tee -a $LOG
	exit 1
fi
### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
STEP=ChromSynB
CHROMSYNSET=pdfwidth=20 orphans=FALSE minregion=0 minbusco=1 ticks=1e7
SEQORDER=SUPER_9,SUPER_10,SUPER_11,SUPER_12,SUPER_13,SUPER_14,SUPER_15,SUPER_16
SYNBASE=B9-16; RESTRICT=$SEQORDER

KEYOUT=$NEWBASE.$SYNBASE.pdf
if [ ! -f "$KEYOUT" ]; then

	# Generate the ChromSyn plot, using hap1 as the focus.
	FOCUS=$PREF1
	SETTINGS="$CHROMSYNSET focus=$FOCUS"
	RUN=$NEWBASE
	BASE=$SYNBASE
	singularity run $SING/depthsizer:v1.9.0.sif Rscript $SLIMSUITE/libraries/r/chromsyn.R $SETTINGS basefile=$RUN.$BASE | tee $RUN.$BASE.log
	touch $STEP.done

else
	echo "[$(date)] File found: $KEYOUT - skipping $STEP step" | tee -a $LOG
fi

# Check for key output file
if [ ! -f "$KEYOUT" ]; then
	echo "[$(date)] Failed to generate $KEYOUT" | tee -a $LOG
	exit 1
fi
### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
STEP=ChromSynC
CHROMSYNSET=pdfwidth=20 orphans=FALSE minregion=0 minbusco=1 ticks=1e7
SEQORDER=SUPER_17,SUPER_18,SUPER_19,SUPER_20,SUPER_21,SUPER_22,SUPER_23,SUPER_24
SYNBASE=C17-24; RESTRICT=$SEQORDER

KEYOUT=$NEWBASE.$SYNBASE.pdf
if [ ! -f "$KEYOUT" ]; then

	# Generate the ChromSyn plot, using hap1 as the focus.
	FOCUS=$PREF1
	SETTINGS="$CHROMSYNSET focus=$FOCUS"
	RUN=$NEWBASE
	BASE=$SYNBASE
	singularity run $SING/depthsizer:v1.9.0.sif Rscript $SLIMSUITE/libraries/r/chromsyn.R $SETTINGS basefile=$RUN.$BASE | tee $RUN.$BASE.log
	touch $STEP.done

else
	echo "[$(date)] File found: $KEYOUT - skipping $STEP step" | tee -a $LOG
fi

# Check for key output file
if [ ! -f "$KEYOUT" ]; then
	echo "[$(date)] Failed to generate $KEYOUT" | tee -a $LOG
	exit 1
fi

####################################### ::: FINISH ::: #############################################
echo "#---------------" | tee -a $LOG
echo "[$(date)]: Run complete" | tee -a $LOG
echo "#---------------" | tee -a $LOG

exit 0
