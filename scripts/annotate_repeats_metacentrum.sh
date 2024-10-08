#!/bin/bash
#PBS -N annotate_repeats
#PBS -l select=1:ncpus=50:mem=512gb:scratch_local=300gb
#PBS -l walltime=172:00:00
#PBS -j oe
#PBS -m bae

## Use approximately 1 cpu per 8-10 GB of RAM

# CONFIGURATION - use absolute paths
SINGULARITY_IMAGE="todo.sif"
GENOME="TODO.fasta"
CUSTOM_DATABASE_TAREAN=""  # optional, keep empty if not used
CUSTOM_DATABASE_REPEATS="" # optional, keep empty if not used
OUTPUT_DIR="TODO"
REPEATMASKER_SENSITIVITY="default"  # default, sensitive, quick
REDUCE_LIBRARY="True"  # True, False, if not defined, default is True
# END OF CONFIGURATION



SINGULARITY_IMAGE_BASE_NAME=$(basename $SINGULARITY_IMAGE)
# copy data and singularity image to scratch


rsync -avt $SINGULARITY_IMAGE $SCRATCHDIR/
rsync -avt $GENOME ${SCRATCHDIR}/genome.fasta


# copy is not empty
if [ -n "$CUSTOM_DATABASE_TAREAN" ]; then
     rsync -avt $CUSTOM_DATABASE_TAREAN ${SCRATCHDIR}/custom_database_tarean.fasta

fi

if [ -n "$CUSTOM_DATABASE_REPEATS" ]; then
      rsync -avt $CUSTOM_DATABASE_REPEATS ${SCRATCHDIR}/custom_database_repeats.fasta
fi

cd $SCRATCHDIR

# Generate config.yaml
# format is:
# genome_fasta: genome.fasta
# output_dir: output_dir
# custom_library: path/to/custom_library.fasta
# tandem_repeat_library: path/to/tandem_repeat_library.fasta

echo "genome_fasta: genome.fasta" > ${SCRATCHDIR}/config.yaml
echo "output_dir: output" >> ${SCRATCHDIR}/config.yaml
# if custom database is not empty
if [ -n "$CUSTOM_DATABASE_REPEATS" ]; then
    echo "custom_library: custom_database_repeats.fasta" >> ${SCRATCHDIR}/config.yaml
fi
# if custom database is not empty
if [ -n "$CUSTOM_DATABASE_TAREAN" ]; then
    echo "tandem_repeat_library: custom_database_tarean.fasta" >> ${SCRATCHDIR}/config.yaml
fi

if [ -n "$REPEATMASKER_SENSITIVITY" ]; then
    echo "repeatmasker_sensitivity: $REPEATMASKER_SENSITIVITY" >> ${SCRATCHDIR}/config.yaml
fi

if [ -n "$REDUCE_LIBRARY" ]; then
    echo "reduce_library: $REDUCE_LIBRARY" >> ${SCRATCHDIR}/config.yaml
fi

tmpdir=$SCRATCHDIR/tmp
mkdir -p $tmpdir
export TMPDIR=$tmpdir

# run the singularity container, specify temp directory variable

singularity run -B $SCRATCHDIR --env TMPDIR=$tmpdir $SINGULARITY_IMAGE_BASE_NAME  -c config.yaml -t $PBS_NCPUS


# copy results back to the $OUTPUT_DIR
mkdir -p $OUTPUT_DIR
rsync -avt $SCRATCHDIR/output/ $OUTPUT_DIR/

# copy also databases and genome to output
rsync -avt $SCRATCHDIR/genome.fasta $OUTPUT_DIR/
if [ -n "$CUSTOM_DATABASE_TAREAN" ]; then
    rsync -avt $SCRATCHDIR/custom_database_tarean.fasta $OUTPUT_DIR/
fi
if [ -n "$CUSTOM_DATABASE_REPEATS" ]; then
    rsync -avt $SCRATCHDIR/custom_database_repeats.fasta $OUTPUT_DIR/
fi

rsync -avt $SCRATCHDIR/config.yaml $OUTPUT_DIR/

# Copy the PBS script itself to the OUTPUT_DIR
rsync -avt $0 $OUTPUT_DIR/pbs_script.sh

# Export environment variables to the output directory
env > $OUTPUT_DIR/env.sh

cp -r .snakemake $OUTPUT_DIR/

# clean up scratch TODO after testing


