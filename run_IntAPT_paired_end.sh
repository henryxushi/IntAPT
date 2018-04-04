BAMLIST=/home/cbil/Henry/IntAPT_demo/demobamlist.txt
INTAPT_PATH=/home/cbil/Henry/IntAPT_v1.1
OUTPUT=/home/cbil/Henry/IntAPT_demo/demo_paired_end
NUM_OF_PROCESSES=5

chmod 777 $INTAPT_PATH/IntAPT
chmod 777 $INTAPT_PATH/tools/samtools
chmod 777 $INTAPT_PATH/tools/processsamS

#please add the boost library and bamtools library to system lib using the following command
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:

$INTAPT_PATH/IntAPT -b $BAMLIST -o $OUTPUT -r p -p $NUM_OF_PROCESSES --IntAPTpath $INTAPT_PATH/tools
