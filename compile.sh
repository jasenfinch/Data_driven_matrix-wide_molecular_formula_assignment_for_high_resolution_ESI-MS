logfile=$(date +"misc/logs/%Y-%m-%d_%H-%M-%S_compile.log")

/usr/lib/R/bin/Rscript 'misc/run.R'  >> $logfile 2>&1