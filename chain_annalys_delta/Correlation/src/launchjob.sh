#!/bin/csh
foreach i (`seq 0 9`)
	foreach j (`seq 0 9`)
		echo $i $j
		clubatch python /home/gpfs/manip/mnt0607/bao/hdumasde/Code/CrossCorrelation/Python/Correlation/main.py $i $j
		sleep 0.1
	end
	sleep 60
end
