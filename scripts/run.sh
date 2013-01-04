#!/bin/tcsh
foreach f ($*)
	echo "Sumbitting $f to TAMNUN..."
	qsub -q all_l_p $f
end

