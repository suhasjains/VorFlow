 for file in 0.*; do awk ' BEGIN{sum = 0; n=0} { if (FNR!=1) { if ($2 > 1./3. - 0.6 && $2 < 1./3. - 0.4) { sum = sum + $4*$4; n = n + 1;}}} END{print FILENAME, sqrt(sum/n)} '  "$file"; done 
