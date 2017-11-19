#!/bin/bash
#
# pwa_omega.bash - bash script for pwa_omega fit jobs
#
# author: richard.t.jones at uconn.edu
# version: november 13, 2016

if [[ ! $1 -gt -1 ]]; then
    echo "Usage: pwa_omega.bash <PROCESS_ID>"
    exit 1
elif [[ $1 -lt 30 ]]; then
    sample=$1
elif [[ $1 -lt 45 ]]; then
    sample=`expr $1 - 30 + 100`
elif [[ $1 -lt 60 ]]; then
    sample=`expr $1 - 45 + 200`
elif [[ $1 -lt 75 ]]; then
    sample=`expr $1 - 60 + 300`
elif [[ $1 -lt 90 ]]; then
    sample=`expr $1 - 75 + 400`
else
    echo "Usage: pwa_omega.bash <PROCESS_ID>"
    exit 1
fi

. setup.sh

cat <<EOI >jobscript.py
#!/bin/env python
import PWAomega
PWAomega.init(1,PWAomega.tbins[$sample],PWAomega.Ebins[$sample])
PWAomega.tool.set_weight_2dhist(PWAomega.gen_weight_2dhist())
PWAomega.multifit(150)
EOI

chmod +x jobscript.py
scl enable python27 ./jobscript.py
ret=$?

rm jobscript.py
exit $ret
