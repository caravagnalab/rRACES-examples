#!/bin/bash

if [ "$#" != "2" ]; then
   echo "Syntax error: $0 <PARTITION_NAME> <USER_NAME>"
   exit 1
fi

cat >local_scratch_delete.sh <<EOC
#!/bin/bash

rm -rf /local_scratch/*

exit 0
EOC

for hostname in $(scontrol show hostname $(sinfo -o "%N" -p ${1}| grep -v NODELIST)); do
  sbatch -A ${2} -p ${1} --nodelist=${hostname} local_scratch_delete.sh
done

rm -f local_scratch_delete.sh

exit 0
