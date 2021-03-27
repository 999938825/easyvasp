dir=$(pwd)
for tmp in $(ls $dir)
do
  if [  -d $tmp  ]; then
      cd $dir/$tmp
      cp ../vaspjob.pbs .
      qsub vaspjob.pbs
      cd $dir
  fi
done
