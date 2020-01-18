#!/bin/bash -ex
MIN=1
MAX=950
file=curve

if [ -f $file ] ; then
  rm $file
fi

for ((i = MIN; i <= MAX; i+=1)); do 
  for ((j = 1; j <= 5; j++)); do 
    ./train $i model_init.txt seq_model_0"$j".txt model_0"$j".txt
  done
  ./test modellist.txt testing_data1.txt result1.txt
  python3 calc_acc.py result1.txt testing_answer.txt >> process
done

