#!/bin/bash

echo "Start Running"

cmndpath="/home/Storage/MG5_Study_Group/Cards/author"
cmscardpath="/home/Storage/MG5_Study_Group/Cards/delphes_card_CMS.dat"
savepath="/home/Storage/MG5_Study_Group/ROOT/author"
filename="SVJ"

i=0
while [ $i != 11 ]
do
# for i in 3 #0 3 9
# do
   echo i=$i
   
#    j=10
#    while [ $j != 11 ]
#    do
    for j in 1 5 10 50 100 150 200 250 300 350  #  1 5 10 100 300 #10 100 300 #1 5 10 #
    do
    
    echo j=$j

       date +"%Y %b %m"
       date +"%r"
       echo "s-channel"
#        python /home/Storage/MG5_Study_Group/Cards/create_cmnd_v2.py $i $j
       python /home/Storage/MG5_Study_Group/Cards/create_cmnd_mimic_author.py $i $j

       nohup ./DelphesPythia8 $cmscardpath $cmndpath/"$filename"_"$i"_"$j".cmnd $savepath/"$filename"_"$i"_"$j".root > $savepath/"$filename"_"$i"_"$j".log &

       date +"%Y %b %m"
       date +"%r"
#        j=$(($j+))
     
   done  
   
   
   date +"%Y %b %m"
   date +"%r"
   i=$(($i+1))

done

echo "Finish"

date
