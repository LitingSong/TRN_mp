#!bin/bash

#Rscript gtrd_meta_tissue_cell.R 

#rm ./macs2_fdr/brain.bed
#rm ./macs2_fdr/lung.bed
#rm ./macs2_fdr/liver.bed
rm ./macs2_fdr/colon.bed
#rm ./macs2_fdr/stomach.bed

######################
while read Line
do
#echo $Line
#ls ./macs2_fdr/${Line}*20.bed|wc -l
cat ./macs2_fdr/${Line}*20.bed  >>./macs2_fdr/hippocampus.bed

done <hippocampus.txt

######################
while read Line
do
#echo $Line
#ls ./macs2_fdr/${Line}*20.bed|wc -l
cat ./macs2_fdr/${Line}*20.bed  >>./macs2_fdr/intestine.bed

done <intestine.txt


######################
while read Line
do
#echo $Line
cat ./macs2_fdr/${Line}*20.bed  >>./macs2_fdr/colon.bed

done <colon.txt

# ######################
# while read Line
# do
# #echo $Line
# cat ./macs2_fdr/${Line}*20.bed  >>./macs2_fdr/lung.bed
# 
# done <lung.txt
# 
# #######################
# while read Line
# do
# #echo $Line
# cat ./macs2_fdr/${Line}*20.bed  >>./macs2_fdr/liver.bed
# 
# done <liver.txt
# 
# 
# #####################################
# 
# 
# while read Line
# do
# #echo $Line
# cat ./macs2_fdr/${Line}*20.bed  >>./macs2_fdr/stomach.bed
# 
# done <stomach.txt
# 
# ######################
# while read Line
# do
# #echo $Line
# #ls ./macs2_fdr/${Line}*20.bed|wc -l
# cat ./macs2_fdr/${Line}*20.bed  >>./macs2_fdr/brain.bed
# 
# done <brain.txt

