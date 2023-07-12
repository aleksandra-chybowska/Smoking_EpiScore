## Run in Terminal 

cd /Cluster_Filespace/Marioni_Group/Ola/Smoking/BayesR/Outputs/Methylation/

## Set up loop for any output file with a *.csv extension 

for i in *.csv
do

## Extract columns with sigma name - sigmaG (variance explained) and sigmaE (error term - unexplained proportion of variance)  

sigma1=$(head -1 $i | sed 's/,/\n/g' | cat -n | grep -n "sigma" | cut -f 1 |  sed 's/:/\n/g' | awk 'NR==1')
sigma2=$(head -1 $i | sed 's/,/\n/g' | cat -n | grep -n "sigma" | cut -f 1 |  sed 's/:/\n/g' | awk 'END{print $NF}')
A=$( echo $i | cut -d"/" -f3)
B=$( echo $A | cut -d_ -f1)
cat $i | cut -d ',' -f $sigma1-$sigma2 > ../../Sigma/Methylation/$A

beta1=$(head -1 $i | sed 's/,/\n/g' | cat -n | grep -n "beta" | cut -f 1 | sed 's/:/\n/g' | awk 'NR==1')
beta2=$(head -1 $i | sed 's/,/\n/g' | cat -n | grep -n "beta" | cut -f 1 | sed 's/:/\n/g' | awk 'END{print $NF}')
A=$( echo $i | cut -d"/" -f3)
B=$( echo $A | cut -d_ -f1)
cat $i | cut -d ',' -f $beta1-$beta2 > ../../Beta/Methylation/$A

comp1=$(head -1 $i | sed 's/,/\n/g' | cat -n | grep -n "comp" | cut -f 1 | sed 's/:/\n/g' | awk 'NR==1')
comp2=$(head -1 $i | sed 's/,/\n/g' | cat -n | grep -n "comp" | cut -f 1 | sed 's/:/\n/g' | awk 'END{print $NF}')
A=$( echo $i | cut -d"/" -f3)
B=$( echo $A | cut -d_ -f1)
cat $i | cut -d ',' -f $comp1-$comp2 > ../../Comp/Methylation/$A

done

head -n+1 ../../Sigma/Methylation/pack_years.csv > ../../Sigma/Methylation/head_pack_years.csv
head -n+1 ../../Beta/Methylation/pack_years.csv > ../../Beta/Methylation/head_beta.csv
head -n+1 ../../Comp/Methylation/pack_years.csv > ../../Comp/Methylation/head_comp.csv

tail -250 ../../Sigma/Methylation/pack_years.csv >>  ../../Sigma/Methylation/head_pack_years_temp.csv
tail -250 ../../Beta/Methylation/pack_years.csv >> ../../Beta/Methylation/head_beta_temp.csv
tail -250 ../../Comp/Methylation/pack_years.csv >> ../../Comp/Methylation/head_comp_temp.csv

cat ../../Sigma/Methylation/head_pack_years.csv ../../Sigma/Methylation/head_pack_years_temp.csv > ../../Sigma/Methylation/pack_years_17833_processed.csv
cat ../../Beta/Methylation/head_beta.csv ../../Beta/Methylation/head_beta_temp.csv > ../../Beta/Methylation/pack_years_17833_processed.csv
cat ../../Comp/Methylation/head_comp.csv ../../Comp/Methylation/head_comp_temp.csv > ../../Comp/Methylation/pack_years_17833_processed.csv

# execute only when I make sure that the processed files are fine. Correct.
rm ../*/*temp.csv
rm ../*/head*.csv
