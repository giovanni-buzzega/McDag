
#!/bin/bash
# Download hiv dataset
if [ ! -d hiv-dataset ] 
then
        echo "Downloading hiv"
	mkdir hiv-dataset
	wget https://webdocs.cs.ualberta.ca/~ghlin/src/WebTools/HIV/42.tgz
	tar -xzvf 42.tgz -C ./hiv-dataset --strip-components=1 && rm 42.tgz
	wget -P hiv-dataset https://webdocs.cs.ualberta.ca/~ghlin/src/WebTools/HIV/CPZ+AF447763.fasta
fi
# Download fish dataset
if [ ! -d fish-dataset ] 
then
        echo "Downloading fish"
	mkdir fish-dataset
	wget -P fish-dataset/ https://afproject.org/media/genome/std/assembled/fish_mito/dataset/assembled-fish_mito.zip
	unzip -d fish-dataset -j fish-dataset/assembled-fish_mito.zip && rm fish-dataset/assembled-fish_mito.zip
fi
# Compile mcdag
if [ ! -e mcdag ]
then
        echo "Compiling McDag"
	g++ -O3 -march=native -o mcdag mcdag.cpp
fi

chars=3000

./mcdag -k 1 -s 4 -n $chars | head -n 1 | tr -d '#$' | tr 'B' 'T' | tr 'D' 'G' > ./random_string
inputfish1='./fish-dataset/NC_009057.fasta'
inputhiv1='./hiv-dataset/*AF005496.fasta'

STRINGhiv1=$(tail -n +2 $inputhiv1 | tr -d '\n')
STRINGfish1=$(tail -n +2 $inputfish1 | tr -d '\n')
STRINGrand1=$(head -n 1 ./random_string | tr -d '\n')

mkdir -p results/hist
# Histogram
echo "computing histograms"
output_file=./results/hist/hist_mix_"$chars"_rand_AF005496.txt
[ -e "$output_file" ] || ./mcdag -r -l -z  $(cut -c 1-$chars <<< $STRINGrand1) $(cut -c 1-$chars <<< $STRINGhiv1) > "$output_file"
output_file=./results/hist/hist_mix_"$chars"_rand_NC_009057_time.txt
[ -e "$output_file" ] || ./mcdag -r -l -z  $(cut -c 1-$chars <<< $STRINGrand1) $(cut -c 1-$chars <<< $STRINGfish1) > "$output_file"
output_file=./results/hist/hist_mix_"$chars"_AF005496_NC_009057.txt
[ -e "$output_file" ] || ./mcdag -r -l -z  $(cut -c 1-$chars <<< $STRINGhiv1) $(cut -c 1-$chars <<< $STRINGfish1) > "$output_file"


sleep 1
./mcdag -k 1 -s 4 -n $chars | head -n 1 | tr -d '#$' | tr 'B' 'T' | tr 'D' 'G' >> ./random_string
inputfish2='./fish-dataset/NC_011179.fasta'
inputhiv2='./hiv-dataset/*K03454.fasta'

STRINGhiv2=$(tail -n +2 $inputhiv2 | tr -d '\n')
STRINGfish2=$(tail -n +2 $inputfish2 | tr -d '\n')
STRINGrand2=$(tail -n 1 ./random_string | tr -d '\n')
output_file=./results/hist/hist_abs_"$chars"_rand.txt
[ -e "$output_file" ] || ./mcdag -r -z  -l  $(cut -c 1-$chars ./random_string) > "$output_file"
output_file=./results/hist/hist_abs_"$chars"_AF005496_K03454.txt
[ -e "$output_file" ] || ./mcdag -r -z  -l  $(cut -c 1-$chars <<< $STRINGhiv1) $(cut -c 1-$chars <<< $STRINGhiv2) > "$output_file"
output_file=./results/hist/hist_abs_"$chars"_NC_009057_NC_011179.txt
[ -e "$output_file" ] || ./mcdag -r -z  -l  $(cut -c 1-$chars <<< $STRINGfish1) $(cut -c 1-$chars <<< $STRINGfish2) > "$output_file"


# hiv curves
echo "Computing curves for hiv dataset"
mkdir -p results/curves

length=$(echo -n $STRINGhiv1 | wc -c)
percentage_hiv_csa () {
	percent=$1
	chars=$((length * percent / 100))

	output_file=./results/curves/curves_hiv_csa_"$percent"_AF005496_K03454.txt
	[ -e "$output_file" ] || ./mcdag -r -z -i $(cut -c 1-$chars <<< $STRINGhiv1) $(cut -c 1-$chars <<< $STRINGhiv2) > "$output_file" 
}
percentage_hiv () {
	percent=$1
	chars=$((length * percent / 100))

	output_file=./results/curves/curves_hiv_mcdag_"$percent"_AF005496_K03454.txt
	[ -e "$output_file" ] || ./mcdag -r -z  $(cut -c 1-$chars <<< $STRINGhiv1) $(cut -c 1-$chars <<< $STRINGhiv2) > "$output_file" 
}


startfrom=1
# do while
while true; do
	[ $(( ( $(nproc) - 1 ) / 2 )) -lt 63 ] && chunk=$(( ( $(nproc) - 1 ) / 2 )) ||  chunk=63
	[ $chunk -lt 1 ] && chunk=1
	[ 63 -lt $(( startfrom + chunk - 1)) ] && endat=63 || endat=$(( startfrom + chunk - 1 ))

	for i in $(seq $startfrom $endat )
	do
		percentage_hiv_csa $i &
		pids[${i}]=$!
		percentage_hiv $i &
	done
	# wait for all pids
	for pid in ${pids[*]}; do
		wait $pid
	done
	
	startfrom=$(( startfrom + chunk ))

	if [ $startfrom -le 63 ]; then
		continue
	else
		break
	fi
done

startfrom=64 
# do while 
while true; do
	[ $(( $(nproc) - 1 )) -lt 100 ] && chunk=$(( $(nproc) - 1 )) ||  chunk=100
	[ $chunk -lt 1 ] && chunk=1
	[ 100 -lt $(( startfrom + chunk - 1)) ] && endat=100 || endat=$(( startfrom + chunk - 1 ))

	for i in $(seq $startfrom $endat )
	do
		percentage_hiv $i &
		pids[${i}]=$!
	done
	# wait for all pids
	for pid in ${pids[*]}; do
		wait $pid
	done
	
	startfrom=$(( startfrom + chunk ))

	if [ $startfrom -le 100 ]; then
		continue
	else
		break
	fi
done

# increasing alphabet size:
echo "Computing curves for increasing alphabet; random strings"
mkdir -p results/sigma

startfrom=2
# do while 
while true; do
	[ $(( $(nproc) - 1 )) -lt 93 ] && chunk=$(( $(nproc) - 1 )) ||  chunk=93
	[ $chunk -lt 1 ] && chunk=1
	[ 93 -lt $(( startfrom + chunk - 1)) ] && endat=93 || endat=$(( startfrom + chunk - 1 ))

	for i in $(seq $startfrom $endat )
	do
		./mcdag  -z -r -s $i -n 3000 > ./results/sigma/curves_rand_sigma_mcdag_"$i".txt & 
		pids[${i}]=$!
	done
	# wait for all pids
	for pid in ${pids[*]}; do
		wait $pid
	done
	
	startfrom=$(( startfrom + chunk ))

	if [ $startfrom -le 93 ]; then
		continue
	else
		break
	fi
done












# build csv
echo building csv files
sed ./results/curves/curves_hiv_csa_1[_.]* -e 's/lcs-1/lcs_minus_one/g' -e 's/d[0-9]/dname/g' | grep -o -- "-\?[0-9]\+\.\?[0-9]*e\?[+-]\?[0-9]*" | tr '\n' '\t' | tr '\t' '\n' | sed -e 's/-\?[0-9]\+\.\?[0-9]*e\?[+-]\?[0-9]*/-1/g' > ./results/curves/default_csa
sed ./results/curves/curves_hiv_mcdag_1[_.]* -e 's/lcs-1/lcs_minus_one/g' -e 's/d[0-9]/dname/g' | grep -o -- "-\?[0-9]\+\.\?[0-9]*e\?[+-]\?[0-9]*" | tr '\n' '\t' | tr '\t' '\n' | sed -e 's/-\?[0-9]\+\.\?[0-9]*e\?[+-]\?[0-9]*/-1/g' > ./results/curves/default_mcdag



echo "|X|	|Y|	num_matches	csa time	csa_all nodes	csa_all edges	csa_maximal time	csa_maximal nodes	csa_maximal edges	#MCS	mcs_min time	mcs_min nodes	mcs_min edges	#lcs-1	#lcs	lcs len	|X|	|Y|	num_matches	csa_filter time	csa_filter nodes	csa_filter edges	csa_filter_codet time	csa_filter_codet nodes	csa_filter_codet edges	mcdag time	mcdag nodes	mcdag edges	#MCS	mcs_min time	mcs_min nodes	mcs_min edges	#(lcs-1)	#lcs	lcs len" > ./results/curves/all_data.csv
for i in $(seq 100)
do 
	[ -e ./results/curves/curves_hiv_csa_${i}[_.]* ] && csa="./results/curves/curves_hiv_csa_${i}[_.]*" || csa="./results/curves/default_csa"
	[ -e ./results/curves/curves_hiv_mcdag_${i}[_.]* ] && mcdag="./results/curves/curves_hiv_mcdag_${i}[_.]*" || mcdag="./results/curves/default_mcdag"
	cat $csa $mcdag   | sed -e 's/lcs-1/lcs_minus_one/g' -e 's/d[0-9]/dname/g' | grep -o -- "-\?[0-9]\+\.\?[0-9]*e\?[+-]\?[0-9]*" | sed -e 's/-1/nan/' | tr '\n' '\t' | tr '.' ',' >> ./results/curves/all_data.csv
	echo "" >> ./results/curves/all_data.csv
done



echo "|X|	|Y|	num_matches	csa_filter time	csa_filter nodes	csa_filter edges	csa_filter_codet time	csa_filter_codet nodes	csa_filter_codet edges	mcdag time	mcdag nodes	mcdag edges	#MCS	mcs_min time	mcs_min nodes	mcs_min edges	#(lcs-1)	#lcs	lcs len" > ./results/sigma/all_data.csv
for i in $(seq 2 93)
do 
	[ -e ./results/sigma/curves_rand_sigma_mcdag_${i}[_.]* ] && mcdag="./results/sigma/curves_rand_sigma_mcdag_${i}[_.]*" || mcdag="./results/curves/default_mcdag"
	cat $mcdag   | sed -e 's/lcs-1/lcs_minus_one/g' -e 's/d[0-9]/dname/g' | grep -o -- "-\?[0-9]\+\.\?[0-9]*e\?[+-]\?[0-9]*" | sed -e 's/-1/nan/' | tr '\n' '\t' | tr '.' ',' >> ./results/sigma/all_data.csv
	echo "" >> ./results/sigma/all_data.csv
done
