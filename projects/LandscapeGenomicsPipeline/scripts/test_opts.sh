DIRECTORY=$(pwd)

while getopts d:s:i:m:b:v:? opt 2>/dev/null
do
  case $opt in
    d) DIRECTORY=$OPTARG;;
    s) SAMPLES=$OPTARG;;
    i) VCFFILE=$OPTARG;;
    m) METADATA=$OPTARG;;
    b) BIOVAR=$OPTARG;;
    v) VARIABLE=$OPTARG
      # Check if -m is set, and if so, modify the behavior of -v
      if [ -n "$METADATA" ]; then
        echo "Option -m is set, so -v behavior is influenced."
      fi
      ;;
    ?) echo "Valid parameters are: [-l, -s]" ;;
  esac
done

if [ -z "$DIRECTORY" ] || [ -z "$SAMPLES" ] || [ -z "$VCFFILE" ]; then
  echo "Error: Missing required parameters. Please provide values."
  exit 1
fi

if [ -n "$DIRECTORY" ]; then
  echo "DIRECTORY is set to: $DIRECTORY"
fi

if [ -n "$SAMPLES" ]; then
  echo "SAMPLES is set to: $SAMPLES"
fi

if [ -n "$VCFFILE" ]; then
  echo "VCFFILE is set to: $VCFFILE"
fi

if [ -n "$METADATA" ]; then
  if [ -n "$BIOVAR" ]; then
    echo "METADATA is set to: $METADATA. Including Parameter BIOVAR: $BIOVAR"
  else
    echo "Error: BIOVAR is missing. Please set the BIOVAR parameter."
    exit 1
  fi
fi

if [ -n "$VARIABLE" ]; then
  echo "VARIABLE is set to: $VARIABLE"
fi


#####


####original run on 8 populations

file="/media/inter/ssteindl/FC/usecaserepo/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/results/2L/Subsampled_2L.recode.af.k10"

SNPs=(431 2588 2964 252 4107 4346 5215 6255 6383)

for num in ${SNPs[@]};
do 
#echo $num
num2=$(expr $num + 1)
#echo $num2
sed "${num2}q;d" $file
done


###result
#2L      1432739            0.265   0.315   0.154   0.19    0.049   0.183   0.32    0.286
#2L      7172247            0.136   0.014   0.033   0.0     0.0     0.035   0.0     0.0
#2L      7997875            0.048   0.079   0.0     0.029   0.065   0.049   0.281   0.119
#2L      1007147            0.0     0.0     0.044   0.102   0.028   0.042   0.0     0.0
#2L      10660174           0.0     0.0     0.0     0.017   0.013   0.0     0.065   0.044
#2L      11366053           0.064   0.0     0.0     0.018   0.071   0.111   0.032   0.026
#2L      13877398           0.286   0.203   0.188   0.068   0.304   0.153   0.286   0.135
#2L      17013869           0.03    0.0     0.0     0.0     0.0     0.0     0.024   0.165
#2L      17355982           0.0     0.0     0.0     0.016   0.0     0.0     0.0     0.0



###run of 160 popualtions
betai="/media/inter/ssteindl/FC/usecaserepo/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/NSATvsBIO1_4/results/2L/BAYPASS/BayPass_summary_betai_reg.out"
awk 'NR>=1 && NR <= 1000 {if ($10 == "15.35252978") {print $2}}' $betai

# Assuming you want to store the numbers in an array named "numbers"
numbers=($(awk 'NR>=1 && NR <= 1000 {if ($10 == "15.35252978") {print $2}}' "$betai"))

file="/media/inter/ssteindl/FC/usecaserepo/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/NSATvsBIO1_4/results/2L/Subsampled_2L.af.h5"
# Print the array elements
for number in "${numbers[@]}"; do
    #num2=$(expr $number + 1)
    sed "${number}q;d" $file >> /media/inter/ssteindl/FC/usecaserepo/uc3-drosophola-genetics/projects/LandscapeGenomicsPipeline/af_analysis_min1.txt
done




