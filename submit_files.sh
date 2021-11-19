
if [ $# -lt 3 ]
then
  echo "Need input directory name: $0 <indir> <startfolder> <endfolder> [nfolder]"
  exit 255
fi

# indir=/home/staff/leeuw179/JEWEL/$1
indir=$1
outdir=$1
firstfile=$2
lastfile=$3
nfile=${4:-10}
#firstfile=5
#lastfile=20

for ((ifile=$firstfile;ifile<=$lastfile;ifile+=$nfile))
do
let  lastfile_job=$ifile+$nfile-1
if [ $lastfile_job -gt $lastfile ]
then
	lastfile_job=$lastfile
fi
 qsub -V -N analyse_files_batch_EF_bkg -q long -o log -e log -v indir=$indir,infile=example.hepmc,firstfile=$ifile,lastfile=$lastfile_job analyze_files_batch_EF.sh

 # qsub -V -N analyse_files_batch_EnergyFlow -q generic -o log -e log -v indir=$indir,infile=example.hepmc,firstfile=$ifile,nfile=$nfile analyze_files_batch_EF.sh
done

