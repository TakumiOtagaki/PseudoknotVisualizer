# TO RUN: bash CLI_test.sh
workdir="/large/otgk/PseudoknotVisualizer"
test_dir=$workdir/test
python="/home/otgk/.conda/envs/pymol/bin/python"
CLI=$workdir/CLI_PseudoknotVisualizer.py


pdb=$test_dir/1KPD.pdb
pdb_id=1KPD
chain=A
model=0
# format: chimera or pymol
# formats=("chimera", "pymol")
# formats=("chimera")
formats=("pymol")

for format in "${formats[@]}"; do
    echo "Processing with $format"
    output=$test_dir/coloring_$pdb_id.$model.$chain.$format.txt
    echo "$python $CLI -i $pdb  -o $output -c $chain -f $format -m $model"
    # $python $CLI -i $pdb  -o $output -c $chain -f $format -m $model
done