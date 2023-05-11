# Preprocessing, denoise, head motion and eddy
dwidenoise $1/DTI/raw.nii $1/DTI/data_denoised.nii
fslroi $1/DTI/data_denoised.nii $1/DTI/nodif_raw.nii.gz 0 1 
bet2 $1/DTI/nodif_raw.nii.gz $1/DTI/nodif_brain.nii.gz -m -f 0.3 getindex.sh 35 $1/DTI/ 
eddy_openmp --imain=$1/DTI/data_denoised.nii --mask=$1/DTI/nodif_brain_mask.nii.gz --bvals=$1/DTI/raw.bval --bvecs=$1/DTI/raw.bvec --acqp=acqparams_one_phase.txt --index=$1/DTI/index.txt --out=$1/DTI/den_eddy --ref_scan_no=0 --ol_nstd=4
rm $1/DTI/den_eddy.eddy_*


# Run bedpost
bet $1/T1/raw.nii $1/T1/t1_brain.nii.gz
fast -t 1 -n 3 -b -B -o $1/T1/t1_brain_fast $1/T1/t1_brain.nii.gz
flirt -in $1/DTI/nodif_brain.nii.gz -ref $1/T1/t1_brain.nii.gz -out $1/DTI/dti_to_t1.nii.gz -omat $1/dti_to_t1.mat -interp nearestneighbour
convert_xfm -omat $1/t1_to_dti.mat -inverse $1/dti_to_t1.mat
flirt -in $1/T1/t1_brain_fast_pve_2.nii.gz -init $1/t1_to_dti.mat -ref $1/DTI/nodif_brain.nii.gz -out $1/DTI/WM_to_dti.nii.gz -applyxfm -interp nearestneighbour
mv $1/DTI/raw.bvec $1/DTI/bvecs
mv $1/DTI/raw.bval $1/DTI/bvals
mv $1/DTI/data_denoised.nii $1/DTI/data.nii
cd $1
bedpostx_gpu DTI/. --nf=3 --fudge=1 --bi=3000 --model=1


# Perform transformations
# Generates WM masks in T1 space
robustfov -i $1/t1/t1.nii -r $1/t1/t1_crop.nii.gz
bet $1/t1/t1_crop.nii.gz $1/t1/t1_brain.nii.gz
fast -t 1 -n 3 -b -B -o $1/t1/t1_brain_fast $1/t1/t1_brain.nii.gz
fslmaths $1/t1/t1_brain_fast_pve_2.nii.gz -bin $1/t1/wm_in_t1_mask.nii.gz

# Creates transformations from T1 to DTI
mkdir $1/transformations
flirt -in $1/dti/nodif_brain.nii.gz -ref $1/t1/t1_brain.nii.gz -out $1/dti/nodif_brain_in_t1.nii.gz -omat $1/transformations/dti_to_t1.mat -interp nearestneighbour
convert_xfm -omat $1/transformations/t1_to_dti.mat -inverse $1/transformations/dti_to_t1.mat

# Creates transformations from MNI to T1
flirt -in MNI152_T1_1mm_brain.nii -ref $1/t1/t1_brain.nii.gz -omat $1/transformations/mni_to_t1_affine.mat
fnirt --in=MNI152_T1_1mm_brain.nii --aff=$1/transformations/mni_to_t1_affine.mat --cout=$1/transformations/mni_to_t1_nonlinear_warp --ref=$1/t1/t1_brain.nii.gz

# Applies transformations
flirt -in $1/t1/wm_in_t1_mask.nii.gz -init $1/transformations/t1_to_dti.mat -ref $1/dti/nodif_brain.nii.gz -out $1/dti/wm_in_dti.nii.gz -applyxfm -interp nearestneighbour
applywarp --ref=$1/t1/t1_brain.nii.gz --in=BN_Atlas_246_1mm.nii --warp=$1/transformations/mni_to_t1_nonlinear_warp.nii.gz --out=$1/t1/bna_in_t1.nii.gz --interp=nn
flirt -in $1/t1/bna_in_t1.nii.gz -init $1/transformations/t1_to_dti.mat -ref $1/dti/nodif_brain.nii.gz -out $1/dti/bna_in_dti.nii.gz -applyxfm -interp nearestneighbour

# Performs tracktography and generates the FN-network
mkdir DTI/track
track -inputmodel bedpostx_dyad -curvethresh 45 -curveinterval 5 -bedpostxminf 0.1 -header $1/DTI/nodif_brain.nii.gz -seedfile $1/DTI/WM_to_dti.nii.gz -bedpostxdir $1/DTI.bedpostX -tracker rk4 -interpolator nn -stepsize 2 -outputfile $1/DTI/track_bedpost_camino
procstreamlines -inputfile $1/DTI/track_bedpost_camino -mintractlength 20 -maxtractlength 250 -header $1/DTI/nodif_brain.nii.gz > $1/DTI/track_bedpost_camino_post




