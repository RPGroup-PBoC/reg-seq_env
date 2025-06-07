// Get indexes for all promoter and gc data
var prom_inds = getAllIndexes(source_data['promoter'], prom);
var gc_inds = getAllIndexes(source_data['growth_condition'], gc);

var overlap_inds = prom_inds.filter(value => gc_inds.includes(value))

var x = new Array(prom_inds.length);

// go through conditions

var rep_indexes = [];
var inds = [];
var reps = [];

for (j=0; j<prom_inds.length; j=j+1){
    var mut_info = smooth_info(prom_inds[j], smooth_selector);
    var cv = calculateCV(mut_info);
    x[j] = cv;

    if (overlap_inds.includes(prom_inds[j])){
        reps.push(source_data['replicate'][prom_inds[j]]);
        rep_indexes.push(prom_inds[j]);
        inds.push(j);
    }
}

// get permutation that sorts array
var sorted_indexes = findSortingPermutation(x);

// assign NaN values to third replicate if not found
if (rep_indexes.length < 3){
    cv_point_gc['x_3'] = [NaN]
    cv_point_gc['y_3'] = [NaN]
}

// find coordinates for points
for (i=0; i<reps.length; i=i+1){
    cv_point_gc['x_'+reps[i]] = [calculateCV(smooth_info(rep_indexes[i], smooth_selector))];
    cv_point_gc['y_'+reps[i]] = [(sorted_indexes.indexOf(inds[i])+1)/prom_inds.length];
}


cv_gc['x'] = sorted_indexes.map(_x=>x[_x]);

cv_patches['x_cv_1'] = [cv_gc['x'][0], cv_gc['x'][0], 0.65, 0.65];
cv_patches['x_cv_3'] = [0.75, 0.75, cv_gc['x'].at(-1), cv_gc['x'].at(-1)];



cds_cv_gc.change.emit();
cds_cv_point_gc.change.emit();





for (i=1;i<4;i=i+1){
    // go through promoters
    let rep_inds = getAllIndexes(source_data['replicate'], ""+i);
    let filteredArray_reps = gc_inds.filter(value => rep_inds.includes(value));
    
    if (filteredArray_reps.length == 0) {
        cv_point_prom['x_'+i] = [NaN];
        cv_point_prom['y_'+i] = [NaN];
        let l = cv_prom['x_'+i].length
        cv_prom['x_'+i] = new Array(l).fill(NaN);
        cv_patches['x_cv_1_'+i] = [NaN, NaN, NaN, NaN];
        cv_patches['x_cv_3_'+i] = [NaN, NaN, NaN, NaN];
    }
    else{
        let x = new Array(filteredArray_reps.length);
        let rep;
        let rep_index;
        let inds
        for (j=0;j<filteredArray_reps.length;j=j+1){
            let mut_info = smooth_info(filteredArray_reps[j], smooth_selector);
            let cv = calculateCV(mut_info);
            
            x[j] = cv;

            if (overlap_inds.includes(filteredArray_reps[j])){
                rep = (source_data['replicate'][filteredArray_reps[j]]);
                rep_index = filteredArray_reps[j];
                inds = j;
            }
        }
       
        let sorted_indexes = findSortingPermutation(x);

        cv_point_prom['x_'+rep] = [calculateCV(smooth_info(rep_index, smooth_selector))];
        cv_point_prom['y_'+rep] = [(sorted_indexes.indexOf(inds)+1)/x.length];

        if (overlap_inds.includes(filteredArray_reps[j])){
            reps.push(source_data['replicate'][filteredArray_reps[j]]);
            rep_indexes.push(filteredArray_reps[j]);
            inds.push(j);
            
        }
        cv_prom['x_'+i] = sorted_indexes.map(_x=>x[_x]);
        cv_patches['x_prom_1_'+i] = [cv_prom['x_'+i][0], cv_prom['x_'+i][0], 0.65, 0.65];
        cv_patches['x_prom_3_'+i] = [0.75, 0.75, cv_prom['x_'+i].at(-1), cv_prom['x_'+i].at(-1)];
    }
}




cds_cv_gc.change.emit();
cds_cv_point_gc.change.emit();
cds_cv_patches.change.emit();
cds_cv_prom.change.emit();
cds_cv_point_prom.change.emit();



function findSortingPermutation(arr) {
    const n = arr.length;
    const permutation = [];
  
    // Create an array of indices and sort it based on the values in the original array
    for (let i = 0; i < n; i++) {
      permutation.push(i);
    }
    permutation.sort((a, b) => arr[a] - arr[b]);
  
    return permutation;
}