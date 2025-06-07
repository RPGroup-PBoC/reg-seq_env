// Get data for footprints
var prom_inds = getAllIndexes(source_data['promoter'], prom);
var gc_inds = getAllIndexes(source_data['growth_condition'], gc);
var filteredArray = prom_inds.filter(value => gc_inds.includes(value));

// Get data for expression shifts
var prom_inds_ex = getAllIndexes(exshift_data['promoter'], prom);
var gc_inds_ex = getAllIndexes(exshift_data['growth_condition'], gc);


let prom_inds_hmm = getAllIndexes(hmm_data['promoter'], prom);
let gc_inds_hmm = getAllIndexes(hmm_data['growth_condition'], gc);
let filteredArray_hmm = prom_inds_hmm.filter(value => gc_inds_hmm.includes(value))


for (var i=1; i<4; i=i+1){
    // check if there is data for replicate
    if (filteredArray.map(x=>source_data['replicate'][x]).includes(i.toString())){
        // get hmm results for replicate
        let reps_inds_hmm = getAllIndexes(hmm_data['replicate'], i.toString());
        let filteredArray_hmm_rep = reps_inds_hmm.filter(value => filteredArray_hmm.includes(value));

       let hmm1;
       let hmm2;
        if (filteredArray_hmm_rep.length > 0){
            hmm2 = hmm_data['state_sequence'][filteredArray_hmm_rep];
            hmm1 = invertBinaryArray(hmm_data['state_sequence'][filteredArray_hmm_rep]);
        }
         // if no hmm were calcuated, record all as noise
        else {
            hmm1 = new Array(160).fill(1);
            hmm2 = new Array(160).fill(0);
        }

        display = mut_info_plot(i.toString(), smooth_selector, display, hmm1, hmm2);

        var rep_inds_ex = getAllIndexes(exshift_data['replicate'], ""+i);
        var filteredArray_ex = prom_inds_ex.filter(value => gc_inds_ex.includes(value))
        var filteredArray_ex = filteredArray_ex.filter(value => rep_inds_ex.includes(value))
        
        ex_display['pos_'+i] = exshift_data['pos'][filteredArray_ex];
        ex_display['base_'+i] = exshift_data['base'][filteredArray_ex];
        ex_display['wt_base_'+i] = exshift_data['wt_base'][filteredArray_ex];
        ex_display['expression_shift_'+i] = exshift_data['expression_shift'][filteredArray_ex];
    }
    
    else {
        display['mut_info_'+i+'_1'] = new Array(display['pos_'+i+'_1'].length).fill(0);
        display['mut_info_'+i+'_2'] = new Array(display['pos_'+i+'_2'].length).fill(0);
        ex_display['expression_shift_'+i] =new Array(ex_display['expression_shift_'+i].length).fill(0);
        ex_display['wt_base_'+i] = new Array(ex_display['wt_base_'+i].length).fill(" ");
    }

    const map1 = new Map();
    for (var ii = 0; ii < ex_display['wt_base_'+i].length; ii=ii+4) {
        var j = (ii)/4 - 115;
        map1.set(j, ex_display['wt_base_'+i][ii]);
    }
    x_axis[i-1].major_label_overrides = map1;
    p[i-1].reset.emit();
    p[i-1].change.emit();  
}




//console.log(display)
data_display.change.emit();
exshift_display.change.emit();
