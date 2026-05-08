// Identify the provided selected class
var growth_condition = gc_selector.value;

// Given promoter selection, restrict selectable growth_conditions to those tested with promoter
var gc_indices = getAllIndexes(data.data['growth_condition'], growth_condition);
var replicate_menu = [];

for (var i=0; i < gc_indices.length; i++) { 
    replicate_menu.push(data.data['replicate'][gc_indices[i]]);
}
replicate_menu.sort()
rep_selector.options =  replicate_menu.filter(onlyUnique);;
rep_selector.value = replicate_menu[0];
rep_selector.change.emit();

// Custom Function Definitions
function getAllIndexes(arr, val) {
    var indices = [], i = -1;
    while ((i = arr.indexOf(val, i+1)) != -1){
        indices.push(i);
    }
    return indices;
}

function onlyUnique(value, index, array) {
    return array.indexOf(value) === index;
  }
  