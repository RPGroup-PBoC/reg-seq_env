// Identify the provided selected class
var promoter = prom_selector.value;

// Update metadata
var meta_inds = getAllIndexes(meta.data['promoter'], promoter);
var direction = meta.data['direction'][meta_inds[0]];
var promoter_seq = meta.data['promoter_seq'][meta_inds[0]];
var genes = meta.data['genes'][meta_inds[0]];
var five_prime = meta.data['five_prime'][meta_inds[0]];
var three_prime = meta.data['three_prime'][meta_inds[0]];

if (promoter.includes("_predicted")) {
    prom_desc.text = '<div style="width:300px; overflow-wrap: break-word;"><b> Genes controlled by putative promoter</b>: <br/>' + genes + '<br/><b>Strand: </b><br/>' + direction + '<br/><b>5\':</b><br/>' + five_prime + '<br/><b>3\':</b><br/>' + three_prime + '</div>';
  } else {
    prom_desc.text = '<div style="width:300px; overflow-wrap: break-word;"><b> Genes controlled by promoter</b>: <br/>' + genes + '<br/><b>Strand: </b><br/>' + direction + '<br/><b>5\':</b><br/>' + five_prime + '<br/><b>3\':</b><br/>' + three_prime + '</div>';
}

regulonDB_desc.text = '<div style="width:700px;"><b> Annotation in RegulonDB</b><br/>';


var regulon_indices = getAllIndexes(regulonDB.data['PROMOTER_NAME'], promoter);
if (regulon_indices.length == 0){
    regulonDB_desc.text += '<br/>No Binding Sites Found';
} else{
    for (var i=0; i < regulon_indices.length; i++) {
        regulonDB_desc.text += '<div style="overflow-wrap: break-word;"><br/><b>' + regulonDB.data['RI_FUNCTION'][regulon_indices[i]] + '</b><br/>Transcription Factor: ' + regulonDB.data['TRANSCRIPTION_FACTOR_NAME'][regulon_indices[i]] + '<br/>Binding Site Position Relative to TSS: ' + regulonDB.data['CENTER_POSITION'][regulon_indices[i]] + '<br/> Binding Site Sequence (Capital Letters): ' + regulonDB.data['RI_SEQUENCE'][regulon_indices[i]] + '<br/> Consensus Sequence: ' + regulonDB.data['CONSENSUS_SEQUENCE'][regulon_indices[i]] + '</div>';
    }
}
regulonDB_desc.text += '</div>'



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
