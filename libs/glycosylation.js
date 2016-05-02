'use strict';
// glycosylation.js
// ========
//
var _ = require('lodash');

module.exports = Glycosylation;

function Glycosylation() {
	
}

Glycosylation.prototype.find_glycosylation_sites = function (sequence) {
    //look for NXS/T!P
    var sequons = [];
    for(let i = 0; i <= sequence.length; i++) {
        let target = sequence.substr(i, 4);
        if(target.length == 4 && target[0] == 'N' && target[3] != 'P' && (target[2] == 'S' || target[2] == 'T')) {
            sequons.push(i);
        }
    }
    return sequons;
}

Glycosylation.prototype.getAdditions = function (residues1, residues2) {
    return _.difference(residues2, residues1); 
}

Glycosylation.prototype.getRemovals = function (residues1, residues2) {
    return _.difference(residues1, residues2); 
}
