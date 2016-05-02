'use strict';
// patches.js
// ========
//
var _ = require('lodash');
var fs = require('fs');

module.exports = Patches;

function Patches(sequence1, sequence2) {
	this.sequence1 = sequence1;
	this.sequence2 = sequence2;
	this.init();
}

Patches.prototype.init = function () {
	this.substitutions = this.set_substitutions();
	this.mutations = this.set_mutations();
    
    this.pdb = JSON.parse(fs.readFileSync('data/1HGD.json', 'utf8'));
}

Patches.prototype.set_substitutions = function () {
	let seq1 = this.sequence1.split("");
	let seq2 = this.sequence2.split("");

	var substitutions = [];
	for (let i = 0; i < seq1.length; i++) {
		if(seq1[i] != seq2[i]) {
			substitutions.push(i);
		}
	}
	return substitutions;
}

Patches.prototype.set_mutations = function () {
	var mutations = [];
	this.substitutions.forEach(function(substitution) {
		if(this.sequence1[substitution]) {
			mutations.push(this.sequence1[substitution] + substitution + this.sequence2[substitution]);
		}
	}, this);
	return mutations;
}

//find the patches contining at least 3 mutations
Patches.prototype.find_patches = function () {
	var fs = require('fs');
	var self = this;
	var patch_list = JSON.parse(fs.readFileSync('data/patch_list_HA1.json', 'utf8'));
	var identified_patches = [];
	_(patch_list).forEach(function(residues, centroid) {
		var residues = residues.map(function(obj){ 
   			return parseInt(obj.replace("A_", ""));
		});
		if(_.intersection(self.substitutions, residues).length >= 3) {
			identified_patches.push({residues: residues, seed: parseInt(centroid.replace("A_", "")), centroid: parseInt(centroid.replace("A_", ""))});
		}
	});
    self.identified_patches = identified_patches;
	return identified_patches;
}

Patches.prototype._dist_to_residue = function (resi_i, residues) {
    var self = this;
    var dists = [];
    _(residues).forEach(function(resi) {
        dists.push(self.pdb[resi_i][resi]);
    });
    return dists[0];
}

Patches.prototype.reduce_patch_size = function (glyco, row) {
    var self = this;
    var dists = [];
    _(glyco).forEach(function(glyco_resn) {
        dists.push(self.pdb[glyco_resn][row.seed]);
    });
    
    var mimimum_patch_size = _.max(dists);
    var patch_residues = [];
    _(row.residues).forEach(function(residue) {
        let dist_to_seed = self.pdb[row.seed][residue];
        if(dist_to_seed < mimimum_patch_size) {
            patch_residues.push( residue );
        }
    });
    row.residues = patch_residues;
    return row;
}

Patches.prototype.remove_glycosylated_patches = function (glycosylated_sites) {
    var self = this;
    
    _(self.identified_patches).forEach(function(row, centroid) {
        let glyco = _.intersection(glycosylated_sites, row.residues);
        if(glyco.length > 0) {
            self.identified_patches[centroid] = self.reduce_patch_size(glyco, row);     //reduce size of patch
        }
    });
    
    //now remove those with less than 3 mutations
    var identified_patches = [];
    _(self.identified_patches).forEach(function(row, centroid) {
        if(_.intersection(self.substitutions, row.residues).length >= 3) {
			identified_patches.push(row);
		}    
    });
    self.identified_patches = identified_patches;
	return self.identified_patches;
}


/*
module.exports = {
	generate_patches: function (request, reply) {
		//console.log(request);
		//calculate patches
		var fs = require('fs');
		var patch_list = JSON.parse(fs.readFileSync('data/patch_list_HA1.json', 'utf8'));

		sequence1 = "XQDFPGNDNSTATLCLGHHAVPNGTLVKTITNDQIEVTNATELVQSSSTGKICNNPHRILDGINCTLIDALLGDPHCDGFQNETWDLFVERSKAFSNCYPYDVPDYASLRSLVASSGTLEFINEGFTWTGVTQNGGSNACKRGPDSGFFSRLNWLYKSGSAYPVLNVTMPNNDNFDKLYIWGVHHPSTDQEQTNLYVQTSGRVTVSTKRSQQTIIPNIGSRPWVRGLSSRISIYWTIVKPGDILVINSNGNLIAPRGYFKMRTGKSSIMRSDAPIGTCISECITPNGSIPNDKPFQNVNKITYGACPKYVKQNTLKLATGMRNVPEKQT"

sequence2 = "XQDLPGNDNSTATLCLGHHAVPNGTLVKTITNDQIEVTNATELVQSSSTGKICDNPHRILDGINCTLIDALLGDPHCDGFQNEKWDLFVERSKAFSNCYPYDVPDYASLRSLVASSGTLEFINEGFNWTGVTQNGGSSACKRGPDNGFFSRLNWLYKLGSTYPVQNVTMPNNDNSDKLYIWGVHHPSTDKEQTDLYVQASGKVTVSTKRNQQTVIPNVGSRPWVRGLSSRVSIYWTIVKPGDILVINSNGNLIAPRGYFKMRTGKSSIMRSDAPIGTCSSECITPNGSIPNDKPFQNVNKITYGACPKYVKQNTLKLATGMRNVPEKQT"
		calculate_substitutions()
		console.log(patch_list);
	},

	calculate_substitutions: function () {
	#set the substitutions
    def set_substitutions(self):
		self.substitutions = [i for i in xrange(len(self.sequence1)) if self.sequence1[i] != self.sequence2[i]]
		return self.substitutions
	},

	calculate_mutations: function () {
    	mutations = []
    	for substitution in self.substitutions :
			mutations.append(self.sequence1[substitution] + str(substitution) + self.sequence2[substitution])
        return = mutations
	}
};
*/