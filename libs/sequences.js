'use strict';
// sequences.js
// ========
//
var _ = require('lodash');

module.exports = Sequences;

function Sequences(sequence1, sequence2) {
	this.sequence1 = sequence1;
	this.sequence2 = sequence2;
	this.ha1_sequence = "QDLPGNDNSTATLCLGHHAVPNGTLVKTITDDQIEVTNATELVQSSSTGKICNNPHRILDGIDCTLIDALLGDPHCDVFQNETWDLFVERSKAFSNCYPYDVPDYASLRSLVASSGTLEFITEGFTWTGVTQNGRSNACKRGPGSGFFSRLNWLTKSGSTYPVLNVTMPNNDNFDKLYIWGIHHPSTNQEQTSLYVQASGRVTVSTRRSQQTIIPNIGSRPWVRGLSSRISIYWTIVKPGDVLVINSNGNLIAPRGYFKMRTGKSSIMRSDAPIDTCISECITPNGSIPNDKPFQNVNKITYGACPKYVKQNTLKLATGMRNVPEKQT";
}

Sequences.prototype.findSimilarity = function (target) {
    var count = 0;
    for(let i = 0; i <= this.ha1_sequence.length; i++) {
        if(this.ha1_sequence[i] == target[i]) {
            count = count + 1;
        }
    }
    return count;
}

Sequences.prototype.getHA1Sequence = function (sequence) {
    
    //generate all variations
    var sequences = [];
    for(let i = 0; i <= sequence.length; i++) {
        let target = sequence.substr(i, 328);
        if(target.length == 328) {
            let similarity = this.findSimilarity(target); 
            sequences.push({target:target, similarity:similarity});
        }
    }
    
    if(sequences.length == 0) {
        return false;
    } else {
        var bestSequence = _.first(_.sortBy(sequences, function(o) { return -o.similarity; }));
        return bestSequence.target;
    }
}

Sequences.prototype.convert = function () {
    //make sure both sequences are HA1 sequences
	this.sequence1 = this.getHA1Sequence(this.sequence1);
	this.sequence2 = this.getHA1Sequence(this.sequence2);
    
    return {status: (this.sequence1 && this.sequence2), sequence1:this.sequence1, sequence2:this.sequence2};
}
