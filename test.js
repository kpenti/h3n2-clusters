'use strict';

var Sequences = require("./libs/sequences.js");
var Glycosylation = require("./libs/glycosylation.js");
var Patches = require("./libs/patches.js");
var Clusters = require("./libs/clusters.js");
var _ = require('lodash');
var request = {};
request.payload = {};
request.payload.sequence1 = "QDLPGNDNSKATLCLGHHAVPNGTLVKTITDDQTEVTNATELVQSSSTGKICNNPHRILDGIDCTLIDALLGDPHCDVFQNETWDLFVERSKAFSNCYPYDVPDYASLRSLVASSGTLEFITEGFTWTEVTENGGSNACKRGPDSGFFSRLNWLTKSGRTYPVLNVTMPNNDNFDKLYIWGVHHPSTDQEQTSLYVQTSGRVTVSTRRSQQTIIPNIGSRPWVRGLSSRISIYWTIVKPGDVLVINSNGNLIAPRGYFKMRTGKSSIMRSDAPIDTCISECITPNGSIPNDKPFQNVNKITYGACPKYVKQNTLKLATGMRNVPEKQT";
request.payload.sequence2 = "QKLPGNDNSTATLCLGHHAVPNGTIVKTITNDQIEVTNATELVQSSSTGEICDSPHQILDGKNCTLIDALLGDPHCDGFQNKKWDLFVERSKAYSNCYPYDVPDYASLRSLVASSGTLEFNNESFNWTGVTQNGTSSACIRRSKNSFFSRLNWLTHLNFKYSALNVTMPNNEQFDKLYIWGVHHPGTDKDQIFLYAQASGRITVSTKRSQQTVIPNIGSRLRVRNIPSRISIYWTIVKPGDILLINSTGNLIAPRGYFKMQSGKSSIMRSDAPIGKCNSKCITPNGSIPNDKPFQNVNRITYGACPRYVKQNTLKLATGMRNVPEKQT";

var sequences = new Sequences(request.payload.sequence1, request.payload.sequence2);
var result = sequences.convert();

var sequence1 = 'X' + result.sequence1;
var sequence2 = 'X' + result.sequence2;

var patches = new Patches(sequence1, sequence2)
var epitopes = patches.find_patches();

var glyco = new Glycosylation();
var glycosylation_sites1 = glyco.find_glycosylation_sites(sequence1);
var glycosylation_sites2 = glyco.find_glycosylation_sites(sequence2);
var additions = glyco.getAdditions(glycosylation_sites1, glycosylation_sites2);
var removals = glyco.getRemovals(glycosylation_sites1, glycosylation_sites2);
 
var unglycosylated_epitopes = patches.remove_glycosylated_patches(glycosylation_sites2);

var clusters = new Clusters(unglycosylated_epitopes)
var clusters_found = clusters.find_clusters();

var results = {};
results.sequence1 = patches.sequence1;
results.sequence2 = patches.sequence2;
results.mutations = patches.mutations;
results.substitutions = patches.substitutions;
results.patches = epitopes;
results.clusters = clusters_found;
results.status = true;

console.log(JSON.stringify(results));            
