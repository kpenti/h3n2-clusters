'use strict';

const Hapi = require('hapi');
var _ = require('lodash');
var Path = require('path');
var Inert = require('inert');

const server = new Hapi.Server();

server.register(Inert, function () {
    server.connection({ port: 3000 });
    server.route({
        method: 'GET',
        path: '/{path*}',
        handler: {
            directory: {
                path: './public',
                listing: true,
                index: true
            }
        }      
    });

    server.route({
        method: 'POST',
        path: '/submission',
        //path: '/{name}',
        handler: function (request, reply) {

            var Sequences = require("./libs/sequences.js");
            var Glycosylation = require("./libs/glycosylation.js");
            var Patches = require("./libs/patches.js");
            var Clusters = require("./libs/clusters.js");
            console.log("rawPayload: " + request.payload.sequence1);
            console.log("rawPayload: " + request.payload.sequence2);

            if(_.isUndefined(request.payload.sequence1) || _.isUndefined(request.payload.sequence2)) {
                var results = {};
                results.status = false;
                results.message = "Please enter two sequences";
                console.log("results 456: " + results);
                reply(results);
            }

            var sequences = new Sequences(request.payload.sequence1, request.payload.sequence2);
            var result = sequences.convert();
            console.log("results 45: " + JSON.stringify(result));
            if(result.status == false) {
                reply(results);
            } else {
                var sequence1 = 'X' + result.sequence1;
                var sequence2 = 'X' + result.sequence2;

                var patches = new Patches(sequence1, sequence2)
                var epitopes = patches.find_patches();
                
                //start glyco
                var glyco = new Glycosylation();
                var glycosylation_sites1 = glyco.find_glycosylation_sites(sequence1);
                var glycosylation_sites2 = glyco.find_glycosylation_sites(sequence2);
                var additions = glyco.getAdditions(glycosylation_sites1, glycosylation_sites2);
                var removals = glyco.getRemovals(glycosylation_sites1, glycosylation_sites2);
                var unglycosylated_epitopes = patches.remove_glycosylated_patches(glycosylation_sites2);
                //end glyco

                var clusters = new Clusters(unglycosylated_epitopes)
                var clusters_found = clusters.find_clusters();

                let results = {};
                results.sequence1 = patches.sequence1;
                results.sequence2 = patches.sequence2;
                results.mutations = patches.mutations;
                results.substitutions = patches.substitutions;
                results.patches = epitopes;
                results.clusters = clusters_found;
                results.glycosylation = {
                    'additions':additions,
                    'removals' :removals,
                    'sequence1':glycosylation_sites1,
                    'sequence2':glycosylation_sites2
                };
                results.status = true;

                reply(results);            
            }
        }
    });


    server.start(function() { console.log('Visit: http://127.0.0.1:3000') });
}); // requires a callback function but can be blank
