'use strict';
// clusters.js
// ========
//
var _ = require('lodash');
var fs = require('fs');

module.exports = Clusters;

function Clusters(patches) {
	this.patches = patches;
	this.init();
}

Clusters.prototype.init = function () {
    this.pdb = JSON.parse(fs.readFileSync('data/1HGD.json', 'utf8'));
}

//now make clusters
Clusters.prototype._calculate_patch_density = function (patch, residue_density) {

	var density_scores = [];
    var sum = 0;
	_(patch.residues).forEach(function(patch_residue) {
        let tmp = {};
        tmp.patch_residue = patch_residue;
        tmp.density = residue_density[patch_residue];
        density_scores.push(tmp);
        sum += tmp.density;
	});
	
	return {'centroid': patch.centroid, 'density': density_scores, 'sum': sum};
    
}

Clusters.prototype._cluster_score = function (centroids, patch_densities) {
	var scores = [];
    _(centroids).forEach(function(centroid) {
		scores.push( _.find(patch_densities, { 'centroid': centroid }).sum );
	});    
	return _.mean(scores);
}

Clusters.prototype._find_core = function (residues) {
    var residues = _.reverse(_.sortBy(residues, 'density'));
	return _.first(residues).patch_residue;
}

Clusters.prototype._clusters_containing_patch = function (seed_i, cluster_patches) {
    var clusters = [];
    _(cluster_patches).forEach(function(patches_j, cluster_j) {        
        if(_.indexOf(_.keys(patches_j), seed_i) > -1) {
            clusters.push(cluster_j);
        }
    });
    clusters = _.uniq(clusters);
    return clusters;    
}

Clusters.prototype._dist_to_residue = function (resi_i, residues) {
    var self = this;
    var dists = []
    _(residues).forEach(function(resi) {
        dists.push(self.pdb[resi_i][resi]);
    });
    //console.log('dists', dists.sort(sortNumber));
    return dists[0];
}
function sortNumber(a,b) {
    return a - b;
}
Clusters.prototype._get_seed_distances = function (seed, cores) {
    var self = this;
    var dists = [];
    _(cores).forEach(function(core) {
        let core_residues = core.split(',');
        dists.push( {'core':core, dist:self._dist_to_residue(seed, core_residues)} );
    });
    return _.sortBy(dists, 'dist');
}

Clusters.prototype._remove_overlapping_patches = function (cluster_i, cluster_patches) {
    var self = this;
    var patches_i = cluster_patches[cluster_i];    
    var patches_to_keep = [];
    var patches_to_remove = [];
    _(patches_i).forEach(function(residues_i, seed_i) {
        //in which clusters does this patch also exist?
        let clusters_to_check = self._clusters_containing_patch(seed_i, cluster_patches);    
        
        
        if(clusters_to_check.length > 1) {
            let seed_dists = self._get_seed_distances(seed_i, clusters_to_check);

            if( cluster_i != _.first(seed_dists).core ) {
                patches_to_remove.push( seed_i );
            }            
        }
	});

    //patches_full = _.sortBy( _.find(cluster_patches, {centroid: (cluster_i)}).residues, function( val ){ return val; } );
    if(patches_to_remove.length > 0) {
		patches_to_keep = {};
        _(patches_i).forEach(function(residues_i, seed_i) {
            if(_.indexOf(patches_to_remove, seed_i) === -1) {
                patches_to_keep[seed_i] = residues_i
            }
        });
	}
	return patches_to_keep;
    
}
/*
function remove_overlapping_patches($cluster_i, $cluster_patches) {
	$patches_i = $cluster_patches[$cluster_i];
	$patches_to_keep = $cluster_patches[$cluster_i];
	$patches_to_remove = [];
	foreach($patches_i as $seed_i => $residues_i) {
		//in which clusters does this patch also exist?
		$clusters_to_check = clusters_containing_patch($seed_i, $cluster_patches);
		if(count($clusters_to_check) > 1) {
			$seed_dists = get_seed_distances($seed_i, $clusters_to_check);
			//if the first cluster is same as current cluster then no need to remove
			reset($seed_dists);
			$first_key = key($seed_dists);
			if($cluster_i != $first_key) {
				$patches_to_remove[] = $seed_i;
			}
		}
	}
	
	if(count($patches_to_remove)) {
		$patches_to_keep = [];
		foreach($patches_i as $seed_i => $residues_i) {
			if(!in_array($seed_i, $patches_to_remove)) {
				$patches_to_keep[$seed_i] = $residues_i;
			}
		}
	}
	return $patches_to_keep;
	
}*/
Clusters.prototype._merge_unique_clusters = function (sets) {
    var stop = false;
    _(sets).forEach(function(set1Values, set1Key) {
        _(sets).forEach(function(set2Values, set2Key) {
            
            if(stop)
				return false;            
            
            if(set1Key == set2Key)
				return;
            
            let matches = _.uniq(_.intersection(set1Values, set2Values));
            if(matches.length > 0) {
                sets[set1Key] = _.uniq(_.merge(set1Values, set2Values)).sort();
                _.unset(sets, set2Key);
                stop = true;
                return false;
            }
            
        });
        if(stop) {
            return false;
        }
    });
    
	return sets;

}

Clusters.prototype._find_closest_core = function (residue, cores) {
	var seed_dists = this._get_seed_distances(residue, cores);
    return _.first(seed_dists).core;
}

Clusters.prototype._get_cluster_seeds = function (clusters) {
    var cluster_residues = [];
    _(clusters).forEach(function(clusterPatches, clusterKey) {
        _(clusterPatches).forEach(function(patchValues, patchKey) {
            cluster_residues.push(patchKey);
        });
	});
    cluster_residues = _.uniq(cluster_residues);
    return cluster_residues;
}

Clusters.prototype.find_clusters = function () {
	
	var self = this;

	//find residue densities
	var residue_densities = [];
	_(this.patches).forEach(function(patch) {
		residue_densities = _.concat(residue_densities, patch.residues);
	});
	residue_densities = _.countBy(residue_densities);

	var patch_densities = [];
	_(this.patches).forEach(function(patch) {
		patch_densities.push( self._calculate_patch_density(patch, residue_densities) );
	});
    
    //now some preliminary clusters
    var clusters = {};
    _(self.patches).forEach(function(patch1) {
        clusters[patch1.centroid] = [];
        _(self.patches).forEach(function(patch2) {
            if(patch1.centroid != patch2.centroid) {
                //if over 8 overlapping residues then it's a likely cluster                
                if(_.intersection(patch1.residues, patch2.residues).length > 8) {
                    clusters[patch1.centroid].push(patch2.centroid);
                }
            }
        });
	});

    //score the clusters
    var cluster_scores = [];
    _(clusters).forEach(function(centroids, clusterKey) {
        cluster_scores.push({centroid: clusterKey, score: self._cluster_score(centroids, patch_densities)});
    });
    var cluster_scores = _.reverse(_.sortBy(cluster_scores, 'score'));
    
    //discard duplicate and small clusters
    var discarded_clusters = [];
    _(cluster_scores).forEach(function(cluster1) {
        _(cluster_scores).forEach(function(cluster2 ) {
            
            if(cluster1.centroid != cluster2.centroid) {
                //if we have matching residues
                if(_.intersection(clusters[cluster1.centroid], clusters[cluster2.centroid]).length > 1) {
                    
                    if(cluster1.score > cluster2.score) {
                        discarded_clusters.push(cluster2.centroid);
                    } else if(cluster2.score > cluster1.score) {
                        discarded_clusters.push(cluster1.centroid);
                    }
                    
                }
            }
            
        });
	});
    
    discarded_clusters = _.uniq(discarded_clusters);
    
    var discarded_clusters = discarded_clusters.map(function(obj){ 
   		return parseInt(obj);
	});
    
    //final clusters
    var final_clusters = [];
    _(clusters).forEach(function(centroids, clusterKey) {
        if(discarded_clusters.indexOf(parseInt(clusterKey)) == -1) {
            final_clusters.push(clusterKey);
        }
    });
    
    //find core residues
    var patch_cores = {};
    var core_residues = [];
    _(self.patches).forEach(function(patch) {
        let pd = _.find(patch_densities, { 'centroid': patch.centroid })
        patch_cores[patch.centroid] = self._find_core(pd.density);
        core_residues = _.concat(core_residues, patch_cores[patch.centroid]);        
    });
    core_residues = _.uniq(core_residues);
    
    var cluster_uniques = {};
    _(self.patches).forEach(function(patch) {
        let intersection = _.intersection(core_residues, patch.residues).sort().join(",");
        if(!_.isArray(cluster_uniques[intersection])) {
            cluster_uniques[intersection] = [];
        }
        cluster_uniques[intersection].push(patch.centroid);
    });
    
    //now we get the largest unique clusters
    var discarded_clusters = [];
    _(cluster_uniques).forEach(function(seeds1, cluster1Key) {
        _(cluster_uniques).forEach(function(seeds2, cluster2Key) {
            
            if(cluster1Key == cluster2Key)
                return;

            let cluster1Data = cluster1Key.split(",");
            let cluster2Data = cluster2Key.split(",");
            
            let matches = _.uniq(_.intersection(cluster1Data, cluster2Data));

            if(matches.length > 0) {
                if(seeds1.length > seeds2.length) {
                    discarded_clusters.push(cluster2Key);
                }

                if(seeds2.length > seeds1.length) {
                    discarded_clusters.push(cluster1Key);
                }
            }
            
            
        });
    });
    
    discarded_clusters = _.uniq(discarded_clusters);
    var uniques_clusters = _.difference(_.keys(cluster_uniques), discarded_clusters);
    
    var uniqueClusters = {};
    _(uniques_clusters).forEach(function(unique_cluster) {
        uniqueClusters[unique_cluster] = unique_cluster.split(',');
    });
    
    var starting_sets = -1;
    while (starting_sets != uniqueClusters.length) {
        starting_sets = uniqueClusters.length
        uniqueClusters = self._merge_unique_clusters(uniqueClusters);
    }
    
    var final_clusters = {};
    _(uniqueClusters).forEach(function(mainCluster, mainKey) {
        final_clusters[mainKey] = [];
        _(cluster_uniques).forEach(function(patchSeeds, clusterKey) {
            let clusterSeeds = clusterKey.split(',');
            let matches = _.uniq(_.intersection(clusterSeeds, mainCluster));
            if(matches.length > 0) {
                final_clusters[mainKey] = _.concat(patchSeeds, final_clusters[mainKey]);
            }
        });
        final_clusters[mainKey] = _.uniq(final_clusters[mainKey]);
        final_clusters[mainKey] = _.sortBy( final_clusters[mainKey], function( val ){ return val; } );
    });
	console.log('final_clusters', final_clusters);
    
    var cluster_patches = {};
    var cluster_residues = {};
    _(final_clusters).forEach(function(patch_seeds, cluster_key) {
        cluster_residues[cluster_key] = [];
        
        _(patch_seeds).forEach(function(patch_seed) {
            let tmp = {};
            tmp.centroid = cluster_key;
            tmp.patch_seed = patch_seed;
            tmp.residues = _.find(self.patches, {seed: parseInt(tmp.patch_seed)}).residues;
            if(_.isUndefined(cluster_patches[tmp.centroid]))
                cluster_patches[tmp.centroid] = {};
            if(_.isUndefined(cluster_patches[tmp.centroid][tmp.patch_seed]))
                cluster_patches[tmp.centroid][tmp.patch_seed] = tmp.residues;
        });
        
    });
    
    //now remove overlapping patches
    var clean_clusters = {};
    _(cluster_patches).forEach(function(patchResidues, centroid) {
        clean_clusters[centroid] = self._remove_overlapping_patches(centroid, cluster_patches);
    });
    
    //now remove clusters with only 1 patch
    var single_clusters = [];
    _(cluster_patches).forEach(function(patches, clusterKey) {
        if(cluster_patches[clusterKey].length == 1) {
            single_clusters[clusterKey] = patches;
        }
    });
    
    var cluster_residues = self._get_cluster_seeds(cluster_patches);
    cluster_residues = _.uniq(cluster_residues);
    
    var cluster_groups = {};
    _(cluster_residues).forEach(function(residue) {
        let closest_core = self._find_closest_core(residue, _.keys(cluster_patches));
        if(_.isUndefined(cluster_groups[closest_core]))
            cluster_groups[closest_core] = [];
        cluster_groups[closest_core].push(parseInt(residue));
    });
    
    return cluster_groups;
}
