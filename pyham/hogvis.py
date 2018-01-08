from __future__ import absolute_import
from __future__ import unicode_literals
from __future__ import print_function
from __future__ import division
from builtins import int
from builtins import range
from future import standard_library
standard_library.install_aliases()
import json
import re
from string import Template
import collections
import logging

logger = logging.getLogger(__name__)

class Hogvis(object):

    """
    Not documented yet.
    
    """
    def __init__(self, newick_str, hog):
        #print("THE HOGVIS IS AT A EXPERIMENTAL STAGE. PLEASE USE IT CAREFULLY.")
        self.hog = hog
        self.newick_str = newick_str
        self.ogs_mapper = self._get_groups_to_level_mapper()
        self.per_species = self._get_per_species_structure()
        self.html_template = self._get_html_template()
        self.xrefs = {gene.unique_id: gene.get_dict_xref() for gene in hog.get_all_descendant_genes()}

        self.renderHTML = self.html_template.safe_substitute({'name': hog.hog_id,
                                                             'species_tree': self.newick_str,
                                                             'xrefs': json.dumps(self.xrefs),
                                                             'per_species': json.dumps(self.per_species)})

    def _get_html_template(self):
            return Template("""<!DOCTYPE html>
                    <html>
                    <head>
                     <title>$name</title>
                     <!-- D3 -->
                     <script src="http://d3js.org/d3.v3.min.js" charset="utf-8"></script>

                     <!-- Underscore -->
                     <script src="http://cdnjs.cloudflare.com/ajax/libs/underscore.js/1.6.0/underscore-min.js"></script>

                     <!-- TnT -->
                     <link rel="stylesheet" href="http://omabrowser.org/static/css/tnt.css" type="text/css" />
                     <script src="http://omabrowser.org/static/js/tnt_old.min.js"></script>
                     <script src="http://omabrowser.org/static/js/hog.js.05869ae311cc8b8e04271cdb1a904fd4"></script>

                      <h1>HOG Viewer for $name</h1>
                      <div id="hog_tree"></div>

                      <script>

                        (function () {
                          var tree = tnt.tree.parse_newick('$species_tree');
                            var query;
                            var per_species = $per_species ;
                            var tooltip_data = $xrefs ;
                            var options = {'show_internal_labels': "true",
                                           'oma_info_url_template': ""};

                            var viz = tnt();
                            var theme = hog_theme();
                            theme (viz, document.getElementById("hog_tree"), query, per_species, tree, tooltip_data, options);
                        }) ();

                      </script>
                    </body>
                    </html>""")

    def _get_groups_to_level_mapper(self):
            hogs = self.hog.get_all_descendant_hogs()
            return _OGLevelMapper(hogs)

    def _get_per_species_structure(self):
            per_species = {}
            re_non_char = re.compile(r'\W')

            for org, genelist in self.hog.get_all_descendant_genes_clustered_by_species().items():
                org = org.taxon.name
                per_species[org] = cur_map = {}
                cur_map[org] = [[] for _ in range(len(genelist))]
                for pos, gene in enumerate(genelist):
                    cur_map[org][pos].append(int(gene.unique_id))
                    parent = gene.parent
                    while parent != self.hog.parent:
                        lev = parent.genome.name
                        lev_repr = lev
                        if re_non_char.search(lev):
                            lev_repr = '{}'.format(lev)
                        if lev_repr not in cur_map:
                            cur_map[lev_repr] = [[] for _ in range(self.ogs_mapper.nr_subhogs_on_level(lev))]
                        if lev != org:
                            cur_map[lev_repr][self.ogs_mapper.index(parent)].append(int(gene.unique_id))
                        parent = parent.parent

            for fspe, dsp in per_species.items():
                for scdsp, lg in dsp.items():
                    for arr in lg:
                        x = lg.index(arr)
                        per_species[fspe][scdsp][x] = list(set(arr))

            return per_species

class _OGLevelMapper(object):
    def __init__(self, ogs):
        self.levels = collections.defaultdict(list)
        self.id2pos = {}
        for og in ogs:
            lev = og.genome.name
            self.levels[lev].append(og)
            pos = len(self.levels[lev]) - 1
            try:
                pos_before = self.id2pos[id(og)]
                if pos != pos_before:
                    logger.warn(
                        "HOG {} with several levels has inconsistent positions."
                        "Lev: {}, pos_before, pos: {}/{}".format(id(og), lev, pos_before, pos))
                    for _ in range(pos_before-pos):
                        self.levels[lev].insert(pos, None)
            except KeyError:
                self.id2pos[id(og)] = pos

    def index(self, og):
        return self.id2pos[id(og)]

    def nr_subhogs_on_level(self, level):
        return len(self.levels[level])
