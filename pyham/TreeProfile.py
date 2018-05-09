from __future__ import absolute_import
from __future__ import unicode_literals
from __future__ import print_function
from __future__ import division
from builtins import str
from future import standard_library
standard_library.install_aliases()
from .abstractgene import HOG
import logging
import json

logger = logging.getLogger(__name__)


class TreeProfile(object):
    """
    Object that map on each node of the Ham species tree the related evolutionary information such as the number of 
    genes the related genome contains or the number of evolutionary events (duplications, genes loss or gain of genes)
    that occurs between the parent node and itself. This can be either applied to the full Ham comparative genomic setup
    or for a specific HOG.

    Attributes:
        | pyham (:obj:`pyham.pyham.Ham`): :obj:`set` of HOG ids used by the FilterOrthoXMLParser.
        | hog (:obj:`pyham.abstractgene.HOG`, optional): If specified, compute TreeProfile on a single HOGs. Defaults is None.
        | treemap (:obj:`ete3.Etree`): Ete3 Etree object containing the taxonomy of interest with annotated node.
    """

    def __init__(self, ham, hog=None):

        self.ham = ham
        self.hog = hog

        if hog is None:
            self.treemap = self.compute_tree_profile_full()

        elif isinstance(hog, HOG):
            self.treemap = self.computeTP_hog(hog)
        else:
            raise TypeError("Invalid argument {} for HOG".format(hog))

    def computeTP_hog(self, hog):
        """  
        Create the treeMap object (ete3 Etree object with annotated nodes) that contained the following features for
        each nodes (representing an AbstractGene):
            - nbr_genes: numbers of HOG/Gene the genome contains. For leaves, it also counts the singletons.
            - dupl: numbers of HOGs that have arose by duplication in between this node and its parent.
            - identical: numbers of HOGs that have stay identical (in term of copy numbers) in between this node and
            its parent.
            - lost: numbers of HOGs that have been lost in between this node and its parent.

        Returns:
            TreeMap
        """

        # copy the required taxonomy using query hog as root level
        treeMap = hog.genome.taxon.copy(method="newick")

        # create a dictionary that map node with related hogs/genes
        levelGroups = {}
        for lvl in treeMap.traverse():
            levelGroups[lvl.name]=[]

        # add all of subhog to the related level in levelGroups
        for subhog in hog.get_all_descendant_hogs():
            levelGroups[subhog.genome.name].append(subhog)

        # add empty extant genome to levelGroups
        for extantGenome in treeMap.get_leaves():
            levelGroups[extantGenome.name] = []

        # add genes to related extant genome in levelGroups
        for species, genes in hog.get_all_descendant_genes_clustered_by_species().items():
            levelGroups[species.name] = genes

        for lvl in treeMap.traverse():
            lvl.add_feature("nbr_genes", len(levelGroups[lvl.name]))
            lvl.add_feature("gain", None)

            if not lvl.is_root():
                cpt_dupl = 0
                cpt_ident = 0
                for h in levelGroups[lvl.name]:
                    if h.arose_by_duplication != False:
                        cpt_dupl += 1
                    else:
                        cpt_ident += 1

                cpt_lost = 0
                for h_up in levelGroups[lvl.up.name]:
                    if len(h_up.children) != 0:
                        found = False
                        for child in h_up.children:
                            if child.genome.taxon.name == lvl.name:
                                found = True
                        if found is False:
                            cpt_lost += 1
                    else:
                        cpt_lost += 1

            else:
                cpt_dupl = None
                cpt_ident = None
                cpt_lost = None

            lvl.add_feature("dupl", cpt_dupl)
            lvl.add_feature("lost", cpt_lost)
            lvl.add_feature("identical", cpt_ident)

        return treeMap

    def compute_tree_profile_full(self):
        """  
        Create the treeMap object (ete3 Etree object with annotated nodes) that contained the following features for
        each nodes (the root node contains only nbr_genes):
            - nbr_genes: numbers of HOG/Gene the genome contains. For leaves, it also counts the singletons.
            - dupl: numbers of HOGs that have arose by duplication in between this node and its parent.
            - lost: numbers of HOGs that have been lost in between this node and its parent.
            - gain: numbers of HOGs that have "emerged" at this node.
            - identical: numbers of HOGs that have stay identical (in term of copy numbers) in between this node and
            its parent.
            
        In order to get all those node informations, the HOGMap between each node and its parent is computed.
            
        Returns:
            TreeMap
        """

        def _add_annot(node, nbr, dupl, lost, gain, identical):
            node.add_feature("nbr_genes", nbr)
            node.add_feature("dupl", dupl)
            node.add_feature("lost", lost)
            node.add_feature("gain", gain)
            node.add_feature("identical", identical)

        treeMap = self.ham.taxonomy.tree.copy(method="newick")

        for node in treeMap.traverse():
            if node.is_root():
                node_genome = self.ham.get_ancestral_genome_by_name(node.name)
                _add_annot(node, len(node_genome.genes), None, None, None, None)

            else:
                node_genome_up = self.ham.get_ancestral_genome_by_name(node.up.name)

                if node.is_leaf():
                    node_genome = self.ham._get_extant_genome_by_name(name=node.name)
                    nbr = node_genome.get_number_genes(singleton=True)

                else:
                    node_genome = self.ham.get_ancestral_genome_by_name(node.name)
                    nbr = node_genome.get_number_genes()
                hogmap = self.ham._get_HOGMap({node_genome, node_genome_up})

                nbr_duplicate = 0
                for gs in hogmap.DUPLICATE.values():
                    nbr_duplicate += len(gs)
                _add_annot(node, nbr, nbr_duplicate, len(hogmap.LOSS), len(hogmap.GAIN), len(hogmap.IDENTICAL.keys()))

        return treeMap

    def export(self, output, layout_function=None, display_internal_histogram=True):

        """  
        Method to export the tree profile object as figure (available format .SVG, .PDF, .PNG).
        
        -- Some magic going there --
        
        Args:
            | output (:obj:`str`): output file name. The extension of output will set the format of the figure (SVG, .PDF, .PNG)
            | layout_function (:obj:`function`, optional): custom layout_fn for ete3 TreeStyle.
            | display_internal_histogram (:obj:`Boolean`, optional): Display internal node as histogram or raw text with numbers. Defaults to True.
        """

        from ete3 import TreeStyle, TextFace, NodeStyle, BarChartFace

        # maximum number of genes per node in this treeMap
        max_genes = max([d for d in self.treemap.traverse()], key=lambda x:x.nbr_genes).nbr_genes

        if self.hog is None:
            _color_scheme = ["#41c1c2","#bdc3c7","#f39c12","#27ae60","#e74c3c"]
            _label_legend = ["Genes","Identicals","Duplicated","Novel","Lost"]
            _values_legend = [max_genes,max_genes,max_genes,max_genes,max_genes]
            w_legend = 50 # todo calculate base on len(_values)
        else:
            _color_scheme = ["#41c1c2", "#bdc3c7", "#f39c12", "#e74c3c"]
            _label_legend = ["Genes", "Identicals", "Duplicated", "Lost"]
            _values_legend = [max_genes, max_genes, max_genes, max_genes]
            w_legend = 40  # todo calculate base on len(_values)

        def _layout(node):

            if self.hog is None:
                 _label = [str(node.nbr_genes),str(node.identical),str(node.dupl),str(node.gain),str(node.lost)]
            else:
                _label = [str(node.nbr_genes), str(node.identical), str(node.dupl), str(node.lost)]

            def _add_face(name_feature, value_feature, cnum=1, pos="branch-right"):
                node.add_face(TextFace("{}: {}".format(name_feature, value_feature)), column=cnum, position=pos)

            def _add_faces(cNbr=1, posNbr="branch-right", cAttr=2, posAtt="branch-right"):

                _add_face("#genes", node.nbr_genes, cnum=cNbr, pos=posNbr)

                if node.identical is not None:
                    _add_face("#Identical", node.identical, cnum=cAttr, pos=posAtt)

                if node.dupl is not None:
                    _add_face("#Duplicated", node.dupl, cnum=cAttr, pos=posAtt)

                if node.gain is not None:
                    _add_face("#Novel", node.gain, cnum=cAttr, pos=posAtt)

                if node.lost is not None:
                    _add_face("#Lost", node.lost, cnum=cAttr, pos=posAtt)

            if node.is_leaf():
                if display_internal_histogram:
                    if self.hog is None:
                        values = [node.nbr_genes,node.identical,node.dupl,node.gain,node.lost]
                        w_plot = 50
                    else:
                        values = [node.nbr_genes, node.identical, node.dupl, node.lost]
                        w_plot = 40
                    node.add_face(BarChartFace(values, deviations=None, width=w_plot, height=25, colors=_color_scheme, labels=_label, min_value=0, max_value=max_genes, label_fsize=6, scale_fsize=6),column=1, position = "branch-right")
                else:
                    _add_faces()

            else:

                if display_internal_histogram:
                    if node.is_root():
                        node.add_face(BarChartFace([node.nbr_genes], deviations=None, width=10, height=25, colors=["#41c1c2"], labels=[str(node.nbr_genes)], min_value=0, max_value=max_genes, label_fsize=6, scale_fsize=6),column=0, position = "branch-bottom")
                    else:
                        if self.hog is None:
                            values = [node.nbr_genes,node.identical,node.dupl,node.gain,node.lost]
                            w_plot = 50
                        else:
                            values = [node.nbr_genes,node.identical,node.dupl,node.lost]
                            w_plot = 40
                        node.add_face(BarChartFace(values, deviations=None, width=w_plot, height=25, colors=_color_scheme, labels=_label, min_value=0, max_value=max_genes, label_fsize=6, scale_fsize=6),column=1, position = "branch-top")



                else:
                    _add_faces(cNbr=0, posNbr="branch-top", cAttr=0, posAtt="branch-bottom")

        ts = TreeStyle()

        if layout_function is not None:
            ts.layout_fn = layout_function
        else:
            ts.layout_fn = _layout
            ts.legend.add_face(BarChartFace(_values_legend, deviations=None, width=w_legend, height=25, colors=_color_scheme, labels=_label_legend, min_value=0, max_value=max_genes, label_fsize=6, scale_fsize=6),column=0)
            ts.legend_position = 3

        self.treemap.render(output,tree_style=ts)


    def export_as_html(self, output ):

        """  
        Method to export the tree profile object as an interactive tool embedded into a html file.


        Args:
            | output (:obj:`str`): output file name. If filename finish by .html, a double click load the file into your default browser.
        """

        def visit_custom(node, data):

            current = {
                "name": node.name,
                "numberGenes": node.nbr_genes,
                "length": 0.01,
                "collapsed": "false",
                "evolutionaryEvents": {
                    "identical": None,
                    "duplicated": None,
                    "gained": None,
                    "lost": None
                }
            }

            if node.is_root():
                current['evolutionaryEvents'] = False

            else:
                current['evolutionaryEvents']["identical"] = node.identical
                current['evolutionaryEvents']["duplicated"] = node.dupl
                current['evolutionaryEvents']["gained"] = node.gain
                current['evolutionaryEvents']["lost"] = node.lost

            if not node.is_leaf():
                current['children'] = []

            for child in node.children:
                current['children'].append(visit_custom(child, data))

            return current

        template = '''
        <!DOCTYPE html>
        <html>
        <head>
            <title>Phylo.io</title>
            <meta charset="UTF-8">
            <script src="https://peterolson.github.com/BigInteger.js/BigInteger.min.js"></script>
            <script type="text/javascript" src="https://cdn.rawgit.com/DessimozLab/phylo-io/7f83193b4efaee7920451ecfb4d917014a4df437/www//js/jquery-2.1.4.min.js"></script>
            <script type="text/javascript" src="https://cdn.rawgit.com/DessimozLab/phylo-io/7f83193b4efaee7920451ecfb4d917014a4df437/www/js/treecompare.js"></script>
            <script type="text/javascript" src="http://underscorejs.org/underscore-min.js"></script>
            <script type="text/javascript" src="https://cdn.rawgit.com/DessimozLab/phylo-io/7f83193b4efaee7920451ecfb4d917014a4df437/www//js/d3.min.js"></script>
            <script type="text/javascript" src="https://cdn.rawgit.com/DessimozLab/phylo-io/7f83193b4efaee7920451ecfb4d917014a4df437/www//js/bootstrap.min.js"></script>
            <script type="text/javascript" src="https://cdn.rawgit.com/DessimozLab/phylo-io/7f83193b4efaee7920451ecfb4d917014a4df437/www//js/FileSaver.min.js"></script>
            <link rel="stylesheet" href="https://use.fontawesome.com/releases/v5.0.8/css/solid.css" integrity="sha384-v2Tw72dyUXeU3y4aM2Y0tBJQkGfplr39mxZqlTBDUZAb9BGoC40+rdFCG0m10lXk" crossorigin="anonymous">
            <link rel="stylesheet" href="https://use.fontawesome.com/releases/v5.0.8/css/fontawesome.css" integrity="sha384-q3jl8XQu1OpdLgGFvNRnPdj5VIlCvgsDQTQB6owSOHWlAurxul7f+JpUOVdAiJ5P" crossorigin="anonymous">
            <link rel="stylesheet" type="text/css" href="https://cdn.rawgit.com/DessimozLab/phylo-io/7f83193b4efaee7920451ecfb4d917014a4df437/www//css/bootstrap.min.css">
            <link rel="stylesheet" type="text/css" href="https://cdn.rawgit.com/DessimozLab/phylo-io/7f83193b4efaee7920451ecfb4d917014a4df437/www//css/bootstrap-theme.min.css">
            <link rel="stylesheet" type="text/css" href="https://cdn.rawgit.com/DessimozLab/phylo-io/7f83193b4efaee7920451ecfb4d917014a4df437/www//css/style.css">
            <style> text {{stroke: none;}}</style>
        </head>
        <body>
        <div id="vis-container1" style="width: 100%; height: 100%;"></div>
        <script type="text/javascript">
            var treecomp = TreeCompare().init({{
                scaleColor: "black",
                showHistogramValues: true,
                showHistogramSummaryValue: false
            }});
            treeData = '{json_data}';
            var tree1 = treecomp.addTree(treeData, undefined, "single");
            treecomp.viewTree(tree1.name, "vis-container1");
            treecomp.addMainLegend(tree1.name);
        </script>
        </body>
        </html>
        '''

        data = {}

        data = visit_custom(self.treemap, data)

        data = json.dumps(data)

        template = template.format(json_data=data)

        with open(output, 'w') as outfile:
            outfile.write(template)
