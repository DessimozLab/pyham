from __future__ import absolute_import
from __future__ import unicode_literals
from __future__ import print_function
from __future__ import division
from future import standard_library
standard_library.install_aliases()
from string import Template
import logging
import json
import datetime
import lxml.etree as etree

logger = logging.getLogger(__name__)

class IHAM(object):

    def __init__(self, newick_str, hog):
        self.hog = hog
        self.newick_str = newick_str
        self.html_template = self._get_html_template()
        self.famdata = json.dumps(self._get_famdata())
        self.orthoxml = OrthoXML_manager(self.hog)
        self.HTML = self.html_template.safe_substitute({'name': hog.hog_id,
                                                        "tree": self.newick_str,
                                                        "orthoxml": self.orthoxml.get_orthoxml_str(),
                                                        "fam_data": self.famdata
                                                        })

    def _get_html_template(self):
        x =Template('''<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <title>HOG $name</title>
    <script src="https://d3js.org/d3.v3.js"></script>

    <!-- Bootstrap CDN -->
    <link rel="stylesheet" href="https://stackpath.bootstrapcdn.com/bootstrap/4.1.0/css/bootstrap.min.css"
          integrity="sha384-9gVQ4dYFwwWSjIDZnLEWnxCjeSWFphJiwGPXr1jddIhOegiu1FwO5qRGvFXOdJZ4" crossorigin="anonymous">
    <script src="https://code.jquery.com/jquery-3.3.1.slim.min.js"
            integrity="sha384-q8i/X+965DzO0rT7abK41JStQIAqVgRVzpbzo5smXKp4YfRvH+8abtTE1Pi6jizo"
            crossorigin="anonymous"></script>
    <script src="https://cdnjs.cloudflare.com/ajax/libs/popper.js/1.14.0/umd/popper.min.js"
            integrity="sha384-cs/chFZiN24E4KMATLdqdvsezGxaGsi4hLGOzlXwp5UZB1LY//20VyM2taTB4QvJ"
            crossorigin="anonymous"></script>
    <script src="https://stackpath.bootstrapcdn.com/bootstrap/4.1.0/js/bootstrap.min.js"
            integrity="sha384-uefMccjFJAIv6A+rW+L4AHf99KvxDjWSu1z9VI8SKNVmz4sk7buKt/6v9KI65qnm"
            crossorigin="anonymous"></script>

    <!-- Roboto font -->
    <link href="https://fonts.googleapis.com/css?family=Roboto" rel="stylesheet">

    <!-- TnT -->
    <link rel="stylesheet" href="https://tntvis.github.io/tnt/build/tnt.css" type="text/css"/>
    <script src="https://tntvis.github.io/tnt/build/tnt.js" charset="utf-8"></script>

    <!-- TnT Tooltip-->
    <link rel="stylesheet" href="https://tntvis.github.io/tnt.tooltip/build/tnt.tooltip.css" type="text/css"/>
    <script src="https://omabrowser.org/static/js/tnt.tooltip.min.js" charset="utf-8"></script>

    <script src="https://dessimozlab.github.io/iHam/iHam.js"></script>
    <link rel="stylesheet" href="https://dessimozlab.github.io/iHam/iHam.css" type="text/css"/>

    <style>
        body {
            font-family: 'Roboto', sans-serif;
        }

        #updating {
            position: absolute;
            left: 20px;
            display: none;
        }

        .alert_remove {
            margin-bottom: 0px;
            padding: 4px;
            display: none;
        }

        .alert-link {
            cursor: pointer;
        }

        #header {
            margin-left: 20px;
            margin-bottom: 20px;
        }

        #header > h3 {
            margin-top: 20px;
            margin-bottom: 20px;
        }

        #menu-bar > div {
            display: inline-block;
        }

        .dropdown-toggle {
            padding-top: 7px;
            padding-bottom: 7px;
        }
    </style>
</head>
<body>

<div id="header">
    <h3>Hierarchical group HOG:0000474 open at level of <span id="current-node"></span></h3>

    <div id="menu-bar">
        
        <div id="gene-tooltips-dropdown" class="dropdown">
            <button class="btn btn-sm btn-outline-dark dropdown-toggle" type="button" data-toggle="dropdown"
                    aria-haspopup="true" aria-expanded="false">
                Show gene tooltips on
            </button>
            <div class="dropdown-menu" aria-labelledby="dropdownMenuButton">
                <a class="dropdown-item active" href="#">Click</a>
                <a class="dropdown-item" href="#">Mouseover</a>
            </div>
        </div>

        <div id="percentage-coverage-selector">
            <button class="btn btn-sm btn-outline-dark" type="button"
                    aria-haspopup="true" aria-expanded="false">Remove columns under <input id="set_min_coverage"
                                                                                           type="number" step="10"
                                                                                           value="0"
                                                                                           min="0" max="100">% of
                species coverage
            </button>
        </div>
    </div>
</div>

<div id="updating">
    Updating...
</div>

<div class="alert alert-info text-center alert_remove"
     role="alert"
>
    Lowly supported hogs have been removed as per settings.
    <a class="alert-link" id="reset_column_filter">Click here to reset.</a>
</div>


<div style="width: 1500px; min-width: 500px;" id="iham"></div>

<script>

    (function (div) {

        const data = {
            "tree": '$tree',
            "orthoxml": `$orthoxml`,
            "fam_data": $fam_data
        }


        var theme = iHam()
                .on('node_selected', function (node) {
                    d3.select('#current-node')
                            .text(node.node_name());
                })
                .orthoxml(data.orthoxml)
                .show_oma_link(false)
                .newick(data.tree)
                .fam_data(data.fam_data)
                .tree_width(330)
                .board_width(530)
                .query_gene({id: 3965})
                .on("updating", function () {
                    d3.select("#updating")
                            .style("display", 'block');
                })
                .on("updated", function () {
                    d3.select("#updating")
                            .style("display", "none");
                })
                .on("hogs_removed", function (what) {
                    if (what.length) {
                        d3.select(".alert_remove")
                                .style("display", "block")
                    } else {
                        d3.select(".alert_remove")
                                .style("display", "none");
                    }
                });

        theme(div);

        // Update the color schemas

        // Update event for gene tooltips
        d3.select("#gene-tooltips-dropdown")
                .selectAll("a")
                .on("click", function () {
                    // Manage state of menu itself
                    d3.select(this.parentNode).selectAll("a").classed("active", false);
                    d3.select(this)
                            .classed("active", true);

                    if (d3.select(this).text() === "Click") {
                        theme.gene_tooltips_on("click");
                    }
                    if (d3.select(this).text() === "Mouseover") {
                        theme.gene_tooltips_on("mouseover");
                    }
                });

        // Set minimum species coverage
        d3.select("#percentage-coverage-selector").select("input")
                .on("input", function () {
                    theme.coverage_threshold(d3.select(this).property("value"));
                });

        // Reset the coverage
        d3.select("#reset_column_filter")
                .on("click", function () {
                    d3.select("#percentage-coverage-selector").select("input").property("value", 0);
                    theme.coverage_threshold(0);
                })
    })(document.getElementById('iham'));


</script>
</body>
</html>''')
        return x

    def _get_famdata(self):

        fd = []

        for gene in self.hog.get_all_descendant_genes():
            data = {}
            data["taxon"] = {"species": gene.genome.name}
            data["xrefid"] = gene.prot_id
            data["sequence_length"] = 0
            data["gc_content"] = 0
            data["protid"] = gene.prot_id
            data["id"] = gene.unique_id

            fd.append(data)

        return fd

class OrthoXML_manager(object):

    def __init__(self, hog):

        self.hog = hog

        origin = "OMA"
        vers = datetime.date.today().strftime("%b %Y")
        origin += " standalone (bottom-up)"

        self.xml = self._create_xml(origin, vers)

        self._add_species_data()

        self._add_groups()

    def _create_xml(self, origin, originVersion):

        xml_core = etree.Element("orthoXML")
        xml_core.set("origin", origin)
        xml_core.set("originVersion", originVersion)
        xml_core.set("version", "0.3")
        xml_core.set("xmlns", "http://orthoXML.org/2011/")
        return xml_core

    def _add_species_data(self):

        sp_data = self.hog.get_all_descendant_genes_clustered_by_species()

        for species, childs in sp_data.items():

            species_xml = etree.Element("species")
            species_xml.set("name", species.name)
            species_xml.set("NCBITaxId", str(-1))
            self.xml.insert(0, species_xml)

            # Add <database> into <species>
            database_xml = etree.SubElement(species_xml, "database")
            database_xml.set("name", "None")
            database_xml.set("version", "None")

            # Add <genes> TAG into <database>
            genes_xml = etree.SubElement(database_xml, "genes")

            # Fill <genes> with <gene>
            for gene_obj in childs:
                gene_xml = etree.SubElement(genes_xml, "gene")
                gene_xml.set("id", str(gene_obj.unique_id))
                gene_xml.set("protId", str(gene_obj.prot_id))

    def _add_groups(self):

        from . import abstractgene as absGene

        def _process_child(child, current_xml):

            if isinstance(child, absGene.Gene):
                generef_xml = etree.SubElement(current_xml, "geneRef")
                generef_xml.set('id', str(child.unique_id))
            else:
                _visit(child, current_xml)

        def _visit(hog, parent):

            if len(hog.children) == 1:
                current_hog_xml = parent

            elif len(hog.duplications) >= 1:

                dup_child = []

                for duplicationNode in hog.duplications:
                    dup_child += duplicationNode.children

                remaining_hog = list(set(hog.children) - set(dup_child))

                if len(remaining_hog) == 0:

                    current_hog_xml = parent

                else:

                    current_hog_xml = etree.SubElement(parent, "orthologGroup")
                    taxon = etree.SubElement(current_hog_xml, "property")
                    taxon.set("name", 'TaxRange')
                    taxon.set("value", hog.genome.name)

                    current_hog_xml.set('id', str(hog.hog_id))

            else:
                current_hog_xml = etree.SubElement(parent, "orthologGroup")
                taxon = etree.SubElement(current_hog_xml, "property")
                taxon.set("name", 'TaxRange')
                taxon.set("value", hog.genome.name)

                current_hog_xml.set('id', str(hog.hog_id))

            processed_child = []
            if len(hog.duplications) > 0:

                for duplicationNode in hog.duplications:

                    paralogGroup = etree.SubElement(current_hog_xml, "paralogGroup")

                    for child in duplicationNode.children:
                        _process_child(child, paralogGroup)
                        processed_child.append(child)

            remaining_hog = list(set(hog.children) - set(processed_child))
            for child in remaining_hog:
                    _process_child(child, current_hog_xml)

        self.groupsxml = etree.SubElement(self.xml, "groups")

        _visit(self.hog, self.groupsxml)

    def get_orthoxml_str(self):
        return etree.tostring(self.xml, encoding='utf8', method='xml', pretty_print=False).decode('utf-8')







