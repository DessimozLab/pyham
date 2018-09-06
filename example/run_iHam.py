#  This is the HAM package
import pyham

#  OPTIONAL: only if you want to have the logger information printed
import logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s %(name)-12s %(levelname)-8s %(message)s")

# Initialise pyHam with a phyloxml tree and orthoXML HOGs
phyloxml_path = "simpleEx.phyloxml"
orthoxml_path3 = "simpleEx.orthoxml"

pyham_analysis = pyham.Ham(phyloxml_path, orthoxml_path3, use_internal_name=True, tree_format='phyloxml')


hog = pyham_analysis.get_hog_by_id('HOG:0355161')

# create the iHam for it and store it into an html file
output_name = "iHam{}.html".format(hog.hog_id)
pyham_analysis.create_iHam(hog=hog,outfile=output_name)
