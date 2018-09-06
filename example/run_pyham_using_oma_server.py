#  This is the HAM package
import pyham

#  OPTIONAL: only if you want to have the logger information printed
import logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s %(name)-12s %(levelname)-8s %(message)s")

# Initialise pyHam with a gene oma id

query = 'HUMAN12'
pyham_analysis = pyham.Ham(query_database=query, use_data_from='oma')


hog = pyham_analysis.get_hog_by_id('HOG:0359282')

# create the iHam for it and store it into an html file
output_name = "iHam{}.html".format(hog.hog_id)
pyham_analysis.create_iHam(hog=hog, outfile=output_name)

# create the iHam for it and store it into an html file
output_name = "TreeProfile{}.html".format(hog.hog_id)
pyham_analysis.create_tree_profile(hog=hog, outfile=output_name)
