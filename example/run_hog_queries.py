# -*- coding: utf-8 -*-
import pyham

# Initialise pyHam with a phyloxml tree and orthoXML HOGs
phyloxml_path = "simpleEx.phyloxml"
orthoxml_path3 = "simpleEx.orthoxml"

pyham_analysis = pyham.Ham(phyloxml_path, orthoxml_path3, use_internal_name=True, tree_format='phyloxml')

# Get the hog of interest
hog = pyham_analysis.get_hog_by_id('HOG:0355161')

question1 = '\n 1 - What is the root level of this HOG, i.e. at which ancestral taxonomic level did this gene originate?'
print(question1)
print('\t' +  hog.genome.name)

question2 = '\n 2 - How many genes are in this gene family?'
print(question2)

hog_genes = hog.get_all_descendant_genes()
print('\t' +  str(len(hog_genes)))

question3 = '\n 3 - Which extant genomes have the most copies of this gene?'
print(question3)

hog_genes_clustered = hog.get_all_descendant_genes_clustered_by_species()
for species, genes in hog_genes_clustered.items():
    print('\t' +  species.name + " -> " +  str(len(genes)))

question4 = '\n 4 - How many genes are in this gene family in humans?'
print(question4)

human_genome = pyham_analysis.get_extant_genome_by_name('Homo sapiens')
print('\t' + str(len(hog_genes_clustered[human_genome])) )

question5 = '\n 5 - How many genes are in the HOG at the Primates level?'
print(question5)

descendant_hog = hog.get_all_descendant_hogs()
primates_genome = pyham_analysis.get_ancestral_genome_by_name('Primates')

number_hog_at_primates = 0

for sub_hog_at_primates in descendant_hog:
    if sub_hog_at_primates.genome == primates_genome:
        number_hog_at_primates +=1

print('\t' + str(number_hog_at_primates) )


