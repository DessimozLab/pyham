import unittest
from pyham import ham, utils, iham
import os, tempfile

def generate_hog_pair(newick_name, orthoxml_name):

    ph, hog, iham_str = get_pyham_data_for_test(newick_name, orthoxml_name)

    file = tempfile.NamedTemporaryFile(mode='w', delete=False)

    with file as f: f.write(iham_str)

    ph2, hog2, iham_str2 = get_pyham_data_for_test(newick_name, file.name, orthoxml_file=True)


    return hog, hog2

def get_pyham_data_for_test(nwk, orthoxml, orthoxml_file = False):

    nwk_path = os.path.join(os.path.dirname(__file__), nwk)
    tree_str = utils.get_newick_string(nwk_path, type="nwk")

    orthoxml_path = orthoxml if orthoxml_file else os.path.join(os.path.dirname(__file__), orthoxml)

    ph = ham.Ham(tree_file=tree_str, hog_file=orthoxml_path, type_hog_file='orthoxml',
                      use_internal_name=True)

    hog = ph.get_list_top_level_hogs()[0]

    iham_str  = iham.OrthoXML_manager(hog).get_orthoxml_str()

    return ph, hog, iham_str

def compare_hog_structure(h1,h2):
    def append_child(current, list):
        list.append(current)
        return list

    v1 = h1.visit([], function_prefix=append_child)
    v2 = h2.visit([], function_prefix=append_child)

    v1 = sorted(v1,  key=lambda x: x.hog_id)
    v2 = sorted(v2, key=lambda x: x.hog_id)

    if len(v1) != len(v2):
        print(len(v1), len(v2))
        return False

    for i in range(len(v1)):
        hv1 = v1[i]
        hv2 = v2[i]

        if hv1.hog_id != hv2.hog_id:
            print(hv1.hog_id, hv2.hog_id)
            return False

        if len(hv1.children) != len(hv2.children):
            print(len(hv1.children),len(hv2.children))
            return False

        if len(hv1.duplications) != len(hv2.duplications):
            print(len(hv1.duplications),len(hv2.duplications))
            return False

        # CHECK FOR COMPOSITION OF GROUP

    return True

class IHAMTest(unittest.TestCase):

    def setUp(self):

        nwk_path = os.path.join(os.path.dirname(__file__), './data/simpleEx.nwk')
        tree_str = utils.get_newick_string(nwk_path, type="nwk")
        orthoxml_path = os.path.join(os.path.dirname(__file__), './data/simpleEx.orthoxml')

        self.ham_analysis = ham.Ham(tree_file=tree_str, hog_file=orthoxml_path, type_hog_file='orthoxml',
                                    use_internal_name=True)

    def test_ok(self):

        hog_3 = self.ham_analysis.get_hog_by_id(1)
        self.ham_analysis.create_iHam(hog=hog_3)
class OrthoXML_managerTest(unittest.TestCase):

    def test_orthoxml_str_simple_case(self):

        hog, hog2 = generate_hog_pair('./data/simpleEx.nwk', './data/simpleEx.orthoxml')

        self.assertTrue(compare_hog_structure(hog, hog2))

    def test_orthoxml_str_para_only_top(self):
        hog, hog2 = generate_hog_pair('./data/paralogs_only_toplevel_og.nwk', './data/paralogs_only_toplevel_og.orthoxml')

        self.assertTrue(compare_hog_structure(hog, hog2))

    def test_orthoxml_str_para_only_inside(self):
        hog, hog2 = generate_hog_pair('./data/paralogs_only_inside_og.nwk',
                                      './data/paralogs_only_inside_og.orthoxml')

        self.assertTrue(compare_hog_structure(hog, hog2))









