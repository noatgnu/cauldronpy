from unittest import TestCase
from cauldron.fuzz_clustering import fuzz_cluster


class Test(TestCase):
    def test_fuzz_cluster(self):
        fuzz_cluster(r"C:\Users\Toan Phung\Downloads\abundance_single-site_None (1).tsv", r"C:\Users\Toan Phung\Downloads\annotation.txt", "test", [3,4])
