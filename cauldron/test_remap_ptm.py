from unittest import TestCase
from remap_ptm import process

class Test(TestCase):
    def test_process(self):
        process(
            "",
            "test/different (39).txt",
            "Peptide",
            "ProteinID",
            "Position.in.peptide",
            "test"
        )
