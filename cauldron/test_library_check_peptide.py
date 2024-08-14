from unittest import TestCase
from cauldron.library_check_peptide import check_data_for_peptide_in_library, load_fasta_library


class Test(TestCase):
    def test_main(self):
        seq_df = load_fasta_library(r"C:\Users\Toan Phung\Downloads\20230102_UniprotSwissProt_Human_Cano+Iso.fasta", 2, 5)
        check_data_for_peptide_in_library(r"C:\Users\Toan Phung\Downloads\report.pr_matrix (4).tsv", "Stripped.Sequence", seq_df, r"test")
