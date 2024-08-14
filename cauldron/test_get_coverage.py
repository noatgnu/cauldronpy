from unittest import TestCase
from get_coverage import process_coverage

class Test(TestCase):
    def test_process_coverage(self):
        input_file = r"C:\Users\Toan Phung\Downloads\LRRK2 in Mouse Brain Lung and MEFs-DN.txt"
        value_columns = [
            "20230507_YL_MouseLung_LysoIP_LRRK2ko_13.raw",
            "20230507_YL_MouseLung_LysoIP_LRRK2ko_14.raw",
            "20230507_YL_MouseLung_LysoIP_LRRK2ko_15.raw",
            "20230507_YL_MouseLung_LysoIP_LRRK2ko_18.raw",
            "20230507_YL_MouseLung_LysoIP_LRRK2ko_20.raw",
            "20230507_YL_MouseLung_LysoIP_LRRK2ko_21.raw",
            "20230507_YL_MouseLung_LysoIP_WT_7.raw",
            "20230507_YL_MouseLung_LysoIP_WT_8.raw",
            "20230507_YL_MouseLung_LysoIP_WT_9.raw",
            "20230507_YL_MouseLung_LysoIP_WT_10.raw",
            "20230507_YL_MouseLung_LysoIP_WT_11.raw",
            "20230507_YL_MouseLung_LysoIP_WT_12.raw",
        ]

        index_column = "Precursor.Id"

        seq_column = "Stripped.Sequence"

        uniprot_acc_column = "Protein.Group"
        process_coverage(input_file, value_columns, index_column, seq_column, uniprot_acc_column, "", "test")