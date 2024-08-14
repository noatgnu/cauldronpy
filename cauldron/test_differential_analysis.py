from unittest import TestCase
from cauldron.differential_analysis import diff_analysis, set_up_R_HOME

class Test(TestCase):
    def test_diff_analysis(self):
        set_up_R_HOME(r'D:\WebstormProjects\cauldron\app\bin\win\R-Portable\App\R-Portable')
        # diff_analysis(
        #     r"test/imputed.data.txt",
        #     r"test",
        #     r"test/annotation.txt",
        #     r"test/comparison.bca.txt",
        #     "Protein.Ids,Precursor.Id",
        #     True,
        # )
        # diff_analysis(
        #     r"C:\Users\Toan Phung\Downloads\tab_full_prot NEW (1).txt",
        #     r"test3",
        #     r"C:\Users\Toan Phung\Downloads\annotation.tab.txt",
        #     r"C:\Users\Toan Phung\Downloads\comparison.tab.txt",
        #     "protein",
        #     True,
        #     normalize="center.median",
        #     impute="RF"
        # )
        diff_analysis(
            r"\\mrc-smb.lifesci.dundee.ac.uk\mrc-group-folder\ALESSI\Toan\JT_SILAC_Secretome_DIA\report-first-pass.pr_matrix_channels_silac.tsv",
            r"test3",
            r"\\mrc-smb.lifesci.dundee.ac.uk\mrc-group-folder\ALESSI\Toan\JT_SILAC_Secretome_DIA\annotation.ratio.txt",
            r"\\mrc-smb.lifesci.dundee.ac.uk\mrc-group-folder\ALESSI\Toan\JT_SILAC_Secretome_DIA\comparison_matrix.txt",
            "Protein.Group,Precursor.Id",
            True,
            normalize="center.median",
            impute="RF",
            aggregate_column="Precursor.Id",
        )