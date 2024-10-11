from unittest import TestCase
from cauldron.msstats import MSstats

class Test(TestCase):
    def test_diff_analysis(self):
        ms = MSstats(r'D:\WebstormProjects\cauldron\app\bin\win\R-Portable\App\R-Portable')
        ms.run_msstats(annotation_file=r"Z:\ALESSI\Toan\RN-AH_TDP43_Mice_Glia_02\annotation.txt", report_file=r"Z:\ALESSI\Toan\RN-AH_TDP43_Mice_Glia_02\report.tsv")