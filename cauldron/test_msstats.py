from unittest import TestCase
from cauldron.msstats import MSstats

class Test(TestCase):
    def test_diff_analysis(self):
        ms = MSstats(r'D:\WebstormProjects\cauldron\app\bin\win\R-Portable\App\R-Portable')