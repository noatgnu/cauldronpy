from unittest import TestCase
from cauldron.correlation_matrix import r_correlation_matrix, setup_r_home
import pandas as pd


class Test(TestCase):
    def test_r_correlation_matrix(self):
        setup_r_home(r'D:\WebstormProjects\cauldron\app\bin\win\R-Portable\App\R-Portable')
        df = pd.read_csv("test/annotation.txt", sep="\t")
        samples = df["Sample"].tolist()
        r_correlation_matrix("test/imputed.data.txt", "test", "Precursor.Id", samples)