from unittest import TestCase

from estimation_plot import est_plot

class TestEstPlot(TestCase):
    def test_est_plot(self):
        est_plot(
            r"test\imputed.data.txt",
            "Precursor.Id",
            ["ALIDPSSGLPNRLPPGAVPPGAR3", "IVSGIITPIHEQWEK3", "SELLPAGWNNNK2"],
            r"test\annotation.txt",
            "test",
            log2=True,
           condition_order=['Pepide-CBQCA_LT-WCL','BCA_LT-IP', 'BCA_LT-MockIP', 'Pepide-CBQCA_LT-IP', 'Pepide-CBQCA_LT-MockIP', 'BCA_LT-WCL']
        )
