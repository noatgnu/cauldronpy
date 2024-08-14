from unittest import TestCase

from cauldron.cv import main


class Test(TestCase):
    def test_diann_main(self):
        main(
            r"D:\PycharmProjects\dQ\data\dq\2be4b930-6cc8-4b37-8eae-1de59751be9b\data\Reports.log.txt",
            r"D:\PycharmProjects\dQ\data\dq\2be4b930-6cc8-4b37-8eae-1de59751be9b\data\Reports.pr_matrix.tsv",
            r"D:\PycharmProjects\dQ\data\dq\2be4b930-6cc8-4b37-8eae-1de59751be9b\data\Reports.pg_matrix.tsv",
            "Intensity",
            r"D:\PycharmProjects\dQ\data\dq\2be4b930-6cc8-4b37-8eae-1de59751be9b\data\annotation.txt", "", "")

