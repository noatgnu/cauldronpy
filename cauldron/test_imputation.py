from unittest import TestCase
from cauldron import imputation
from cauldron.cv import extract_filename


# example data for testing imputation



class Test(TestCase):
    def test_impute(self):
        """
        Test the impute function from cauldron/imputation.py
        """
        log_file = r"D:\PycharmProjects\dQ\data\dq\2be4b930-6cc8-4b37-8eae-1de59751be9b\data\Reports.log.txt"
        extracted_filename = extract_filename(log_file)
        imputation.impute(r"D:\PycharmProjects\dQ\data\dq\2be4b930-6cc8-4b37-8eae-1de59751be9b\data\Reports.pr_matrix.tsv",
                          r"test",
                            ",".join(extracted_filename))


