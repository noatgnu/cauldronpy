from unittest import TestCase

from cauldron.cv import extract_filename
from cauldron.normalization import normalize

class Test(TestCase):
    def test_normalize(self):
        log_file = r"D:\PycharmProjects\dQ\data\dq\2be4b930-6cc8-4b37-8eae-1de59751be9b\data\Reports.log.txt"
        extracted_filename = extract_filename(log_file)
        normalize(r"D:\PycharmProjects\dQ\data\dq\2be4b930-6cc8-4b37-8eae-1de59751be9b\data\Reports.pr_matrix.tsv",
                  "test",
                  ",".join(extracted_filename),

                  )
