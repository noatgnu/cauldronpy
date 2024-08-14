from unittest import TestCase

from cauldron.cv import extract_filename
from cauldron.pca import pca_

class Test(TestCase):
    def test_pca_(self):
        log_file = r"D:\PycharmProjects\dQ\data\dq\2be4b930-6cc8-4b37-8eae-1de59751be9b\data\Reports.log.txt"
        extracted_filename = extract_filename(log_file)
        pca_(r"test/imputed.data.txt", "test", extracted_filename, 2, True)
