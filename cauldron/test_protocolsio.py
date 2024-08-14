from unittest import TestCase
from cauldron.protocolsio import extract_json

class Test(TestCase):
    def test_extract_json(self):
        extract_json(50519)
