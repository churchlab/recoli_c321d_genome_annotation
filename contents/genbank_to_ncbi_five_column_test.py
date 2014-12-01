"""
Tests for genbank_to_ncbi_five_column.py.
"""

import unittest

from genbank_to_ncbi_five_column import _handle_go_term


class TestGenbankToNCBIFiveColumn(unittest.TestCase):

    def test_handle_go_term(self):
        """Tests getting the line to represent a GO term.
        """
        KEY = 'GO_process'
        VALUE = 'GO:0009088 - threonine biosynthetic process'
        EXPECTED_LINE_VALUE = (
                '\t\t\tGO_process\tthreonine biosynthetic process|0009088')
        self.assertEqual(EXPECTED_LINE_VALUE, _handle_go_term(KEY, VALUE))

    def test_handle_go_term__multiple(self):
        """Handle multiple GO terms in a single line.
        """
        # TODO: Implement.
        pass


if __name__ == '__main__':
    unittest.main()
