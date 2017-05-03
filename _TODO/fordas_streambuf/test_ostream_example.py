# vim: ai ts=4 sts=4 et sw=4 ft=python

import io 
import unittest

class TestOStreamExample(unittest.TestCase):
    def test_ostream(self):
        import ostream_example

        iob = io.BytesIO()
        ostream_example.testprint_noflush(iob)
        ostream_example.testprint(iob)
        assert iob.getvalue() == b"testprint_noflushtestprint\n"

    def test_istream(self):
        import ostream_example
        assert ostream_example.testparse(io.BytesIO(b"4")) == 4

if __name__ == '__main__':
    unittest.main()
