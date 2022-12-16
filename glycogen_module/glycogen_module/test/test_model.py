from unittest import TestCase
from glycogen_module.core import Model


class TestModel(TestCase):
    def test_hashable(self):
        m = Model.FLEXIBLE_LOCATION_MODEL
        self.assertEqual(hash(m), hash(m.name))
