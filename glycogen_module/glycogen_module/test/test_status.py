from unittest import TestCase
from glycogen_module.core import Status


class TestStatus(TestCase):
    def test_naming(self):
        gs = Status.GS_SUBSTRATE
        self.assertEqual(gs.value, 0)

        gbe = Status.GBE_SUBSTRATE
        self.assertEqual(gbe.value, 1)

        gp = Status.GP_SUBSTRATE
        self.assertEqual(gp.value, 2)

        gde = Status.GDE_SUBSTRATE
        self.assertEqual(gde.value, 3)

        expected = "0: GS_SUBSTRATE (substrate for elongation)\n1: GBE_SUBSTRATE (substrate for branching)\n2: GP_SUBSTRATE (substrate for reduction)\n3: GDE_SUBSTRATE (substrate for debranching)"
        self.assertEqual(Status.help(), expected)
