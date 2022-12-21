import json
from unittest import TestCase
from glycogen_module.core import Status


class TestStatus(TestCase):
    def test_naming(self):
        gs = Status.GS_SUBSTRATE
        self.assertEqual(gs.value, "GS_SUBSTRATE")

        gbe = Status.GBE_SUBSTRATE
        self.assertEqual(gbe.value, "GBE_SUBSTRATE")

        gp = Status.GP_SUBSTRATE
        self.assertEqual(gp.value, "GP_SUBSTRATE")

        gde = Status.GDE_SUBSTRATE
        self.assertEqual(gde.value, "GDE_SUBSTRATE")

        expected = "GS_SUBSTRATE: substrate for elongation\nGBE_SUBSTRATE: substrate for branching\nGP_SUBSTRATE: substrate for reduction\nGDE_SUBSTRATE: substrate for debranching"
        self.assertEqual(Status.help(), expected)

    def test_serialize(self):
        s = Status.GBE_SUBSTRATE

        print(s)
        json_dumps_s = json.dumps(s)
        print(json_dumps_s)
        json_loads_s = json.loads('"GBE_SUBSTRATE"')
        string_s = Status("GBE_SUBSTRATE")

        self.assertEqual(json.loads(json_dumps_s) , string_s)
        self.assertEqual(Status(string_s), Status.GBE_SUBSTRATE)
