from unittest import TestCase
from glycogen_module.core import Chain, Status
import json

class TestChain(TestCase):
    def test_serializable(self):
        c = Chain(0, [Status.GBE_SUBSTRATE, Status.GP_SUBSTRATE], 0, None, [], [], [[0, 0, 11.87], [0, 0, 12.87], [0, 0, 13.87]])
        expected = '{"id": 0, "substrate_status": ["GBE_SUBSTRATE", "GP_SUBSTRATE"], "generation": 0, "mother_id": null, "daughters_ids": [], "daughters_positions": [], "glucose_positions": [[0, 0, 11.87], [0, 0, 12.87], [0, 0, 13.87]]}'
        s = c._to_serializable()
        print(s)
        json_dumps_s = json.dumps(s)
        print(json_dumps_s)

        self.assertEqual(expected, json_dumps_s)

    def test_import_from_json(self):
        c = Chain(0, [Status.GBE_SUBSTRATE, Status.GP_SUBSTRATE], 0, None, [], [], [[0, 0, 11.87], [0, 0, 12.87], [0, 0, 13.87]])
        s = json.dumps(c._to_serializable())
        print(s)
        dct = json.loads(s)
        new_chain = Chain.from_dct(dct)
        print(new_chain)
        self.assertEqual(c, new_chain)
