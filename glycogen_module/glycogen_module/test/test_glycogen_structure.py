import json
from pathlib import Path
from unittest import TestCase
from glycogen_module.core import GlycogenStructure, Model
import pytest


class TestGlycogenStructure(TestCase):
    SCRIPTDIR = Path(__file__).parent.resolve()
    PARAMS_PATH = Path(f"{SCRIPTDIR}/testdata/parameters_1.json")
    """ testdata/parameters_1.json:
    {"GBE": 1.0,"GS": 0.2,"GDE": 0.0,"GP": 0.0,"l_gbe_spacing": 1,"l_gbe_leftover": 3,"l_gbe_transferred": 3,"model_for_gbe": "flexible_location_model","radius_of_glucose_sphere": 5.4,"radius_of_gn_core": 11.87}
    """
    EXPECTED_PARAMS_1 = {"gbe": 1.0, "gs": 0.2, "gde": 0.0, "gp": 0.0, "l_gs_min": 2, "l_gbe_spacing": 1, "l_gbe_leftover": 3,
                         "l_gbe_transferred": 3, "model_for_gbe": "flexible_location_model", "radius_of_glucose_sphere": 5.4, "radius_of_gn_core": 11.87}

    EXPECTED_INITIAL_GLUCOSE_LOCATIONS_C0 = [[0, 0, 11.87], [
        0, 0, 12.87], [0, 0, 13.87], [0, 0, 14.87], [0, 0, 15.87]]
    EXPECTED_INITIAL_GLUCOSE_LOCATIONS_C1 = [
        [0, 0, -11.87], [0, 0, -12.87], [0, 0, -13.87], [0, 0, -14.87], [0, 0, -15.87]]

    OUTFILE_TMP_JSON = Path(f"{SCRIPTDIR}/testdata/tmp_out.json")
    JSON_SNAPSHOT = Path(f"{SCRIPTDIR}/testdata/params_and_chains.json")

    def test_init(self):
        gs = GlycogenStructure(
            gbe=1.0,
            gs=0.2,
            gde=0.0,
            gp=0.0,
            l_gs_min=2,
            l_gbe_spacing=1,
            l_gbe_leftover=3,
            l_gbe_transferred=3,
            model_for_gbe=Model.FLEXIBLE_LOCATION_MODEL,
            radius_of_glucose_sphere=5.4,
            radius_of_gn_core=11.87
        )
        self.assertEqual(len(gs.chains), 2)
        for c in gs.chains:
            print(c.glucose_positions)

        self.assertEqual(gs.chains[0].glucose_positions,
                         TestGlycogenStructure.EXPECTED_INITIAL_GLUCOSE_LOCATIONS_C0)
        self.assertEqual(gs.chains[1].glucose_positions,
                         TestGlycogenStructure.EXPECTED_INITIAL_GLUCOSE_LOCATIONS_C1)

    def test_from_dct(self):
        with open(TestGlycogenStructure.PARAMS_PATH, 'r') as infile:
            d = json.load(infile)
        gs = GlycogenStructure.from_dct(d)
        print(gs)

    def test_from_json_file(self):
        gs = GlycogenStructure.from_json_file(
            TestGlycogenStructure.PARAMS_PATH)
        print(gs)

        self.assertEqual(
            gs.gbe, TestGlycogenStructure.EXPECTED_PARAMS_1["gbe"])
        self.assertEqual(
            gs.gde,  TestGlycogenStructure.EXPECTED_PARAMS_1["gde"])
        self.assertEqual(gs.gp,  TestGlycogenStructure.EXPECTED_PARAMS_1["gp"])
        self.assertEqual(gs.gs,  TestGlycogenStructure.EXPECTED_PARAMS_1["gs"])

        self.assertEqual(
            gs.l_gs_min, TestGlycogenStructure.EXPECTED_PARAMS_1["l_gs_min"])
        self.assertEqual(
            gs.l_gbe_spacing, TestGlycogenStructure.EXPECTED_PARAMS_1["l_gbe_spacing"])
        self.assertEqual(
            gs.l_gbe_leftover, TestGlycogenStructure.EXPECTED_PARAMS_1["l_gbe_leftover"])
        self.assertEqual(
            gs.l_gbe_transferred, TestGlycogenStructure.EXPECTED_PARAMS_1["l_gbe_transferred"])
        self.assertEqual(
            gs.model_for_gbe, TestGlycogenStructure.EXPECTED_PARAMS_1["model_for_gbe"])
        self.assertEqual(gs.radius_of_glucose_sphere,
                         TestGlycogenStructure.EXPECTED_PARAMS_1["radius_of_glucose_sphere"])
        self.assertEqual(
            gs.radius_of_gn_core, TestGlycogenStructure.EXPECTED_PARAMS_1["radius_of_gn_core"])

    def test_find_chains_for_gs(self):
        expected = [0, 1]
        gs = GlycogenStructure.from_json_file(
            TestGlycogenStructure.PARAMS_PATH)
        print(gs)
        found = gs.find_chains_for_gs()
        self.assertEqual(found, expected)

    def test_find_chains_for_gbe(self):
        expected = []
        gs = GlycogenStructure.from_json_file(
            TestGlycogenStructure.PARAMS_PATH)
        print(gs)
        found = gs.find_chains_for_gbe()
        self.assertEqual(found, expected)

    def test_find_chains_for_gde(self):
        expected = []
        gs = GlycogenStructure.from_json_file(
            TestGlycogenStructure.PARAMS_PATH)
        print(gs)
        found = gs.find_chains_for_gde()
        self.assertEqual(found, expected)

    def test_find_chains_for_gp(self):
        expected = [0, 1]
        gs = GlycogenStructure.from_json_file(
            TestGlycogenStructure.PARAMS_PATH)
        print(gs)
        found = gs.find_chains_for_gp()
        self.assertEqual(found, expected)

    def test_act_gs(self):
        """Example from tutorial: act_gs() will try to add a glucose unit on top of a chain. """

        g = GlycogenStructure.from_json_file(
            TestGlycogenStructure.PARAMS_PATH)
        print(g)
        fixed_glucose_start = sum([len(c.glucose_positions) for c in g.chains])
        print(fixed_glucose_start)
        for _ in range(30):
            g.act_gs()
            print(g.chains)
        # for c in g.chains:
        #    print(f"{c.id} {c.glucose_positions}")
        fixed_glucose = sum([len(c.glucose_positions) for c in g.chains])
        self.assertEqual(fixed_glucose-fixed_glucose_start, 30)

    def test_get_num_glucose_fixed(self):
        """Example from tutorial"""

        g = GlycogenStructure.from_json_file(
            TestGlycogenStructure.PARAMS_PATH)
        fixed_glucose_start = g.get_num_glucose_fixed()
        print(fixed_glucose_start)
        for _ in range(30):
            g.act_gs()
        self.assertEqual(g.get_num_glucose_fixed()-fixed_glucose_start, 30)

    def test_gbe(self):
        # TODO: this is a very weak and lazy test...
        # should be stricter with forced parameters
        g = GlycogenStructure.from_json_file(
            TestGlycogenStructure.PARAMS_PATH)
        for _ in range(30):
            g.act_gs()
        init_num_chains = g.get_num_chains()

        print(f"Initial number of chains: {init_num_chains}")
        while len(g.find_chains_for_gbe()) > 0:
            g.act_gbe_flexible_model()

        result_num_chains = g.get_num_chains()
        print(f"Final number of chains: {result_num_chains}")

        self.assertGreater(result_num_chains, init_num_chains)

    def test_act_gp(self):
        # TODO: this is a very weak and lazy test...
        # should be stricter with forced parameters
        g = GlycogenStructure.from_json_file(
            TestGlycogenStructure.PARAMS_PATH)
        for _ in range(30):
            g.act_gs()

        while len(g.find_chains_for_gbe()) > 0:
            g.act_gbe_flexible_model()

        before = g.get_num_glucose_fixed()
        while len(g.find_chains_for_gp()) > 0:
            g.act_gp()

        result = g.get_num_glucose_fixed()
        self.assertLess(result, before)

    def test_act_gde(self):
        # TODO: this is a very weak and lazy test...
        # should be stricter with forced parameters
        g = GlycogenStructure.from_json_file(
            TestGlycogenStructure.PARAMS_PATH)
        for _ in range(30):
            g.act_gs()

        while len(g.find_chains_for_gbe()) > 0:
            g.act_gbe_flexible_model()

        while len(g.find_chains_for_gp()) > 0:
            g.act_gp()

        num_glucose_before = g.get_num_glucose_fixed()
        num_chains_before = g.get_num_chains()
        while len(g.find_chains_for_gde()) > 0:
            g.act_gde()

        num_chains_after = g.get_num_chains()
        num_glucose_after = g.get_num_glucose_fixed()

        difference_num_chains = num_chains_before-num_chains_after
        expected_num_glucose = num_glucose_before-difference_num_chains

        # Since a debranching release one glucose,
        # the final number of glucose should be the number of glucose after the reduction,
        # proceed by GP, minus the number of debranched chains.
        self.assertEqual(num_glucose_after, expected_num_glucose)
        self.assertLess(num_chains_after, num_chains_before)

    def test_export_parameters_back_to_dct(self):
        with open(TestGlycogenStructure.PARAMS_PATH, 'r') as infile:
            d = json.load(infile)
        s = GlycogenStructure.from_dct(d)
        print(s)
        result = s.get_parameters_as_dct()
        new_obj = GlycogenStructure.from_dct(result)
        print(new_obj)
        self.assertEqual(new_obj.__dict__, s.__dict__)

    def test_export_to_file(self):
        with open(TestGlycogenStructure.PARAMS_PATH, 'r') as infile:
            d = json.load(infile)
        s = GlycogenStructure.from_dct(d)
        s.act_gs()
        s.export_to_file(TestGlycogenStructure.OUTFILE_TMP_JSON, format='json')

        #with open(TestGlycogenStructure.OUTFILE_TMP_JSON, 'r') as tmp:
        #    print(tmp.read())

        r = GlycogenStructure.from_json_file(
            TestGlycogenStructure.OUTFILE_TMP_JSON, no_init=True)

        self.assertEqual(s.chains, r.chains)


    def test_import_from_json_no_init(self):
        with open(TestGlycogenStructure.PARAMS_PATH, 'r') as infile:
            d = json.load(infile)
        s1 = GlycogenStructure.from_dct(d)

        s2 = GlycogenStructure.from_json_file(
            TestGlycogenStructure.JSON_SNAPSHOT, no_init=True)
        s3 = GlycogenStructure.from_json_file(
            TestGlycogenStructure.JSON_SNAPSHOT, no_init=False)
        print(s1)
        print(s2)
        print(s1.chains)
        print(s2.chains)
        # no_init=True keeps the chains from the snapshot
        self.assertNotEqual(s1.chains, s2.chains)
        # no_init=False leads to overwriting the chains from the snapshot
        self.assertEqual(s1.chains, s3.chains)

        # self.assertTrue(False)
