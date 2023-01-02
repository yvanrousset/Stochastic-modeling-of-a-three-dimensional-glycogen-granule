import json
from pathlib import Path
from unittest import TestCase
from glycogen_module.core import GlycogenStructure, Model
import pytest


class TestPrivMethods(TestCase):
    SCRIPTDIR = Path(__file__).parent.resolve()
    PARAMS_PATH = Path(f"{SCRIPTDIR}/testdata/parameters_1.json")

    def test_act_gbe_once(self):
        g = GlycogenStructure(gs=0, gbe=0, gde=0, gp=0, l_gs_min=2, l_gbe_spacing=2, l_gbe_leftover=2, l_gbe_transferred=2,
                              model_for_gbe='flexible_location', radius_of_glucose_sphere=5.4, radius_of_gn_core=11.87)

        print(g.l_gbe_spacing, g.l_gbe_transferred, g.l_gbe_leftover)
        for _ in range(2):
            g.act_gs(chain_id=0)

        # amount of glucose should stay the same
        num_glucose_before = g.get_num_glucose_fixed()
        num_chains_before = g.get_num_chains()
        for c in g.chains:
            print(
                f"After 2x GS:{c.id} {c.daughters_ids}{c.daughters_positions} {c.get_num_glucose_positions()} \n{c.glucose_positions}")

        g.act_gbe_flexible_model(chain_id=0, force_branch_index=2,
                                 force_cleave_index=4, force_theta=-3.0, force_phi=1.5)

        for c in g.chains:
            print(
                f"After 1x GBE: {c.id} {c.daughters_ids}{c.daughters_positions} {c.get_num_glucose_positions()}\n{c.glucose_positions}")
        num_chains_after = g.get_num_chains()
        num_glucose_after = g.get_num_glucose_fixed()

        print(g.chains)
        self.assertEqual(num_glucose_before, num_glucose_after)

        self.assertEqual(num_chains_after, num_chains_before+1)

        self.assertEqual(g.get_chain_by_id(0).get_num_glucose_positions(
        ), 5)
        self.assertEqual(g.get_chain_by_id(g.get_chain_by_id(0).daughters_ids[0]).get_num_glucose_positions(
        ), 2)


class TestMiscMethods(TestCase):
    SCRIPTDIR = Path(__file__).parent.resolve()
    PARAMS_PATH = Path(f"{SCRIPTDIR}/testdata/parameters_1.json")
    PREMADE_G = Path(f"{SCRIPTDIR}/testdata/test_a_b_zero.json")
    PREMADE_G_7_CHAINS = Path(f"{SCRIPTDIR}/testdata/test_7chains.json")
    PREMADE_TAB1 = Path(
        f"{SCRIPTDIR}/testdata/g_N2000_gamma_0.2_1.0_exported_tab1_seed123.json")

    def test_ab_ratio_zero(self):
        g = GlycogenStructure.from_json_file(
            TestMiscMethods.PREMADE_G, no_init=True)
        print(g)
        with pytest.raises(Exception) as e:
            print(g.get_a_b_ratio())
            self.assertRaises(expected_exception=Exception)
        assert str(
            e.value) == "Cannot return A:B ratio, no branches with daughters present."

    def test_ab_ratio_7chains(self):
        g = GlycogenStructure.from_json_file(
            TestMiscMethods.PREMADE_G_7_CHAINS, no_init=True)
        print(g)
        self.assertEqual(g.get_a_b_ratio(), 0.4)

    def test_ab_ratio_tab1(self):
        g = GlycogenStructure.from_json_file(
            TestMiscMethods.PREMADE_TAB1, no_init=True)
        print(g)
        expected = 0.828125
        self.assertEqual(g.get_a_b_ratio(), expected)

    def test_avg_chain_length_init(self):
        g = GlycogenStructure.from_json_file(
            TestMiscMethods.PARAMS_PATH)
        print(g)
        result = g.get_avg_chain_length()
        expected = 5.0
        self.assertEqual(result, expected)

    def test_avg_chain_length(self):
        g = GlycogenStructure.from_json_file(
            TestMiscMethods.PREMADE_G_7_CHAINS, no_init=True)
        print(g)
        result = g.get_avg_chain_length()
        expected = 8.571428571428571
        self.assertEqual(result, expected)

    def test_avg_chain_length_tab1(self):
        g = GlycogenStructure.from_json_file(
            TestMiscMethods.PREMADE_TAB1, no_init=True)
        print(g)
        result = g.get_avg_chain_length()
        expected = 8.547008547008547
        self.assertEqual(result, expected)

    def test_branching_degree(self):
        g = GlycogenStructure.from_json_file(
            TestMiscMethods.PREMADE_G_7_CHAINS, no_init=True)
        print(g)
        result = g.get_branching_degree()
        expected = 0.11320754716981132
        self.assertEqual(result, expected)

    def test_branching_degree_tab1(self):
        g = GlycogenStructure.from_json_file(
            TestMiscMethods.PREMADE_TAB1, no_init=True)
        print(g)
        result = g.get_branching_degree()
        expected = 0.1319365798414496
        self.assertEqual(result, expected)

    def test_last_gen_index_tab1(self):
        g = GlycogenStructure.from_json_file(
            TestMiscMethods.PREMADE_TAB1, no_init=True)
        print(g)
        result = g.get_last_gen_index()
        expected = 9
        self.assertEqual(result, expected)

    def test_occupancy_tab1(self):
        g = GlycogenStructure.from_json_file(
            TestMiscMethods.PREMADE_TAB1, no_init=True)
        print(g)
        result = g.get_occupancy()
        expected = 0.2958106823848145
        self.assertAlmostEqual(result, expected)

    def test_radius_tab1(self):
        g = GlycogenStructure.from_json_file(
            TestMiscMethods.PREMADE_TAB1, no_init=True)
        print(g)
        result = g.get_radius(unit='nm')
        expected = 8.016721661194044
        self.assertAlmostEqual(result, expected)

