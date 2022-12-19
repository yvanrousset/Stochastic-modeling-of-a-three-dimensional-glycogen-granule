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
