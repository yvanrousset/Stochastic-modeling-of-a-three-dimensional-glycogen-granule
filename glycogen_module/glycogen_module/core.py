from dataclasses import dataclass
from datetime import datetime
from enum import Enum
import json

import logging
import math
from pathlib import Path
import random
import sys

from glycogen_module.utils import angle3Dchain
import numpy as np
logger = logging.getLogger(__name__)


class Model(str, Enum):
    """Enum representing the branching model. """
    FLEXIBLE_LOCATION_MODEL = "flexible_location_model"
    STRICT = "strict"

    def __hash__(self) -> int:
        return hash(self.name)

    def __str__(self) -> str:
        return self.value


class Status(str, Enum):
    """Enum representing the status of a chain, showing which enzymatic reactions can take place."""
    GS_SUBSTRATE = "GS_SUBSTRATE"
    GBE_SUBSTRATE = "GBE_SUBSTRATE"
    GP_SUBSTRATE = "GP_SUBSTRATE"
    GDE_SUBSTRATE = "GDE_SUBSTRATE"

    @classmethod
    def help(cls) -> str:
        """Show help message."""
        msg = []
        for i, m in zip([Status.GS_SUBSTRATE, Status.GBE_SUBSTRATE, Status.GP_SUBSTRATE, Status.GDE_SUBSTRATE], ["elongation", "branching", "reduction", "debranching"]):
            msg.append(f"{cls(i).name}: substrate for {m}")
        msg = "\n".join(msg)
        print(msg)
        return (msg)


@dataclass
class Chain:
    """Represents a chain of glucose units."""

    id: int
    substrate_status: list[Status]
    generation: int
    mother_id: int
    daughters_ids: list[int]
    daughters_positions: list[int]
    glucose_positions: list[int]

    def get_num_glucose_positions(self):
        return len(self.glucose_positions)

    def get_num_daughters(self):
        return len(self.daughters_ids)

    def _to_serializable(self):
        return self.__dict__

    @classmethod
    def from_dct(cls, dct: dict):
        return cls(
            id=dct['id'],
            substrate_status=[Status(v) for v in dct['substrate_status']],
            generation=dct['generation'],
            mother_id=dct['mother_id'],
            daughters_ids=dct['daughters_ids'],
            daughters_positions=dct['daughters_positions'],
            glucose_positions=dct['glucose_positions']
        )


class GlycogenStructure:
    """Represents a glycogen structure."""

    # 0.65nm The radius of the circular cross-section
    RHO = 0.65

    # 0.24 nm Each glucose residue contributes to the chain length by l = 0.24 nm
    L = 0.24  # appears to be used while plotting

    # Approximate radius  ρGN ≈ 2.85 nm
    RHO_GN = 2.85

    # we propose that monomers are described as overlapping spheres with radius ρ = 0.65 nm, equal to that of the helix.
    MONOMER_RADIUS = 0.65

    # This model currently only handles max. 2 initial chains
    MAX_INITIAL_CHAINS = 2  # This model currently only handles max. 2 initial chains

    # Fixed to 1 in this model
    # unitless distance between two glucoses units organised in a single helix TODO: is this not * L for the coordinates?
    DISTANCE_BETWEEN_GLUCOSE_UNITS = 1  # TODO: rename this to [...]MULTIPLIER

    @classmethod
    def from_json_file(cls, file_path: Path | str, no_init=False):
        """ If no_init is set to False, the initial chains will not be created. Use no_init when you are importing a snapshot of a simulation. """
        try:
            with open(file=file_path, mode='r') as infile:
                result = json.load(infile)

        except Exception as e:
            print(
                f"{e} - Could not load json from file with path '{file_path}'", file=sys.stderr)
            raise
        return cls.from_dct(parameters_dct=result, no_init=no_init)

    @classmethod
    def from_dct(cls, parameters_dct: dict, no_init=False):
        """ If no_init is set to True, the initial chains will not be created. Use no_init when you are importing a snapshot of a simulation. """
        lowercase_converted_dct = {
            k.lower(): v for k, v in parameters_dct.items()}

        available_models = [member.value for member in Model]
        try:
            parameters_dct['model_for_gbe'] in available_models
        except ValueError as e:
            print(
                f"{e}: {parameters_dct['model_for_gbe']} is not a valid model. Available models are: {available_models}", file=sys.stderr)
            raise ValueError

        if no_init:
            chains = []
            tmp = lowercase_converted_dct['chains']
            for t in tmp:
                new_chain = Chain.from_dct(t)
                chains.append(new_chain)
        else:
            chains = None
        return cls(
            gs=lowercase_converted_dct['gs'],
            gbe=lowercase_converted_dct['gbe'],
            gde=lowercase_converted_dct['gde'],
            gp=lowercase_converted_dct['gp'],

            l_gs_min=lowercase_converted_dct['l_gs_min'],

            l_gbe_spacing=lowercase_converted_dct['l_gbe_spacing'],
            l_gbe_transferred=lowercase_converted_dct['l_gbe_transferred'],
            l_gbe_leftover=lowercase_converted_dct['l_gbe_leftover'],
            model_for_gbe=Model(
                lowercase_converted_dct['model_for_gbe']),
            # Paper: ??
            # distance_between_units=lowercase_converted_dct['distance_between_units'],
            # used in scan_excluded_volume
            radius_of_glucose_sphere=lowercase_converted_dct['radius_of_glucose_sphere'],
            radius_of_gn_core=lowercase_converted_dct['radius_of_gn_core'],
            chains=chains
        )

    def _create_initial_chains(self, num_chains: int):
        chains = []
        if num_chains > GlycogenStructure.MAX_INITIAL_CHAINS:
            raise Exception(
                "Not Implemented: Only two inital chains are allowed in this version of the model.")
        for i in range(num_chains):
            if i == 0:
                new_chain = Chain(id=i,
                                  substrate_status=[
                                      Status.GP_SUBSTRATE, Status.GS_SUBSTRATE],
                                  generation=0,
                                  mother_id=None,
                                  daughters_ids=[],
                                  daughters_positions=[],
                                  glucose_positions=[
                                      # TODO: look at act_gs for the distance parameter...
                                      [0, 0, self.radius_of_gn_core],
                                      [0, 0, self.radius_of_gn_core+1 * \
                                          GlycogenStructure.DISTANCE_BETWEEN_GLUCOSE_UNITS],
                                      [0, 0, self.radius_of_gn_core+2 * \
                                          GlycogenStructure.DISTANCE_BETWEEN_GLUCOSE_UNITS],
                                      [0, 0, self.radius_of_gn_core+3 * \
                                          GlycogenStructure.DISTANCE_BETWEEN_GLUCOSE_UNITS],
                                      [0, 0, self.radius_of_gn_core+4 * \
                                          GlycogenStructure.DISTANCE_BETWEEN_GLUCOSE_UNITS]
                                  ])
                chains.append(new_chain)
            elif i == 1:
                new_chain = Chain(id=i,
                                  substrate_status=[
                                      Status.GP_SUBSTRATE, Status.GS_SUBSTRATE],
                                  generation=0,
                                  mother_id=None,
                                  daughters_ids=[],
                                  daughters_positions=[],
                                  glucose_positions=[
                                      # TODO: look at act_gs for the distance parameter
                                      [0, 0, - self.radius_of_gn_core],
                                      [0, 0, - self.radius_of_gn_core-1 * \
                                          GlycogenStructure.DISTANCE_BETWEEN_GLUCOSE_UNITS],
                                      [0, 0, - self.radius_of_gn_core-2 * \
                                          GlycogenStructure.DISTANCE_BETWEEN_GLUCOSE_UNITS],
                                      [0, 0, - self.radius_of_gn_core-3 * \
                                          GlycogenStructure.DISTANCE_BETWEEN_GLUCOSE_UNITS],
                                      [0, 0, - self.radius_of_gn_core-4 * \
                                          GlycogenStructure.DISTANCE_BETWEEN_GLUCOSE_UNITS]
                                  ])
                chains.append(new_chain)
                return chains

    def __init__(self, gs: float,
                 gbe: float,
                 gde: float,
                 gp: float,
                 l_gs_min: int,
                 l_gbe_spacing: int,
                 l_gbe_transferred: int,
                 l_gbe_leftover: int,
                 model_for_gbe: Model,
                 radius_of_glucose_sphere: float,
                 radius_of_gn_core: float,
                 chains: list[Chain] = None
                 ) -> None:
        """
        - gs:float GS activity
        - gbe:float GBE activity
        - gde:float GDE activity
        - gp:float GP activity

        - l_gs_min:int Minimal chain length of required primer by GS (was size_spec_gys )

        #- l_gp: int # TODO?
        #- l_gde: int # TODO?

        #- l_gbe_min:int Minimal chain length for the substrate such that GBE is able to bind TODO
        - l_gbe_spacing:int New α-1,6 branching point must not be closer than l_gbe_spacing from either the first glucose of the chain, or an α-1,6 branching point.
        - l_gbe_transferred:int After binding, GBE cleaves at least L GBE transferred glucose units.
        - l_gbe_leftover:int A new branching point must not be closer than L_GBE_leftover from the non-reducing end of the substrate chain, which is the original position of cleavage.

        - model_for_gbe:Model Branching model,
        #- distance_between_units:float TODO: sth like how many times radius of glucose sphere
        #- monomer_radius:float TODO: Paper: rho
        - radius_of_glucose_sphere: float
        - radius_gn_core The radius of the glycogenin core in multiples of rise per turn

        """
        logger.debug(f"Init GlycogenStructure")

        self.gs = gs
        self.gbe = gbe
        self.gde = gde
        self.gp = gp
        self.l_gs_min = l_gs_min
        self.l_gbe_spacing = l_gbe_spacing
        self.l_gbe_transferred = l_gbe_transferred
        self.l_gbe_leftover = l_gbe_leftover
        self. model_for_gbe = model_for_gbe
        self.radius_of_glucose_sphere = radius_of_glucose_sphere
        self.radius_of_gn_core = radius_of_gn_core

        # by default create inital chains
        if chains is None:
            self.chains: list[Chain] = self._create_initial_chains(
                GlycogenStructure.MAX_INITIAL_CHAINS)
        # if chains are passed as parameter, do not create inital chains, but use provided chains
        else:
            self.chains = chains

    def find_chains_for_gs(self) -> list[int]:
        return [c.id for c in self.chains if Status.GS_SUBSTRATE in c.substrate_status]

    def find_chains_for_gbe(self) -> list[int]:
        return [c.id for c in self.chains if Status.GBE_SUBSTRATE in c.substrate_status]

    def find_chains_for_gde(self) -> list[int]:
        return [c.id for c in self.chains if Status.GDE_SUBSTRATE in c.substrate_status]

    def find_chains_for_gp(self) -> list[int]:
        return [c.id for c in self.chains if Status.GP_SUBSTRATE in c.substrate_status]

    def find_chains_for_gde(self) -> list[int]:
        return [c.id for c in self.chains if Status.GDE_SUBSTRATE in c.substrate_status]

    def get_chain_by_id(self, id: int) -> Chain:
        """Returns Chain object if chain with given id exists in the GlycogenStructure instance, otherwise returns None """
        try:
            return [c for c in self.chains if c.id == id][0]
        except IndexError as e:
            logger.debug(f"No chain with id {id} present")
            return None
        # TODO: switch to dict for efficiency or return by index if things always stay in order

    def get_num_glucose_fixed(self) -> int:
        return sum([len(c.glucose_positions) for c in self.chains])

    def get_num_chains(self) -> int:
        return len(self.chains)

    def get_last_chain_index(self) -> int:
        # TODO: not sure if necessary to find max, maybe just return last index...
        return max([c.id for c in self.chains])

    def create_a_chain(self, new_chain_index: int, status: Status | None = None, generation: int | None = None):  # TODO
        """Appends a new chain onto the GlycogenStructure instance.
        Raises Exception when chain with given id is already present.
         """
        if self.get_chain_by_id(new_chain_index):
            raise Exception(
                f"Chain with id {new_chain_index} already present. ")
        new_chain = Chain(id=new_chain_index, substrate_status=status, generation=generation,
                          mother_id=None, daughters_ids=[], daughters_positions=[], glucose_positions=[])

        self.chains.append(new_chain)

    def update_chain_status(self, chain_id: int) -> None:
        """needs:
        - l_gs_min
        - l_gbe_spacing
        - l_gbe_leftover
        - l_gbe_transferred
         """
        chain = self.get_chain_by_id(chain_id)
        N = len(chain.glucose_positions)
        # why does it check for the positions not ids like below?
        num_pos_daughters = len(chain.daughters_positions)
        if num_pos_daughters == 0:
            ind = 0
        else:
            # last index +1... last branching index(?)
            ind = chain.daughters_positions[-1]+1

        logger.debug(
            f"Update chain status : N:{N}, ind: {ind}, self.l_gs_min: {self.l_gs_min}, N-ind: {N-ind} ")

        if N-ind < self.l_gs_min:  # size spec glycogen synthase
            # was '1' Only GS substrate
            chain.substrate_status = [Status.GS_SUBSTRATE]

        # as long as l_gs_min and no daughters
        elif N-ind == self.l_gs_min and len(chain.daughters_ids) == 0:
            # was '4' The chain is really short AND does not have any branching points (rather: any daughter chains?).
            # It is therefore subject both to elongation (enzyme GS) and debranching (enzyme GDE) Substrate for : GS, GDE
            chain.substrate_status = [
                Status.GS_SUBSTRATE, Status.GDE_SUBSTRATE]  # does that make any sense? How will it be debranched without branching points?

        elif N-ind <= self.l_gbe_spacing + self.l_gbe_transferred + self.l_gbe_leftover:
            # From Tutorial (https://github.com/yvanrousset/Stochastic-modeling-of-a-three-dimensional-glycogen-granule/blob/main/tutorial_synthesis_and_degradation.ipynb):
            # The chains can be seen as of 'intermediary' size.
            # This means that they are long enough to be degraded (enzyme GP),
            # too long to be debranched (GDE), and too short to be branched (GBE).
            # Substrate for: GS, GP
            chain.substrate_status = [
                Status.GS_SUBSTRATE, Status.GP_SUBSTRATE]  # was '2'
        else:
            # status = 3 :
            # The chains are "long".
            # This means that they are long enough to be branched (GBE), # fixed GDE->GBE
            # but also to be elongated and degraded (GS and GP),
            # and again, too long to be debranched (GDE). Substrate for: GS, GP, GBE.
            chain.substrate_status = [
                Status.GS_SUBSTRATE, Status.GP_SUBSTRATE, Status.GBE_SUBSTRATE]  # was '3'

        return NotImplemented

    def scan_excluded_volume(self, position_for_candidate, chain_id) -> bool:
        logger.debug(
            f"radius of glucose sphere: {self.radius_of_glucose_sphere}")
        x, y, z = position_for_candidate
        logger.debug(f"x**2 + y**2 + z**2: {x**2 + y**2 + z**2}")
        if x**2 + y**2 + z**2 < self.radius_of_gn_core**2:
            logger.debug(f"was < {self.radius_of_gn_core}")
            return False

        for chain in self.chains:
            for (pos_x, pos_y, pos_z) in chain.glucose_positions:

                if (x-pos_x)**2 + (y-pos_y)**2 + (z-pos_z)**2 < self.radius_of_glucose_sphere**2:
                    # if it's not the same chain
                    if chain_id != chain.id:
                        return False
        return True

    def act_gs(self, chain_id=None):
        """Elongation"""
        if chain_id is None:
            chain_id = random.choice(self.find_chains_for_gs())
            logger.debug(f"Randomly selected chain {chain_id} for act_gs")
        c = self.get_chain_by_id(chain_id)
        if Status.GS_SUBSTRATE not in c.substrate_status:
            raise Exception("Cannot use chain {chain_id} for GS")

        logger.debug(f"calling angle3Dchain on {c.glucose_positions}")
        theta, phi = angle3Dchain(c.glucose_positions)

        # last glucose position
        x_nre, y_nre, z_nre = c.glucose_positions[-1]
        x = x_nre+GlycogenStructure.DISTANCE_BETWEEN_GLUCOSE_UNITS * \
            math.cos(
                phi)*math.cos(theta)                   # new x monomer position
        y = y_nre+GlycogenStructure.DISTANCE_BETWEEN_GLUCOSE_UNITS * \
            math.cos(
                phi)*math.sin(theta)                   # new y monomer position
        z = z_nre+GlycogenStructure.DISTANCE_BETWEEN_GLUCOSE_UNITS * \
            math.sin(phi)
        logger.debug(f"x,y,z: {x,y,z}")
        new_glucose_position_candidate = [x, y, z]
        if self.scan_excluded_volume(new_glucose_position_candidate, chain_id) is True:
            c.glucose_positions.append(new_glucose_position_candidate)
            self.update_chain_status(c.id)
            return True
        else:
            # return False(?)
            raise Exception(
                f"Could not attach {new_glucose_position_candidate } to chain {chain_id}")

    def act_gbe_flexible_model(self, chain_id: int = None, force_cleave_index: int | None = None, force_branch_index: int | None = None, force_theta: float | None = None, force_phi: float | None = None):
        """Branching using 'flexible model'.
        Use force_branch_index with an id to try and forch a branching point, similarly for force_cleave_index"""
        eps = 5.4  # TODO: can that be made a constant?
        available_for_gbe = self.find_chains_for_gbe()
        if not available_for_gbe:
            raise Exception(
                f"Could not find a chain suitable for act_gbe_flexible_model().")

        # get a random suitable chain when no id is provided
        if not chain_id:
            chain_id = random.choice(available_for_gbe)
            logger.debug(
                f"Randomly selected chain {chain_id}  during act_gbe_flexible_model().")
        if not chain_id in available_for_gbe:
            raise Exception(
                f"Chain {chain_id} is not suitable for act_gbe_flexible_model().")

        chain = self.get_chain_by_id(chain_id)
        # selection of the position to cleave:
        N = len(chain.glucose_positions)

        # TODO: refactor, this is used multiple times (see also: update_chain_status())
        if len(chain.daughters_positions) == 0:
            ind = 0
        else:
            ind = chain.daughters_positions[-1]+1

        start = ind + self.l_gbe_spacing + self.l_gbe_leftover
        end = N-self.l_gbe_transferred
        sub_list_to_cleave = [index for index in range(start, end)]

        index_where_to_cleave = random.choice(sub_list_to_cleave)
        logger.debug(
            f"Randomly chose index_where_to_cleave as  {index_where_to_cleave} on chain {chain_id}")
        if force_cleave_index is not None:
            index_where_to_cleave = force_cleave_index
            logger.debug(
                f"Force overwriting index_where_to_cleave as {force_cleave_index} on chain {chain_id}: {index_where_to_cleave}")
            if index_where_to_cleave not in sub_list_to_cleave:
                raise Exception(
                    f"Cannot use {index_where_to_cleave}, available: {sub_list_to_cleave}")

        N_new = N-index_where_to_cleave-1
        start = ind+self.l_gbe_spacing
        end = N-N_new-self.l_gbe_leftover
        sub_list_to_branch = [index for index in range(start, end)]
        index_where_to_branch = random.choice(sub_list_to_branch)
        logger.debug(
            f"Randomly chose index_where_to_branch as  {index_where_to_branch} on chain {chain_id}")
        if force_branch_index is not None:
            index_where_to_branch = force_branch_index
            logger.debug(
                f"Force overwriting index_where_to_branch as {force_branch_index} on chain {chain_id}: {index_where_to_branch}")
            if index_where_to_branch not in sub_list_to_branch:
                raise Exception(
                    f"Cannot use {index_where_to_branch}, available: {sub_list_to_branch}")

        chain_length = N-(index_where_to_cleave+1)
        alpha0, beta0 = angle3Dchain(chain.glucose_positions)
        theta = random.uniform(-math.pi, math.pi)
        logger.debug(
            f"Random theta {theta} (chain: {chain_id})")
        if force_theta is not None:
            theta = force_theta
            logger.debug(
                f"Force overwriting theta as {force_theta} on chain {chain_id}")
            if theta < -math.pi or theta > math.pi:
                raise Exception(
                    f"theta needs to be within bounds {-math.pi,math.pi}, was: {theta}")

        phi = random.uniform(-math.pi/2, math.pi/2)
        logger.debug(
            f"Random phi {phi} (chain: {chain_id})")
        if force_phi is not None:
            phi = force_phi
            logger.debug(
                f"Force overwriting phi as {force_phi} on chain {chain_id}")
            if phi < -math.pi/2 or phi > math.pi/2:
                raise Exception(
                    f"phi needs to be within bounds {-math.pi/2,math.pi/2}, was: {phi}")

        x0, y0, z0 = chain.glucose_positions[index_where_to_branch]
        l = GlycogenStructure.DISTANCE_BETWEEN_GLUCOSE_UNITS  # self.distance_between_units
        # position on x axis of the first monomer on the new chain
        xstart = x0+(l+eps)*math.cos(phi+beta0)*math.cos(theta+alpha0)
        ystart = y0+(l+eps)*math.cos(phi+beta0)*math.sin(theta+alpha0)
        zstart = z0+(l+eps)*math.sin(phi+beta0)

        # position on x axis of the last monomer on the new chain
        xend = xstart+(chain_length-1)*l*math.cos(phi+beta0) * \
            math.cos(theta+alpha0)
        # position on y axis of the last monomer on the new chain
        yend = ystart+(chain_length-1)*l*math.cos(phi+beta0) * \
            math.sin(theta+alpha0)
        zend = zstart+(chain_length-1)*l*math.sin(phi+beta0)

        X = np.linspace(xstart, xend, chain_length)
        Y = np.linspace(ystart, yend, chain_length)
        Z = np.linspace(zstart, zend, chain_length)

        cond = 0
        for k in range(chain_length):
            # TODO: is it correct that this doesn't get a chain id ???
            if self.scan_excluded_volume([X[k], Y[k], Z[k]], chain_id=None) == False:
                cond = 1
                # why not just return False (?)
        if cond == 0:
            # creation of the new chain:
            index = self.get_last_chain_index()
            new_index = index + 1
            self.create_a_chain(new_index)

            new_chain = self.get_chain_by_id(new_index)
            new_chain.generation = chain.generation+1

            for k in range(chain_length):
                new_chain.glucose_positions.append([X[k], Y[k], Z[k]])
                logger.debug(
                    f"Appending {[X[k],Y[k],Z[k]]} on new chain with index {new_index}")

            new_chain.mother_id = chain.id  # current chain is mother of new chain

            chain.daughters_ids.append(new_index)
            chain.daughters_positions.append(index_where_to_branch)
            # remove everything after index_where_to_cleave from current chain:
            # logger.debug(
            #    f"should remove glucose positions after index_where_to_cleave = {index_where_to_cleave}: before {chain.glucose_positions}")
            # ????
            chain.glucose_positions = chain.glucose_positions[:index_where_to_cleave+1]
            # logger.debug(
            #    f"should remove glucose positions after index_where_to_cleave = {index_where_to_cleave}: after {chain.glucose_positions}")

            self.update_chain_status(new_chain.id)
            self.update_chain_status(chain.id)
            return True
        else:
            return False

    def act_gp(self, chain_id: int = None) -> bool:
        if chain_id is None:
            chain_id = random.choice(self.find_chains_for_gp())
            logger.debug(f"Randomly selected chain {chain_id} for act_gp")
        c = self.get_chain_by_id(chain_id)
        if Status.GP_SUBSTRATE not in c.substrate_status:
            raise Exception("Cannot use chain {chain_id} for GP")
        else:
            if len(c.glucose_positions) > 0:
                logger.debug(
                    f"act_gp: {c.id} has {len(c.glucose_positions)}, will have {len(c.glucose_positions[:-1])}")
                c.glucose_positions = c.glucose_positions[:-1]
                logger.debug(
                    f"Update status of chain {chain_id} after act_gp old status was: {c.substrate_status}")
                self.update_chain_status(c.id)
                logger.debug(
                    f"Update status of chain {chain_id} after act_gp new status is : {c.substrate_status}")
                return True
            else:
                return False

    def act_gde(self, chain_id: int | None = None) -> bool:
        logger.debug(f"Chains for act_gde {self.find_chains_for_gde}")

        if not self.find_chains_for_gde():
            # TODO: make sure act_[gp,gs,gbe] behave the same
            logging.debug(
                f"Returning False as there were no chains for act_gde")
            return False
        if chain_id is None:
            chain_id = random.choice(self.find_chains_for_gde())
            logger.debug(f"Randomly selected chain {chain_id} for act_gde")

        c = self.get_chain_by_id(chain_id)
        if c.mother_id is None:
            raise Exception(
                f"Error during act_gde: Cannot de-branch chain {c.id} as it has no mother chain.")
        if Status.GDE_SUBSTRATE not in c.substrate_status:
            raise Exception(f"Cannot use chain {chain_id} for GDE")
        N = len(c.glucose_positions)

        for _ in range(N-1):
            self.act_gs(c.mother_id)

        mother = self.get_chain_by_id(c.mother_id)
        logger.debug(f"Mother was{c.mother_id}")

        mother.daughters_ids = mother.daughters_ids[:-1]
        mother.daughters_positions = mother.daughters_positions[:-1]
        logger.debug(f"Number of chains before {len(self.chains)}")
        self.chains = [chain for chain in self.chains if chain.id != c.id]
        logger.debug(f"Number of chains after {len(self.chains)}")

        return True

    def __repr__(self):
        suitable_chains = self.get_all_possible_chains_for_reactions()
        gbe_c = len(suitable_chains["gbe"])
        gs_c = len(suitable_chains["gs"])
        gde_c = len(suitable_chains["gde"])
        gp_c = len(suitable_chains["gp"])
        s = f"<GlycogenStructure@{id(self)} with {self.get_num_chains()} chains, {self.get_num_glucose_fixed()} glucose positions, {gbe_c} chains for GBE, {gs_c} chains for GS, {gde_c} chains for GDE, {gp_c} chains for GP >"

        return s

    def export_to_file(self, filename: Path | str, format: str = 'json'):
        available_formats = ['json']
        params_dct = self.get_parameters_as_dct()
        chains_dct = self.get_chains_as_dct()
        result_dct = params_dct | chains_dct  # merge
        if not format in available_formats:
            raise Exception(
                f"Cannot export to {format}. Available formats are {available_formats} .")
        if format == 'json':
            with open(Path(filename), 'w') as f:
                json.dump(result_dct, f, indent=2)

    def show_enzymes(self):
        header = f"{'Glycogen Branching Enzyme (GBE)':^30} {'Glycogen Synthase (GS)':^30} {'Glycogen Debranching Enzyme (GDE)':^30} {'Glycogen Phosphorylase (GP)':^30}"
        res = f"{header}\n{self.gbe:^30} {self.gs:^30} {self.gde:^30} {self.gp:^30}"
        return res

    def get_all_possible_chains_for_reactions(self) -> dict[str, int]:
        """Return mapping of [gbe, gs, gde, gp] to the number of chains that are suitable for them respecively"""
        gbe_c = self.find_chains_for_gbe()
        gs_c = self.find_chains_for_gs()
        gde_c = self.find_chains_for_gde()
        gp_c = self.find_chains_for_gp()
        return {"gbe": gbe_c, "gs": gs_c, "gde": gde_c, "gp": gp_c}

    def show_possible_reactions(self):
        suitable_chains = self.get_all_possible_chains_for_reactions()
        gbe_c = len(suitable_chains["gbe"])
        gs_c = len(suitable_chains["gs"])
        gde_c = len(suitable_chains["gde"])
        gp_c = len(suitable_chains["gp"])

        header = f"{'#chains for GBE':^30} {'#chains for GS':^30} {'#chains for GDE':^30} {'#chains for GP':^30}"
        res = f"{header}\n{gbe_c:^30} {gs_c:^30} {gde_c:^30} {gp_c:^30}"
        return res

    def get_parameters_as_dct(self):
        import inspect
        result = dict()
        for k in inspect.getfullargspec(GlycogenStructure.__init__).args:
            if k == 'self':
                continue
            if isinstance(self.__getattribute__(k), Model):
                result[k] = str(self.__getattribute__(k))
            else:
                result[k] = self.__getattribute__(k)
        return result

    def get_chains_as_dct(self):
        serializable = []
        for c in self.chains:
            serializable.append(c._to_serializable())
        return {'chains': serializable}
