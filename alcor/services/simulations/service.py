import logging
import random
import uuid
from collections import namedtuple
from functools import partial
from types import GeneratorType
from typing import (Dict,
                    Tuple,
                    Callable)

import numpy as np
import pandas as pd
import pydevd
from sqlalchemy.orm.session import Session

from alcor.models import (STAR_PARAMETERS_NAMES,
                          Group,
                          Star)
from alcor.models.simulation import Parameter
from alcor.models.star import GalacticDiskType
from alcor.services.simulations.coordinates import set_coordinates
from alcor.services.simulations.equatorial_coordinates import \
    assign_equatorial_coordinates
from alcor.services.simulations.luminosities import get_white_dwarfs
from alcor.services.simulations.magnitudes import assign_estimated_values
from alcor.services.simulations.polar import assign_polar_coordinates
from alcor.services.simulations.proper_motions import assign_proper_motions
from alcor.services.simulations.sphere_stars import generate_stars
from alcor.services.simulations.tracks import (read_cooling,
                                               read_table)
from alcor.services.simulations.velocities import (set_velocities,
                                                   halo_stars_velocities)
from alcor.types import (GridParametersInfoType,
                         CSVParametersInfoType,
                         UnitRangeGeneratorType,
                         GaussianGeneratorType)
from alcor.utils import validate_header
from . import grid

logger = logging.getLogger(__name__)

# TODO: where to put this?
VelocityVector = namedtuple('VelocityVector', ['u', 'v', 'w'])


def run(*,
        geometry: str,
        precision: int,
        grid_parameters_info: GridParametersInfoType,
        csv_parameters_info: CSVParametersInfoType,
        session: Session) -> None:
    for parameters_values in grid.parameters_values(
            parameters_info=grid_parameters_info,
            precision=precision):
        group = Group(id=uuid.uuid4())
        white_dwarfs = run_simulation(parameters_values=parameters_values,
                                      geometry=geometry)

        header = white_dwarfs.columns
        validate_header(header,
                        possible_columns_names=STAR_PARAMETERS_NAMES)

        parameters = [Parameter(group_id=group.id,
                                name=name,
                                value=str(value))
                      for name, value in parameters_values.items()]

        white_dwarfs['group_id'] = group.id
        white_dwarfs.to_sql(name=Star.__tablename__,
                            con=session.get_bind(),
                            if_exists='append',
                            index=False)

        session.add(group)
        session.add_all(parameters)
        session.commit()


def run_simulation(
        *,
        parameters_values: Dict[str, float],
        geometry: str,
        da_tracks_path: str = 'input_data/da_cooling.hdf5',
        db_tracks_path: str = 'input_data/db_cooling.hdf5',
        one_table_path: str = 'input_data/one_wds_tracks.hdf5',
        da_colors_path: str = 'input_data/da_colors.hdf5',
        db_colors_path: str = 'input_data/db_colors.hdf5',
        da_metallicities: Tuple[int] = (1, 10, 30, 60),
        db_metallicities: Tuple[int] = (1, 10, 60),
        interest_parameters: Tuple[str, ...] = ('cooling_time',
                                                'effective_temperature',
                                                'luminosity'),
        colors: Tuple[str, ...] = ('color_u',
                                   'color_b',
                                   'color_v',
                                   'color_r',
                                   'color_i',
                                   'color_j'),
        time_bins_count: int = 5000,
        max_stars_count: int = 6000000,
        burst_formation_factor: float = 5.,
        formation_rate_parameter: float = 25.,
        generator: Callable[[float, float], float] = random.uniform,
        array_generator: GeneratorType = np.random.uniform,
        unit_range_generator: UnitRangeGeneratorType = np.random.rand,
        signs_generator: Callable[[int], np.ndarray] = partial(
                np.random.choice, [-1, 1]),
        gaussian_generator: GaussianGeneratorType = np.random.normal,
        min_mass: float = 0.04,
        max_mass: float = 50.,
        mass_relation_parameter: float = 1.,
        chandrasekhar_limit: float = 1.4,
        solar_metallicity: float = 0.01,
        subsolar_metallicity: float = 0.001,
        scale_length: float = 3.5,
        solar_galactocentric_distance: float = 8.5,
        thin_disk_scale_height: float = 0.25,
        thick_disk_scale_height: float = 0.9,
        halo_core_radius: float = 5.,
        thin_disk_velocity_std: VelocityVector = VelocityVector(
                u=32.4, v=23., w=18.1),
        thick_disk_velocity_std: VelocityVector = VelocityVector(
                u=50., v=56., w=34.),
        peculiar_solar_velocity: VelocityVector = VelocityVector(
                u=-11., v=-12., w=-7.),
        lsr_velocity: float = -220.,
        oort_constant_a: float = 14.4,
        oort_constant_b: float = -12.8,
        max_carbon_oxygen_core_wd_mass: float = 1.14,
        db_to_da_fraction: float = 0.2) -> pd.DataFrame:
    if geometry == 'cones':
        raise NotImplementedError('Only spherical geometry available now')

    # TODO: take this out of the function and out of loop
    da_cooling_tracks = read_cooling(path=da_tracks_path,
                                     metallicities=da_metallicities,
                                     interest_parameters=interest_parameters)
    db_cooling_tracks = read_cooling(path=db_tracks_path,
                                     metallicities=db_metallicities,
                                     interest_parameters=interest_parameters)
    one_table = read_table(path=one_table_path,
                           interest_parameters=interest_parameters + colors)
    da_colors = read_table(path=da_colors_path,
                           interest_parameters=colors + ('luminosity',))
    db_colors = read_table(path=db_colors_path,
                           interest_parameters=colors + ('luminosity',))

    main_sequence_stars = generate_stars(
            max_stars_count=max_stars_count,
            time_bins_count=time_bins_count,
            thin_disk_age=parameters_values['thin_disk_age'],
            thick_disk_age=parameters_values['thick_disk_age'],
            halo_age=parameters_values['halo_age'],
            halo_stars_formation_time=parameters_values[
                'halo_stars_formation_time'],
            burst_age=parameters_values['burst_time'],
            thick_disk_formation_rate_exponent=parameters_values[
                'thick_disk_star_formation_exponent'],
            thick_disk_stars_fraction=parameters_values[
                'thick_disk_stars_fraction'],
            halo_stars_fraction=parameters_values['halo_stars_fraction'],
            initial_mass_function_exponent=parameters_values[
                'initial_mass_function_exponent'],
            sector_radius_kpc=parameters_values['radius'],
            burst_formation_factor=burst_formation_factor,
            formation_rate_parameter=formation_rate_parameter,
            mass_reduction_factor=parameters_values['mass_reduction_factor'],
            generator=generator,
            min_mass=min_mass,
            max_mass=max_mass)

    # TODO: this is also calculated in generate_stars
    max_age = max(parameters_values['thin_disk_age'],
                  parameters_values['thick_disk_age'],
                  parameters_values['halo_age'])

    white_dwarfs = get_white_dwarfs(
            stars=main_sequence_stars,
            max_galactic_structure_age=max_age,
            mass_relation_parameter=mass_relation_parameter,
            chandrasekhar_limit=chandrasekhar_limit,
            max_mass=max_mass,
            solar_metallicity=solar_metallicity,
            subsolar_metallicity=subsolar_metallicity)

    white_dwarfs = assign_polar_coordinates(
            stars=white_dwarfs,
            sector_radius=parameters_values['radius'],
            scale_length=scale_length,
            solar_galactocentric_distance=solar_galactocentric_distance,
            thin_disk_scale_height=thin_disk_scale_height,
            thick_disk_scale_height=thick_disk_scale_height,
            halo_core_radius=halo_core_radius,
            generator=array_generator,
            unit_range_generator=unit_range_generator,
            signs_generator=signs_generator)

    set_velocities(stars=white_dwarfs,
                   thin_disk_velocity_std=thin_disk_velocity_std,
                   thick_disk_velocity_std=thick_disk_velocity_std,
                   peculiar_solar_velocity=peculiar_solar_velocity,
                   solar_galactocentric_distance=solar_galactocentric_distance,
                   oort_constant_a=oort_constant_a,
                   oort_constant_b=oort_constant_b,
                   generator=gaussian_generator)

    set_coordinates(
            stars=white_dwarfs,
            solar_galactocentric_distance=solar_galactocentric_distance)

    halo_stars_mask = (white_dwarfs['galactic_disk_type']
                       == GalacticDiskType.halo)

    (white_dwarfs.loc[halo_stars_mask, 'u_velocity'],
     white_dwarfs.loc[halo_stars_mask, 'v_velocity'],
     white_dwarfs.loc[halo_stars_mask, 'w_velocity']) = halo_stars_velocities(
            galactic_longitudes=(
                white_dwarfs.loc[halo_stars_mask, 'galactic_longitude']
                .values),
            thetas_cylindrical=(
                white_dwarfs.loc[halo_stars_mask, 'theta_cylindrical'].values),
            peculiar_solar_velocity=peculiar_solar_velocity,
            lsr_velocity=lsr_velocity,
            spherical_velocity_component_sigma=lsr_velocity / np.sqrt(2.),
            generator=gaussian_generator)

    sin_longitude = np.sin(white_dwarfs['galactic_longitude'])
    cos_longitude = np.cos(white_dwarfs['galactic_longitude'])
    sin_latitude = np.sin(white_dwarfs['galactic_latitude'])
    cos_latitude = np.cos(white_dwarfs['galactic_latitude'])

    assign_proper_motions(white_dwarfs,
                          sin_longitude=sin_longitude,
                          cos_longitude=cos_longitude,
                          sin_latitude=sin_latitude,
                          cos_latitude=cos_latitude)

    assign_equatorial_coordinates(white_dwarfs,
                                  sin_latitude=sin_latitude,
                                  cos_latitude=cos_latitude)

    # pydevd.settrace('dockerhost', port=20111)
    return assign_estimated_values(
            white_dwarfs,
            max_carbon_oxygen_core_wd_mass=max_carbon_oxygen_core_wd_mass,
            db_to_da_fraction=db_to_da_fraction,
            da_cooling_sequences=da_cooling_tracks,
            da_color_table=da_colors,
            db_cooling_sequences=db_cooling_tracks,
            db_color_table=db_colors,
            one_color_table=one_table)
