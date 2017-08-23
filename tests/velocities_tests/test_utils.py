from alcor.models import Star
from alcor.models.velocities import (LepineCaseUVCloud,
                                     LepineCaseUWCloud,
                                     LepineCaseVWCloud)
from alcor.services.velocities.utils import star_cloud


def test_star_cloud(star: Star) -> None:
    cloud = star_cloud(star)

    assert isinstance(cloud, (LepineCaseUVCloud,
                              LepineCaseUWCloud,
                              LepineCaseVWCloud))
    assert cloud.group_id == star.group_id
