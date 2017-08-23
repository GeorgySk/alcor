from alcor.models import Star
from alcor.models.velocities import (LepineCaseUVCloud,
                                     LepineCaseUWCloud,
                                     LepineCaseVWCloud)
from alcor.services.velocities.utils import lepine_cloud


def test_lepine_cloud(star: Star) -> None:
    cloud = lepine_cloud(star)

    assert isinstance(cloud, (LepineCaseUVCloud,
                              LepineCaseUWCloud,
                              LepineCaseVWCloud))
    assert cloud.group_id == star.group_id
