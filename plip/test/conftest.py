import pytest

from plip.basic import config


@pytest.fixture(autouse=True)
def restore_config():
    """Undo changes to the global plip.basic.config made by a test."""
    saved = {key: value for key, value in vars(config).items() if key.isupper()}
    yield
    for key, value in saved.items():
        setattr(config, key, value)
