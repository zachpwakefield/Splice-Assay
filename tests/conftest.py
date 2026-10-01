import matplotlib
import pytest

matplotlib.use("Agg")

import splice_assay as sa  # noqa: E402
from splice_assay import example  # noqa: E402


@pytest.fixture(scope="session")
def tables():
    return example.make()


@pytest.fixture(scope="session")
def ds(tables):
    return sa.Dataset.from_tables(**tables)


@pytest.fixture(scope="session")
def results(ds):
    return sa.analyse(ds)


@pytest.fixture(scope="session")
def gtf_path(tmp_path_factory):
    p = tmp_path_factory.mktemp("gtf") / "annotation.gtf"
    p.write_text(example.gtf_text())
    return p
