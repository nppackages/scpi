from pathlib import Path

import numpy as np
import pandas as pd

from scpi_pkg.scdata import scdata
from scpi_pkg.scdataMulti import scdataMulti
from scpi_pkg.scest import scest
from scpi_pkg.scpi import scpi


def test_scdata_multi_time_effect_uses_treatment_date_label():
    data = pd.read_csv(Path(__file__).resolve().parents[2] / "scpi_germany.csv")
    data["status"] = 0
    data.loc[(data["country"] == "West Germany") & (data["year"] >= 1991), "status"] = 1
    data.loc[(data["country"] == "Italy") & (data["year"] >= 1992), "status"] = 1

    covs_adj = {
        "Italy": ["constant", "trend"],
        "West Germany": [["constant", "trend"], ["constant", "trend"]],
    }

    result = scdataMulti(
        df=data,
        id_var="country",
        treatment_var="status",
        outcome_var="gdp",
        time_var="year",
        features={
            "Italy": ["gdp", "trade"],
            "West Germany": ["gdp", "infrate"],
        },
        constant={"Italy": True, "West Germany": False},
        cointegrated_data=True,
        cov_adj=covs_adj,
        effect="time",
        verbose=False,
    )

    assert result.effect == "time"
    assert result.A.shape == (125, 1)
    assert result.Y_post.shape == (25, 1)


def test_scpi_unit_effect_feature_designs_stay_conformable():
    data = pd.read_csv(Path(__file__).resolve().parents[2] / "scpi_germany.csv")
    data["status"] = 0
    data.loc[(data["country"] == "West Germany") & (data["year"] >= 1991), "status"] = 1
    data.loc[(data["country"] == "Italy") & (data["year"] >= 1992), "status"] = 1

    result = scdataMulti(
        df=data,
        id_var="country",
        treatment_var="status",
        outcome_var="gdp",
        time_var="year",
        features={"features": ["gdp", "trade"]},
        cov_adj={"cov_adj": [["constant"], ["trend"]]},
        constant=False,
        cointegrated_data=True,
        effect="unit",
        verbose=False,
    )

    np.random.seed(8894)
    intervals = scpi(result, sims=10, cores=1, verbose=False)

    assert len(intervals.CI_in_sample) == len(result.P)


def test_scpi_reuses_scest_output_single_unit():
    data = pd.read_csv(Path(__file__).resolve().parents[2] / "scpi_germany.csv")
    donor_units = sorted(set(data["country"]) - {"West Germany"})

    result = scdata(
        df=data,
        id_var="country",
        time_var="year",
        outcome_var="gdp",
        period_pre=np.arange(1960, 1991),
        period_post=np.arange(1991, 1998),
        unit_tr="West Germany",
        unit_co=donor_units,
        constant=True,
        cointegrated_data=True,
        verbose=False,
    )
    estimate = scest(result, w_constr={"name": "simplex"})
    w_bounds = np.zeros((len(estimate.Y_post_fit), 2))

    direct = scpi(
        result,
        w_constr={"name": "simplex"},
        w_bounds=w_bounds,
        sims=10,
        cores=1,
        e_method="gaussian",
        verbose=False,
    )
    reused = scpi(
        estimate,
        w_bounds=w_bounds,
        sims=10,
        cores=1,
        e_method="gaussian",
        verbose=False,
    )

    np.testing.assert_allclose(reused.b.to_numpy(), direct.b.to_numpy(), atol=1e-8, rtol=1e-8)
    np.testing.assert_allclose(reused.CI_in_sample.to_numpy(), direct.CI_in_sample.to_numpy(), atol=1e-8, rtol=1e-8)
    np.testing.assert_allclose(reused.CI_all_gaussian.to_numpy(), direct.CI_all_gaussian.to_numpy(), atol=1e-8, rtol=1e-8)


def test_scpi_reuses_scest_output_multi_unit():
    data = pd.read_csv(Path(__file__).resolve().parents[2] / "scpi_germany.csv")
    data["status"] = 0
    data.loc[(data["country"] == "West Germany") & (data["year"] >= 1991), "status"] = 1
    data.loc[(data["country"] == "Italy") & (data["year"] >= 1992), "status"] = 1

    result = scdataMulti(
        df=data,
        id_var="country",
        treatment_var="status",
        outcome_var="gdp",
        time_var="year",
        features={"features": ["gdp", "trade"]},
        cov_adj={"cov_adj": [["constant"], ["trend"]]},
        constant=False,
        cointegrated_data=True,
        effect="unit",
        verbose=False,
    )
    estimate = scest(result, w_constr={"name": "simplex"})
    w_bounds = np.zeros((len(estimate.Y_post_fit), 2))

    direct = scpi(
        result,
        w_constr={"name": "simplex"},
        w_bounds=w_bounds,
        sims=10,
        cores=1,
        e_method="gaussian",
        verbose=False,
    )
    reused = scpi(
        estimate,
        w_bounds=w_bounds,
        sims=10,
        cores=1,
        e_method="gaussian",
        verbose=False,
    )

    np.testing.assert_allclose(reused.b.to_numpy(), direct.b.to_numpy(), atol=1e-8, rtol=1e-8)
    np.testing.assert_allclose(reused.CI_in_sample.to_numpy(), direct.CI_in_sample.to_numpy(), atol=1e-8, rtol=1e-8)
    np.testing.assert_allclose(reused.CI_all_gaussian.to_numpy(), direct.CI_all_gaussian.to_numpy(), atol=1e-8, rtol=1e-8)


def test_scpi_reuse_matches_direct_simulated_inference_with_aligned_rng():
    data = pd.read_csv(Path(__file__).resolve().parents[2] / "scpi_germany.csv")
    data["status"] = 0
    data.loc[(data["country"] == "West Germany") & (data["year"] >= 1991), "status"] = 1
    data.loc[(data["country"] == "Italy") & (data["year"] >= 1992), "status"] = 1

    result = scdataMulti(
        df=data,
        id_var="country",
        treatment_var="status",
        outcome_var="gdp",
        time_var="year",
        features={"features": ["gdp", "trade"]},
        cov_adj={"cov_adj": [["constant"], ["trend"]]},
        constant=False,
        cointegrated_data=True,
        effect="unit",
        verbose=False,
    )

    np.random.seed(8894)
    direct = scpi(
        result,
        w_constr={"name": "simplex"},
        sims=10,
        cores=1,
        e_method="gaussian",
        verbose=False,
    )

    np.random.seed(8894)
    estimate = scest(result, w_constr={"name": "simplex"})
    reused = scpi(
        estimate,
        sims=10,
        cores=1,
        e_method="gaussian",
        verbose=False,
    )

    np.testing.assert_allclose(reused.b.to_numpy(), direct.b.to_numpy(), atol=1e-8, rtol=1e-8)
    np.testing.assert_allclose(reused.CI_in_sample.to_numpy(), direct.CI_in_sample.to_numpy(), atol=1e-8, rtol=1e-8)
    np.testing.assert_allclose(reused.CI_all_gaussian.to_numpy(), direct.CI_all_gaussian.to_numpy(), atol=1e-8, rtol=1e-8)
