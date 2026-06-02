from pathlib import Path

import pandas as pd

from scpi_pkg.scdataMulti import scdataMulti


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
