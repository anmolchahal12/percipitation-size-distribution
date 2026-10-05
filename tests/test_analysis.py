"""Numerical and input-validation checks, using synthetic measurements."""
import tempfile
import unittest
from pathlib import Path
import numpy as np
import pandas as pd
from ppt_size_distribution import empirical_cdf, load_data, summarize


class AnalysisTests(unittest.TestCase):
    def frame(self):
        return pd.DataFrame({"aging_time_min": np.repeat([15, 60, 165], 3),
                             "diameter_nm": [10, 20, 30, 20, 30, 40, 30, 40, 50]})

    def test_sample_statistics(self):
        stats = summarize(self.frame())
        np.testing.assert_allclose(stats.mean_nm, [20, 30, 40])
        np.testing.assert_allclose(stats.sample_std_nm, [10, 10, 10])
        np.testing.assert_allclose(stats.median_nm, [20, 30, 40])
        self.assertEqual(stats.n.tolist(), [3, 3, 3])

    def test_tied_ecdf(self):
        x, y = empirical_cdf(np.array([20, 10, 20, 30]))
        np.testing.assert_array_equal(x, [10, 20, 30])
        np.testing.assert_allclose(y, [.25, .75, 1])

    def test_invalid_input(self):
        for value in [0, -10, np.nan, np.inf, "bad"]:
            with self.subTest(value=value), tempfile.TemporaryDirectory() as tmp:
                frame = self.frame().astype(object)
                frame.loc[0, "diameter_nm"] = value
                path = Path(tmp) / "data.csv"
                frame.to_csv(path, index=False)
                with self.assertRaises(ValueError):
                    load_data(path)

    def test_missing_condition(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "data.csv"
            self.frame().iloc[:6].to_csv(path, index=False)
            with self.assertRaises(ValueError):
                load_data(path)


if __name__ == "__main__":
    unittest.main()
