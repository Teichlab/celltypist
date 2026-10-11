import unittest

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy import sparse

from celltypist.samples import downsample_adata


class DownsampleTests(unittest.TestCase):
    def test_each_preserves_seeded_cell_order(self):
        # Labels are interleaved and deliberately not in sorted order.
        labels = list('bcabcabcaabc')
        for dtype in ('object', 'category'):
            for sparse_input in (False, True):
                with self.subTest(dtype=dtype, sparse_input=sparse_input):
                    values = np.arange(60, dtype=np.float32).reshape(12, 5)
                    adata = AnnData(
                        sparse.csr_matrix(values) if sparse_input else values,
                        obs=pd.DataFrame({'label': pd.Series(labels, dtype=dtype)}),
                    )
                    kwargs = dict(mode='each', by='label', random_state=3)
                    expected = np.array([9, 5, 3, 6, 11, 7])
                    np.testing.assert_array_equal(
                        downsample_adata(adata, n_cells=2, **kwargs), expected,
                    )
                    selected = downsample_adata(adata, n_cells=2, return_index=False, **kwargs)
                    np.testing.assert_array_equal(selected.obs_names, adata.obs_names[expected])
                    np.testing.assert_array_equal(
                        selected.X.toarray() if sparse_input else selected.X,
                        values[expected],
                    )
                    self.assertFalse(selected.is_view)
                    np.testing.assert_array_equal(
                        downsample_adata(adata, n_cells=20, **kwargs),
                        [9, 5, 2, 8, 3, 6, 0, 10, 11, 7, 4, 1],
                    )

    def test_each_does_not_sample_missing_labels(self):
        adata = AnnData(np.zeros((6, 1)), obs=pd.DataFrame({'label': [2, np.nan, 1, 2, np.nan, 1]}))
        np.testing.assert_array_equal(
            downsample_adata(adata, mode='each', n_cells=10, by='label', random_state=3),
            [5, 2, 3, 0],
        )
        # NaT can leave repeated values in NumPy's sorted object output.
        labels = pd.to_datetime(['2020-01-02', '2020-01-01', None, '2020-01-02'], utc=True)
        adata = AnnData(np.zeros((4, 1)), obs=pd.DataFrame({'label': labels}))
        np.random.seed(0)
        expected = np.concatenate([
            np.random.choice(np.where(adata.obs['label'] == label)[0],
                             min(10, np.sum(adata.obs['label'] == label)), replace=False)
            for label in np.unique(adata.obs['label'])
        ])
        expected_state = np.random.get_state()
        np.testing.assert_array_equal(
            downsample_adata(adata, mode='each', n_cells=10, by='label', random_state=0), expected,
        )
        actual_state = np.random.get_state()
        np.testing.assert_array_equal(actual_state[1], expected_state[1])
        self.assertEqual(actual_state[2:], expected_state[2:])


if __name__ == '__main__':
    unittest.main()
