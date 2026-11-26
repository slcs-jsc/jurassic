"""Regression tests: completed runs must serve actual images for all plot tabs."""
from __future__ import annotations
import time
from pathlib import Path


def _new_run(client, headers, kernel=True):
    spec = {'geometry': 'nadir', 'dataset': 'airs', 'gases': ['CO2'],
            'channels': ['667.7820'], 'output': 'bbt', 'rayds': 10,
            'raydz': 0.1, 'refrac': True, 'kernel': kernel}
    result = client.post('/api/runs', json=spec, headers=headers)
    assert result.status_code == 202, result.json
    uid = result.json['id']
    for _ in range(120):
        run = client.get('/api/runs/' + uid).json
        if run['status'] in ('completed', 'failed'):
            assert run['status'] == 'completed', run
            return uid
        time.sleep(0.1)
    raise AssertionError('Timed out waiting for simulation')


def test_completed_run_serves_all_plot_tabs(client, headers, app):
    uid = _new_run(client, headers)
    for metric in ('radiance', 'transmittance', 'atmosphere', 'kernels'):
        for theme in ('light', 'dark'):
            response = client.get(f'/api/runs/{uid}/plot/{metric}.png?theme={theme}')
            assert response.status_code == 200, (metric,theme,response.get_data(as_text=True)[:200])
            assert response.mimetype == 'image/png'
            assert response.data.startswith(b'\x89PNG\r\n\x1a\n')
    assert client.get(f'/api/runs/{uid}/plot/atmosphere.png?download=1').headers.get('Content-Disposition','').startswith('attachment;')
    assert client.get(f'/api/runs/{uid}/plot/../secret.png').status_code in (404, 308)
    assert client.get(f'/api/runs/{uid}/plot/radiance.png?theme=unknown').status_code == 400


def test_legacy_run_regenerates_missing_images(client, headers, app):
    uid = _new_run(client, headers, kernel=False)
    folder=Path(app.config['INSTANCE_DIR'])/'runs'/uid/'plots'
    (folder/'atmosphere.png').unlink(missing_ok=True)
    result = client.get(f'/api/runs/{uid}/plot/atmosphere.png')
    assert result.status_code==200
    assert (folder/'atmosphere.png').is_file()
    assert client.get(f'/api/runs/{uid}/plot/kernels.png').status_code==404
