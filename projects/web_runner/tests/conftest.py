from __future__ import annotations
import os
import textwrap
from pathlib import Path
import pytest
from jurassic_web import create_app


@pytest.fixture
def fake_checkout(tmp_path):
    root = tmp_path / 'jurassic'
    src = root / 'src'
    src.mkdir(parents=True)
    data = root / 'tests' / 'data'
    data.mkdir(parents=True)
    (data / 'airs_667.7820_CO2.tab').write_text('# lookup table\n1 2 3 4\n')
    (data / 'airs_667.7820.filt').write_text('667 1\n')
    example = root / 'projects' / 'examples' / 'nadir'
    example.mkdir(parents=True)
    (example / 'nadir.ctl').write_text('TBLBASE = ../../../tests/data/airs\nNG = 1\nEMITTER[0] = CO2\nND = 1\nNU[0] = 667.7820\nWRITE_BBT = 1\n')
    (src / 'jurassic.c').write_text('/* fake */')
    code = textwrap.dedent('''\
        #!/usr/bin/env python3
        from pathlib import Path
        import os, sys
        name = Path(sys.argv[0]).name
        args = sys.argv[1:]
        if name == 'climatology':
            Path(args[1]).write_text('# $1 = time\\n# $2 = altitude [km]\\n# $3 = lon\\n# $4 = lat\\n# $5 = pressure\\n# $6 = temp\\n# $7 = CO2\\n# $8 = extinction\\n0 0 0 0 1017 285 0.0004 0\\n')
        elif name in ('nadir','limb','zenith'):
            Path(args[1]).write_text('0 700 0 0 0 0 -8 0 0 0 0 0\\n')
        elif name == 'formod':
            cols = ''.join(f'# ${i} = col {i}\\n' for i in range(1,13))
            Path(args[3]).write_text(cols+'0 700 0 0 0 0 -8 0 0 0 251 0.52\\n')
        elif name == 'kernel':
            Path(args[3]).write_text('# kernel\\n0 667.7820 0 0 0 0 0 TEMPERATURE 0 0 0 0 0.1\\n')
        print(name,'finished')
    ''')
    for exe in ('climatology','nadir','limb','zenith','formod','kernel'):
        path=src/exe
        path.write_text(code)
        path.chmod(0o755)
    return root


@pytest.fixture
def app(fake_checkout, tmp_path):
    return create_app({'TESTING': True, 'JURASSIC_ROOT': fake_checkout,
                       'INSTANCE_DIR':tmp_path/'instance', 'RUN_TIMEOUT':10})


@pytest.fixture
def client(app):
    return app.test_client()


@pytest.fixture
def headers(client):
    client.get('/')
    with client.session_transaction() as sess:
        token = sess['csrf']
    return {'X-CSRF-Token': token}
