import pytest
from jurassic_web.configuration import import_form, make_control, normalize, parse_ctl
from jurassic_web.catalog import datasets


def test_import_only_safe_keys():
    data=import_form('NG = 1\nEMITTER[0] = CO2\nND = 1\nNU[0] = 667.7820\nTBLBASE = /tmp/unsafe\nWRITE_BBT = 1\n')
    assert data['gases']==['CO2']
    assert data['channels']==['667.7820']
    assert 'tblbase' not in data
    assert data['output']=='bbt'


def test_generated_control_uses_absolute_lookup_path(app):
    ds=datasets(app)['airs']
    params=normalize({'geometry':'nadir','gases':['CO2'],'channels':['667.7820'], 'output':'bbt','rayds':5,'raydz':.1,'refrac':True,'kernel':False},ds)
    ctl=make_control(app.config['JURASSIC_ROOT'], params, ds)
    assert 'NG = 1' in ctl
    assert 'EMITTER[0] = CO2' in ctl
    assert 'RAYDS = 5' in ctl
    assert 'TBLBASE = /' in ctl
    assert 'WRITE_BBT = 1' in ctl


def test_bad_parameters_rejected(app):
    ds=datasets(app)['airs']
    base={'geometry':'nadir','gases':['CO2'],'channels':['667.7820'],'output':'bbt', 'rayds':10}
    with pytest.raises(ValueError): normalize({**base,'rayds':0},ds)
    with pytest.raises(ValueError): normalize({**base,'channels':['700.0000']},ds)
    with pytest.raises(ValueError): normalize({**base,'gases':['H2O']},ds)
    with pytest.raises(ValueError): normalize({**base,'output':'system'},ds)
