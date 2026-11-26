import io
import time
import zipfile


def spec():
    return {'geometry':'nadir','dataset':'airs','gases':['CO2'],'channels':['667.7820'],'output':'bbt',
            'rayds':10,'raydz':.1,'refrac':True,'kernel':False}


def test_dashboard_and_system(client):
    assert client.get('/').status_code==200
    system=client.get('/api/system').json
    assert system['ready']
    assert client.get('/api/catalog').json['airs']['channels']['667.7820']['gases']['CO2']['data']


def test_csrf_rejected(client):
    r=client.post('/api/config/preview',json=spec())
    assert r.status_code==403


def test_preview_and_run(client,headers,app):
    p=client.post('/api/config/preview',json=spec(),headers=headers)
    assert p.status_code==200
    assert 'NG = 1' in p.json['text']
    res=client.post('/api/runs',json=spec(),headers=headers)
    assert res.status_code==202
    id=res.json['id']
    record=None
    for _ in range(60):
        record=client.get('/api/runs/'+id).json
        if record['status'] in ('completed','failed'):break
        time.sleep(0.15)
    assert record['status']=='completed',record
    assert client.get('/api/runs/'+id+'/chart/radiance').json['traces'][0]['y']==[251]
    assert client.get('/api/runs/'+id+'/table').json['rows'][0][-2:]==[251,0.52]
    assert client.get('/api/runs/'+id+'/download/rad.tab').status_code==200
    assert client.get('/api/runs/'+id+'/download/rad.csv').status_code==200
    assert client.get('/api/runs/'+id+'/download/all_results.zip').status_code==200
    assert client.get('/api/runs/'+id+'/download/../../etc/passwd').status_code in (404,308)


def test_import_configuration(client,headers):
    upload={'configuration':(io.BytesIO(b'NG = 1\nEMITTER[0] = CO2\nND = 1\nNU[0] = 667.7820\n'), 'test.ctl')}
    r=client.post('/api/config/import',data=upload,headers=headers)
    assert r.status_code==200
    assert r.json['values']['gases']==['CO2']


def test_register_dataset(client,headers):
    archive=io.BytesIO()
    with zipfile.ZipFile(archive,'w') as z:
        z.writestr('custom_900.0000_CO2.tab', '# table\n1 2 3 4\n')
        z.writestr('custom_900.0000.filt','900 1\n')
    archive.seek(0)
    r=client.post('/api/catalog/upload',data={'label':'Lab Spectroscopy','dataset':(archive,'data.zip')},headers=headers)
    assert r.status_code==201, r.json
    cat=client.get('/api/catalog').json
    assert r.json['id'] in cat
